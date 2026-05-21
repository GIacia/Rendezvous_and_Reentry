clear; clc; close all;

fprintf('================================================\n');
fprintf(' Unmanned Rendezvous Mission 6-DOF Simulator\n');
fprintf('================================================\n\n');

%% 1. Initialization
sys = Mission_Config();

% Budget Tracking Table
Budget = table('Size',[0 4], 'VariableTypes',{'string','double','double','double'}, ...
    'VariableNames',{'Phase','DeltaV_ms','Fuel_Consumed_kg','Remaining_Mass_kg'});

m_current = sys.Chaser_Mass_Init;

% 체이서 초기 상태 
x_insert = [sys.Re+sys.h_insert; 0; 0];
% v_insert = sqrt(sys.mu / (sys.Re + sys.h_insert));
v_insert = sqrt(sys.mu / (sys.Re + sys.h_insert) * (1 - sys.J2 * (sys.Re / (sys.Re + sys.h_insert))^2 * (3*x_insert(3)^2/(sys.Re+sys.h_insert)^2 - 1)));  % we want the total mechanical energy to be same as no J2 condition.
X_chaser_init = [x_insert;  0; 0; v_insert;  0;0;0;1;  0;0;0;  m_current];
target_radius = sys.Re + sys.h_wait;

% 타겟 위성 상태 세팅 
% v_t0 = sqrt(sys.mu / (sys.Re + sys.h_target));
x_target = [0; 0; sys.Re+sys.h_target];
v_t0 = sqrt(sys.mu / (sys.Re + sys.h_target) * (1 - sys.J2 * (sys.Re / (sys.Re + sys.h_target))^2 * (3*x_target(3)^2/(sys.Re+sys.h_target)^2 - 1)));  % we want the total mechanical energy to be same as no J2 condition.
X_target = [x_target; -v_t0; 0; 0];

%% 2. Phase 1: Phasing & Homing 
fprintf('[Phase 1] 3-DOF 물리 시뮬레이션 시작 (J2 섭동 포함)...\n');

% 변수 세팅
custom_params = struct();
custom_params.TOF = 10000; % LAMBERT 모드용 예비값
custom_params.target_pos = [0; 0; sys.Re + sys.h_wait]; % LAMBERT 모드용 예비값

% Prameters for "CUSTOM IMPULSE" Mode: Phase_angle, delta_v, gamma value from external python program
custom_params.phase_angle = deg2rad(3.8242);        % [rad] 예: deg2rad(3.88)
custom_params.delta_v     = 57.9955;                % [m/s] 예: 57.2
custom_params.gamma       = deg2rad(0.076877);      % [rad] 예: deg2rad(-0.0153)
custom_params.phase_angle_unit = "rad"; % "rad" or "deg"
custom_params.gamma_unit       = "rad"; % "rad" or "deg"

% Phase 1 종료 목표: Target 기준 LVLH 상대위치 [m]
custom_params.desired_rel_lvlh = [0; -5000; 0];

% 이벤트 탐지/전파 설정
custom_params.time_step = 10;           % [s] 기본 propagation step
custom_params.dt_phase = 5;            % [s] phase-angle 탐지용 coarse step
custom_params.dt_capture = 5;          % [s] capture 탐지용 coarse step
custom_params.max_wait = 2*pi * sys.mu^(-0.5) / ((sys.Re + 300e3)^(-1.5) - (sys.Re + 500e3)^(-1.5));
custom_params.max_capture_time = 30000; % [s] impulse 이후 최대 추적 시간
custom_params.phase_tol = 1e-7;         % [rad] phase-angle event tolerance
custom_params.event_time_tol = 1e-3;    % [s] event refinement tolerance
custom_params.capture_pos_tol = 0.5;    % [m] desired_rel_lvlh 도달 판정 허용오차
custom_params.min_capture_time = 0;     % [s] 너무 이른 local minimum을 무시하고 싶으면 키우기

% 원하는 모드로 문자열만 변경하세요: "CUSTOM_IMPULSE", "HOHMANN", "LAMBERT"
phasing_mode = "CUSTOM_IMPULSE";

[X_chaser, dV_p1, fuel_p1, hist_p1, X_target] = Phasing_Propagator(sys, X_chaser_init, target_radius, phasing_mode, custom_params, X_target, 1);
m_current = X_chaser(14);

Budget = [Budget; {"Phase 1: "+phasing_mode, dV_p1, fuel_p1, m_current}];

h_vec_p1 = cross(X_target(1:3), X_target(4:6));
i_u_p1 = X_target(1:3)/norm(X_target(1:3));
k_u_p1 = h_vec_p1/norm(h_vec_p1);
j_u_p1 = cross(k_u_p1, i_u_p1);
C_I2L_p1 = [i_u_p1'; j_u_p1'; k_u_p1'];
rel_p1_lvlh = C_I2L_p1 * (X_chaser(1:3) - X_target(1:3));
miss_p1 = norm(rel_p1_lvlh - custom_params.desired_rel_lvlh);

fprintf('   Phase 1 종료 LVLH 상대위치: [%+.3f, %+.3f, %+.3f] m\n', rel_p1_lvlh(1), rel_p1_lvlh(2), rel_p1_lvlh(3));
fprintf('   desired_rel_lvlh 오차: %.6f m\n', miss_p1);
fprintf('   Phase 1 종료 고도: %.2f km\n', (norm(X_chaser(1:3)) - sys.Re)/1000);

%% 3. Phase 2: +R-bar Approach & Berthing (Waypoint-Impulse Simulation)
fprintf('\n[Phase 2] Starting waypoint-impulsive +R-bar approach...\n');

% -------------------------------------------------------------------------
% Phase 2 concept
%   - S2:  -V-bar 5 km hold point, inherited from Phase 1 when possible.
%   - S2 -> S3: CW/Hill-targeted free-flight arc to a +R-bar hold point.
%   - S3 -> S4: small impulsive hops down the +R-bar corridor.
%   - Every leg is executed as: compute LVLH delta-V -> apply impulse ->
%     propagate the full nonlinear Env_EOM with zero continuous thrust.
% -------------------------------------------------------------------------
p2 = struct();
p2.dt = 1.0;                         % [s] propagation step during proximity ops
p2.S2 = [0; -5000; 0];               % [m] -V-bar hold point after Phase 1
p2.S3 = [500; 0; 0];                 % [m] +R-bar hold point
p2.S4 = [30; 0; 0];                  % [m] final +R-bar approach point
p2.initial_S2_tol = 50.0;            % [m] if Phase 1 is farther, insert cleanup leg to S2
p2.arc_count = 6;                    % number of S2->S3 targeting legs
p2.tof_initial_s2 = 1200;            % [s] cleanup transfer time to S2, only used if needed
p2.tof_arc_total = 3600;             % [s] total time for S2->S3 arc
p2.rbar_hops = [430 370 310 255 205 160 120 85 55 30];  % [m] +R-bar hop targets
p2.tof_hop = 300;                    % [s] transfer time per +R-bar hop
p2.capture_pos_tol = 0.25;           % [m] terminal position tolerance at S4
p2.max_terminal_refines = 4;         % extra S4 targeting attempts if nonlinear drift remains
p2.tof_terminal_refine = 180;        % [s] time per final retargeting leg
p2.Isp_fallback_s = 220;             % [s] used only if Mission_Config has no Isp field

% Current actual relative state after Phase 1
[rel0_lvlh, vrel0_lvlh] = rel_state_lvlh_local(X_chaser, X_target);
fprintf('   Phase 2 initial LVLH rel-pos = [%+.2f, %+.2f, %+.2f] m\n', ...
        rel0_lvlh(1), rel0_lvlh(2), rel0_lvlh(3));
fprintf('   Phase 2 initial LVLH rel-vel = [%+.4f, %+.4f, %+.4f] m/s\n', ...
        vrel0_lvlh(1), vrel0_lvlh(2), vrel0_lvlh(3));

% Build S2 -> S3 cycloidal-like corridor waypoints in the R/V plane.
% The shape is only a reference corridor. Each leg itself is physically flown
% by an impulse followed by nonlinear propagation.
u_arc = linspace(0, 1, p2.arc_count + 1);
arc_waypoints = zeros(3, numel(u_arc));
for ii = 1:numel(u_arc)
    u = u_arc(ii);
    arc_waypoints(:,ii) = [p2.S3(1) * (1 - cos(0.5*pi*u)); ...
                           p2.S2(2) * cos(0.5*pi*u); ...
                           0];
end

phase2_targets = [];
phase2_tofs = [];
phase2_names = strings(1,0);

if norm(rel0_lvlh - p2.S2) > p2.initial_S2_tol
    phase2_targets = [phase2_targets, p2.S2];
    phase2_tofs = [phase2_tofs, p2.tof_initial_s2];
    phase2_names(end+1) = "cleanup_to_S2";
    fprintf('   Phase 1 endpoint is %.1f m away from S2, so a physical cleanup leg is inserted.\n', norm(rel0_lvlh - p2.S2));
end

for ii = 2:size(arc_waypoints,2)
    phase2_targets = [phase2_targets, arc_waypoints(:,ii)];
    phase2_tofs = [phase2_tofs, p2.tof_arc_total / p2.arc_count];
    phase2_names(end+1) = "S2_to_S3_arc_" + string(ii-1);
end

for ii = 1:numel(p2.rbar_hops)
    phase2_targets = [phase2_targets, [p2.rbar_hops(ii); 0; 0]];
    phase2_tofs = [phase2_tofs, p2.tof_hop];
    phase2_names(end+1) = "Rbar_hop_" + string(ii);
end

% Initialize Phase 2 history containers
hist_pos = rel0_lvlh;
hist_mass = X_chaser(14);
hist_p2 = struct();
hist_p2.pos = X_chaser(1:3);
hist_p2.time = 0;
hist_p2.mass = X_chaser(14);
phase2_time = 0;
dV_p2 = 0;
fuel_p2 = 0;

% Execute each waypoint leg
for seg = 1:size(phase2_targets,2)
    r_goal = phase2_targets(:,seg);
    tof_seg = phase2_tofs(seg);

    [r_rel, v_rel] = rel_state_lvlh_local(X_chaser, X_target);
    n_now = target_mean_motion_local(X_target);
    dv_lvlh = cw_delta_v_to_waypoint_local(r_rel, v_rel, r_goal, tof_seg, n_now);

    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;

    fprintf('   %-18s: target [%+7.1f,%+7.1f,%+6.1f] m, TOF %5.0f s, dV %8.4f m/s\n', ...
            char(phase2_names(seg)), r_goal(1), r_goal(2), r_goal(3), tof_seg, norm(dv_lvlh));

    [X_chaser, X_target, seg_hist] = propagate_pair_free_local(X_chaser, X_target, tof_seg, p2.dt, sys, phase2_time);
    phase2_time = seg_hist.time(end);

    hist_pos = [hist_pos, seg_hist.rel];
    hist_mass = [hist_mass, seg_hist.mass];
    hist_p2.pos = [hist_p2.pos, seg_hist.pos];
    hist_p2.time = [hist_p2.time, seg_hist.time];
    hist_p2.mass = [hist_p2.mass, seg_hist.mass];

    [r_arrive, v_arrive] = rel_state_lvlh_local(X_chaser, X_target);
    fprintf('      arrival error: pos %.3f m, rel-vel %.4f m/s\n', norm(r_arrive - r_goal), norm(v_arrive));
end

% Nonlinear/J2 residual cleanup at S4. This is still done by impulses and
% free propagation, not by overwriting the state.
for refine = 1:p2.max_terminal_refines
    [r_rel, v_rel] = rel_state_lvlh_local(X_chaser, X_target);
    pos_err = norm(r_rel - p2.S4);
    if pos_err <= p2.capture_pos_tol
        break;
    end

    n_now = target_mean_motion_local(X_target);
    dv_lvlh = cw_delta_v_to_waypoint_local(r_rel, v_rel, p2.S4, p2.tof_terminal_refine, n_now);
    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;

    fprintf('   terminal_refine_%d : remaining pos %.3f m, TOF %5.0f s, dV %8.4f m/s\n', ...
            refine, pos_err, p2.tof_terminal_refine, norm(dv_lvlh));

    [X_chaser, X_target, seg_hist] = propagate_pair_free_local(X_chaser, X_target, p2.tof_terminal_refine, p2.dt, sys, phase2_time);
    phase2_time = seg_hist.time(end);

    hist_pos = [hist_pos, seg_hist.rel];
    hist_mass = [hist_mass, seg_hist.mass];
    hist_p2.pos = [hist_p2.pos, seg_hist.pos];
    hist_p2.time = [hist_p2.time, seg_hist.time];
    hist_p2.mass = [hist_p2.mass, seg_hist.mass];
end

% Final braking impulse: cancel residual LVLH relative velocity at S4.
[r_final, v_final] = rel_state_lvlh_local(X_chaser, X_target);
dv_stop_lvlh = -v_final;
[X_chaser, fuel_stop] = apply_impulse_lvlh_local(X_chaser, X_target, dv_stop_lvlh, sys, p2);
dV_p2 = dV_p2 + norm(dv_stop_lvlh);
fuel_p2 = fuel_p2 + fuel_stop;
hist_mass(end) = X_chaser(14);
hist_p2.mass(end) = X_chaser(14);

[r_final, v_final] = rel_state_lvlh_local(X_chaser, X_target);
fprintf('   terminal braking dV: %.4f m/s\n', norm(dv_stop_lvlh));
fprintf('   Phase 2 final LVLH rel-pos = [%+.3f, %+.3f, %+.3f] m\n', r_final(1), r_final(2), r_final(3));
fprintf('   Phase 2 final LVLH rel-vel = [%+.5f, %+.5f, %+.5f] m/s\n', v_final(1), v_final(2), v_final(3));
fprintf('   Phase 2 total dV: %.4f m/s, fuel: %.4f kg, elapsed: %.2f min\n', dV_p2, fuel_p2, phase2_time/60);

m_current = X_chaser(14);
Budget = [Budget; {"Phase 2: Waypoint +R-bar", dV_p2, fuel_p2, m_current}];

%% 4. Phase 3: Re-entry (500 km -> 200 km)
fprintf('\n[Phase 3] 3-DOF 재진입 시뮬레이션 시작 (J2 섭동 및 노이즈 역추진)...\n');

% 고추력으로 대기권 인터페이스에 하강
reentry_mode = "HOHMANN"; 
fprintf('   선택된 기동 방식: %s\n', reentry_mode);

target_reentry_r = sys.Re + sys.h_reentry;

% Phase 2 Berthing 직후의 X_chaser 상태를 그대로 입력하여 하강
[X_chaser, dV_p3, fuel_p3, hist_p3, X_target] = Phasing_Propagator(sys, X_chaser, target_reentry_r, reentry_mode, custom_params, X_target, 3);
m_current = X_chaser(14);

Budget = [Budget; {"Phase 3: "+reentry_mode, dV_p3, fuel_p3, m_current}];
fprintf('   도달 고도: %.2f km (목표: 200.00 km)\n', (norm(X_chaser(1:3)) - sys.Re)/1000);

%% 5. Summary and Output
disp(' ');
disp('=== Final Mission Delta-V & Mass Budget ===');
disp(Budget);
fprintf('Total Delta-V required: %.2f m/s\n', sum(Budget.DeltaV_ms));
fprintf('Total Fuel consumed:    %.2f kg\n', sum(Budget.Fuel_Consumed_kg));
fprintf('Remaining Dry Mass:     %.2f kg\n', m_current);

% Plotting R-bar trajectory
figure('Name','Proximity Operations','Color','w');
plot(hist_pos(2,:), hist_pos(1,:), 'b-', 'LineWidth', 2); hold on;
plot(phase2_targets(2,:), phase2_targets(1,:), 'ko', 'MarkerSize', 4);
plot(0,0,'r^','MarkerSize',10,'MarkerFaceColor','r');
grid on; xlabel('V-bar (m)'); ylabel('R-bar (m)');
title('R-bar Approach Trajectory (LVLH)');
legend('Chaser Trajectory', 'Commanded Waypoints', 'Target');
set(gca, 'XDir', 'reverse'); % Flight direction to the left

%% 6. 종합 임무 프로파일 시각화 (Mission Profile Visualization)
fprintf('\n시각화 대시보드를 생성 중입니다...\n');

figure('Name', 'Mission Comprehensive Dashboard', 'Color', 'w', 'Position', [100, 100, 1400, 800]);

% --- (1) 3D Mission Trajectory (ECI Frame) ---
subplot(2, 2, [1, 3]);
% 지구 그리기
[XE, YE, ZE] = sphere(50);
surf(XE*sys.Re/1e3, YE*sys.Re/1e3, ZE*sys.Re/1e3, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on; grid on; axis equal;

% 궤적 그리기 (단위: km)
plot3(hist_p1.pos(1,:)/1e3, hist_p1.pos(2,:)/1e3, hist_p1.pos(3,:)/1e3, 'g-', 'LineWidth', 1.5); % Phase 1 (상승)
plot3(hist_p2.pos(1,:)/1e3, hist_p2.pos(2,:)/1e3, hist_p2.pos(3,:)/1e3, 'b-', 'LineWidth', 1.5); % Phase 2 (근접운용)
plot3(hist_p3.pos(1,:)/1e3, hist_p3.pos(2,:)/1e3, hist_p3.pos(3,:)/1e3, 'r-', 'LineWidth', 1.5); % Phase 3 (재진입)

% 타겟 궤도 (500km) 참조선 그리기
theta = linspace(0, 2*pi, 100);
plot3((sys.Re+sys.h_target)/1e3*cos(theta), zeros(1,100), (sys.Re+sys.h_target)/1e3*sin(theta), 'k--', 'LineWidth', 1);

title('3D Mission Trajectory (Earth Centered Inertial)');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
legend('Earth', 'Phase 1: Phasing (Ascent)', 'Phase 2: +R-bar Prox Ops', 'Phase 3: De-orbit (Descent)', 'Target Orbit (500km)', 'Location', 'best');
view(45, 30); % 보기 좋은 각도 설정

% --- (2) Altitude vs Time Profile ---
subplot(2, 2, 2);
alt_p1 = vecnorm(hist_p1.pos) - sys.Re; % 중심으로부터의 거리 - 지구반지름 = 고도
alt_p2 = vecnorm(hist_p2.pos) - sys.Re;
alt_p3 = vecnorm(hist_p3.pos) - sys.Re;
time_p1 = hist_p1.time / 3600; % 시간(초)을 시간(Hour)으로 변환
time_p2 = hist_p2.time / 3600;
time_p3 = hist_p3.time / 3600;

plot(time_p1, alt_p1/1e3, 'g-', 'LineWidth', 2); hold on;
plot(time_p2 + time_p1(end), alt_p2/1e3, 'b-', 'LineWidth', 2);
% 재진입 시간은 Phasing + Prox Ops 직후에 이어지도록 offset 설정
plot(time_p3 + time_p1(end) + time_p2(end), alt_p3/1e3, 'r-', 'LineWidth', 2); 
yline(sys.h_insert/1e3, 'k:', 'Insertion (300km)');
yline(sys.h_target/1e3, 'k--', 'Target (500km)');
yline(sys.h_reentry/1e3, 'k:', 'Re-entry (200km)');

title('Altitude Profile');
xlabel('Mission Time (Hours)'); ylabel('Altitude (km)');
legend('Phasing Maneuver', '+R-bar Prox Ops', 'Re-entry Maneuver');
grid on;

% --- (3) Mass Depletion Profile ---
subplot(2, 2, 4);
plot(time_p1, hist_p1.mass, 'g-', 'LineWidth', 2); hold on;
plot(time_p2 + time_p1(end), hist_p2.mass, 'b-', 'LineWidth', 2);
plot(time_p3 + time_p1(end) + time_p2(end), hist_p3.mass, 'r-', 'LineWidth', 2);

title('Spacecraft Mass Depletion (Fuel Consumption)');
xlabel('Mission Time (Hours)'); ylabel('Mass (kg)');
legend('Fuel used in Phasing', 'Fuel used in Prox Ops', 'Fuel used in Re-entry');
grid on;

%% Local helper functions for Phase 2 waypoint-impulse proximity operations
function [r_rel, v_rel, C_I2L] = rel_state_lvlh_local(X_chaser, X_target)
    r_c = X_chaser(1:3);
    v_c = X_chaser(4:6);
    r_t = X_target(1:3);
    v_t = X_target(4:6);

    h_vec = cross(r_t, v_t);
    i_u = r_t / norm(r_t);          % +R-bar, radial outward
    k_u = h_vec / norm(h_vec);      % +H-bar, orbit normal
    j_u = cross(k_u, i_u);          % +V-bar, along-track
    C_I2L = [i_u'; j_u'; k_u'];

    r_rel = C_I2L * (r_c - r_t);
    w_eci = h_vec / norm(r_t)^2;
    w_lvlh = C_I2L * w_eci;
    v_rel = C_I2L * (v_c - v_t) - cross(w_lvlh, r_rel);
end

function n = target_mean_motion_local(X_target)
    r_t = X_target(1:3);
    v_t = X_target(4:6);
    h_vec = cross(r_t, v_t);
    n = norm(h_vec) / norm(r_t)^2;
end

function dv_lvlh = cw_delta_v_to_waypoint_local(r0, v0, rf, tof, n)
    % CW/Hill single-impulse targeting: choose post-impulse v0+dv so that
    % the linearized relative trajectory reaches rf after tof.
    nt = n * tof;
    s = sin(nt);
    c = cos(nt);

    Phi_rr = [4 - 3*c,        0, 0; ...
              6*(s - nt),     1, 0; ...
              0,              0, c];

    Phi_rv = [s/n,              2*(1-c)/n,       0; ...
              -2*(1-c)/n,      (4*s - 3*nt)/n,  0; ...
              0,               0,               s/n];

    if rcond(Phi_rv) < 1e-10
        error('Phase 2 CW targeting became ill-conditioned. Change the segment TOF away from singular transfer times.');
    end

    v_req = Phi_rv \ (rf - Phi_rr*r0);
    dv_lvlh = v_req - v0;
end

function [X_chaser, fuel_used] = apply_impulse_lvlh_local(X_chaser, X_target, dv_lvlh, sys, p2)
    [~, ~, C_I2L] = rel_state_lvlh_local(X_chaser, X_target);
    dv_eci = C_I2L' * dv_lvlh;
    dv_mag = norm(dv_eci);

    X_chaser(4:6) = X_chaser(4:6) + dv_eci;

    m0 = X_chaser(14);
    Isp = impulse_isp_local(sys, p2);
    g0 = 9.80665;
    if Isp > 0 && dv_mag > 0
        m1 = m0 * exp(-dv_mag/(Isp*g0));
        fuel_used = m0 - m1;
        X_chaser(14) = m1;
    else
        fuel_used = 0;
    end
end

function Isp = impulse_isp_local(sys, p2)
    candidates = {'Isp_Impulsive', 'Isp_RCS', 'Isp_Thruster', 'Isp'};
    Isp = NaN;
    for ii = 1:numel(candidates)
        if isfield(sys, candidates{ii}) && isnumeric(sys.(candidates{ii})) && isscalar(sys.(candidates{ii}))
            Isp = sys.(candidates{ii});
            break;
        end
    end
    if isnan(Isp)
        Isp = p2.Isp_fallback_s;
    end
end

function [X_chaser, X_target, hist] = propagate_pair_free_local(X_chaser, X_target, tof, dt, sys, t0)
    n_steps = ceil(tof/dt);
    hist.rel = zeros(3, n_steps);
    hist.pos = zeros(3, n_steps);
    hist.mass = zeros(1, n_steps);
    hist.time = zeros(1, n_steps);

    t_elapsed = 0;
    for kk = 1:n_steps
        dt_step = min(dt, tof - t_elapsed);
        t_abs = t0 + t_elapsed;

        % Target propagation: 6-DOF wrapper with zero force/torque.
        X_t_state = [X_target; zeros(7,1); sys.Target_Mass];
        k1_t = Env_EOM(t_abs,             X_t_state,               [0;0;0], [0;0;0], sys, false);
        k2_t = Env_EOM(t_abs+dt_step/2,   X_t_state+k1_t*dt_step/2,[0;0;0], [0;0;0], sys, false);
        k3_t = Env_EOM(t_abs+dt_step/2,   X_t_state+k2_t*dt_step/2,[0;0;0], [0;0;0], sys, false);
        k4_t = Env_EOM(t_abs+dt_step,     X_t_state+k3_t*dt_step,  [0;0;0], [0;0;0], sys, false);
        X_target = X_target + (dt_step/6)*(k1_t(1:6) + 2*k2_t(1:6) + 2*k3_t(1:6) + k4_t(1:6));

        % Chaser free-flight propagation after the impulse.
        k1 = Env_EOM(t_abs,             X_chaser,             [0;0;0], [0;0;0], sys, true);
        k2 = Env_EOM(t_abs+dt_step/2,   X_chaser+k1*dt_step/2,[0;0;0], [0;0;0], sys, true);
        k3 = Env_EOM(t_abs+dt_step/2,   X_chaser+k2*dt_step/2,[0;0;0], [0;0;0], sys, true);
        k4 = Env_EOM(t_abs+dt_step,     X_chaser+k3*dt_step,  [0;0;0], [0;0;0], sys, true);
        X_chaser = X_chaser + (dt_step/6)*(k1 + 2*k2 + 2*k3 + k4);

        if norm(X_chaser(7:10)) > 0
            X_chaser(7:10) = X_chaser(7:10) / norm(X_chaser(7:10));
        end

        t_elapsed = t_elapsed + dt_step;
        [r_rel, ~] = rel_state_lvlh_local(X_chaser, X_target);
        hist.rel(:,kk) = r_rel;
        hist.pos(:,kk) = X_chaser(1:3);
        hist.mass(kk) = X_chaser(14);
        hist.time(kk) = t0 + t_elapsed;
    end
end
