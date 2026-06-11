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
%   - S2 -> S3: one retrograde/-V-bar impulse, then free cycloidal drift.
%       * v0 is selected from |delta_R|max = 4*v0/n.
%       * S3 is not pre-fixed. S3 is the first point where V-bar coordinate
%         crosses zero after the cycloidal drift starts.
%   - S3 -> S4: short two-impulse R-bar hops. Each hop uses CW targeting for
%     the departure impulse, nonlinear Env_EOM propagation, then a small
%     braking impulse to hold before the next hop.
%   - No state overwriting: all corrections are impulses + free propagation.
% -------------------------------------------------------------------------
p2 = struct();
p2.dt = 1.0;                         % [s] propagation step during proximity ops
p2.S2 = [0; -5000; 0];               % [m] -V-bar hold point after Phase 1
p2.S3 = [NaN; NaN; NaN];             % [m] determined by first V-bar crossing
p2.S4_R_abs = 30;                    % [m] final stand-off distance on the same R-bar side as S3
p2.S4 = [NaN; NaN; NaN];             % [m] determined after S3 is detected
p2.initial_S2_tol = 50.0;            % [m] if Phase 1 is farther, insert cleanup leg to S2
p2.tof_initial_s2 = 1200;            % [s] cleanup transfer time to S2, only used if needed
p2.delta_R_cycloid = 400;            % [m] maximum radial excursion magnitude for the cycloid
p2.vbar_burn_sign = -1;              % -1 means apply the S2 impulse in the -V-bar direction
p2.vbar_cross_tol = 1.0;             % [m] acceptable error for detecting first V-bar crossing
p2.max_cycloid_orbits = 4.0;         % safety limit for S2->S3 free drift
p2.rbar_hop_count = 8;               % number of S3->S4 R-bar hops
p2.tof_hop = 300;                    % [s] transfer time per R-bar hop
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

% Keep a list of meaningful Phase 2 points for plotting.
phase2_targets = [];
phase2_target_times = [];      % [s] elapsed time from Phase 2 start
phase2_tofs = [];
phase2_names = strings(1,0);

% -------------------------------------------------------------------------
% 2-0. Optional cleanup to S2, then brake to make S2 a true hold point.
% -------------------------------------------------------------------------
if norm(rel0_lvlh - p2.S2) > p2.initial_S2_tol
    tof_seg = p2.tof_initial_s2;
    [r_rel, v_rel] = rel_state_lvlh_local(X_chaser, X_target);
    n_now = target_mean_motion_local(X_target);
    dv_lvlh = cw_delta_v_to_waypoint_local(r_rel, v_rel, p2.S2, tof_seg, n_now);

    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;

    fprintf('   cleanup_to_S2     : target [%+7.1f,%+7.1f,%+6.1f] m, TOF %5.0f s, dV %8.4f m/s\n', ...
            p2.S2(1), p2.S2(2), p2.S2(3), tof_seg, norm(dv_lvlh));

    [X_chaser, X_target, seg_hist] = propagate_pair_free_local(X_chaser, X_target, tof_seg, p2.dt, sys, phase2_time);
    phase2_time = seg_hist.time(end);

    hist_pos = [hist_pos, seg_hist.rel];
    hist_mass = [hist_mass, seg_hist.mass];
    hist_p2.pos = [hist_p2.pos, seg_hist.pos];
    hist_p2.time = [hist_p2.time, seg_hist.time];
    hist_p2.mass = [hist_p2.mass, seg_hist.mass];
end

% S2 hold trim: cancel any remaining LVLH relative velocity before the
% cycloidal free drift. This is necessary because the 4*v0/n relation assumes
% a clean S2 starting condition except for the deliberate -V-bar impulse.
[r_s2, v_s2] = rel_state_lvlh_local(X_chaser, X_target);
dv_hold_s2_lvlh = -v_s2;
if norm(dv_hold_s2_lvlh) > 1e-6
    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_hold_s2_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_hold_s2_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;
    hist_mass(end) = X_chaser(14);
    hist_p2.mass(end) = X_chaser(14);
    fprintf('   S2_hold_trim      : residual rel-vel canceled, dV %8.4f m/s\n', norm(dv_hold_s2_lvlh));
end

% Record S2 waypoint time after optional cleanup and hold trim.
phase2_targets = [phase2_targets, p2.S2];
phase2_target_times = [phase2_target_times, phase2_time];
phase2_names(end+1) = "S2";

% -------------------------------------------------------------------------
% 2-1. S2 -> S3 natural cycloidal drift using one -V-bar impulse.
% -------------------------------------------------------------------------
n_now = target_mean_motion_local(X_target);
v0_cycloid = n_now * p2.delta_R_cycloid / 4;
dv_cycloid_lvlh = [0; p2.vbar_burn_sign * v0_cycloid; 0];
expected_delta_R_code = 4 * dv_cycloid_lvlh(2) / n_now;

[X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_cycloid_lvlh, sys, p2);
dV_p2 = dV_p2 + norm(dv_cycloid_lvlh);
fuel_p2 = fuel_p2 + fuel_seg;
hist_mass(end) = X_chaser(14);
hist_p2.mass(end) = X_chaser(14);

fprintf('   S2_to_S3_cycloid  : single %+.4f m/s V-bar impulse, |delta_R|max %.1f m', ...
        dv_cycloid_lvlh(2), p2.delta_R_cycloid);
fprintf(' (current LVLH signed delta_R %.1f m)\n', expected_delta_R_code);

% Propagate freely until the chaser first reaches the target R-bar line,
% i.e., until V-bar coordinate crosses zero.
max_cycloid_time = p2.max_cycloid_orbits * 2*pi / n_now;
cycloid_elapsed = 0;
[r_prev, ~] = rel_state_lvlh_local(X_chaser, X_target);
vbar_prev = r_prev(2);
initial_vbar_sign = sign(vbar_prev);
if initial_vbar_sign == 0
    initial_vbar_sign = -1;
end
crossed_rbar = false;

while cycloid_elapsed < max_cycloid_time
    [X_chaser, X_target, seg_hist] = propagate_pair_free_local(X_chaser, X_target, p2.dt, p2.dt, sys, phase2_time);
    phase2_time = seg_hist.time(end);
    cycloid_elapsed = cycloid_elapsed + p2.dt;

    hist_pos = [hist_pos, seg_hist.rel];
    hist_mass = [hist_mass, seg_hist.mass];
    hist_p2.pos = [hist_p2.pos, seg_hist.pos];
    hist_p2.time = [hist_p2.time, seg_hist.time];
    hist_p2.mass = [hist_p2.mass, seg_hist.mass];

    [r_now, v_now] = rel_state_lvlh_local(X_chaser, X_target);
    vbar_now = r_now(2);

    if abs(vbar_now) <= p2.vbar_cross_tol || sign(vbar_now) ~= initial_vbar_sign
        crossed_rbar = true;
        break;
    end
    vbar_prev = vbar_now;
end

if ~crossed_rbar
    error('Phase 2 cycloidal drift did not reach the R-bar line within %.2f orbits. Check S2 offset, delta_R_cycloid, and sign convention.', p2.max_cycloid_orbits);
end

[r_s3, v_s3] = rel_state_lvlh_local(X_chaser, X_target);
p2.S3 = [r_s3(1); 0; 0];
phase2_targets = [phase2_targets, p2.S3];
phase2_target_times = [phase2_target_times, phase2_time];
phase2_names(end+1) = "S3";
phase2_tofs = [phase2_tofs, cycloid_elapsed];

fprintf('      S3 detected at first V-bar crossing after %.2f min\n', cycloid_elapsed/60);
fprintf('      S3 actual LVLH rel-pos = [%+.3f, %+.3f, %+.3f] m\n', r_s3(1), r_s3(2), r_s3(3));
fprintf('      S3 actual LVLH rel-vel = [%+.5f, %+.5f, %+.5f] m/s\n', v_s3(1), v_s3(2), v_s3(3));

% Brake at S3 before starting the controlled R-bar hop sequence. Without this,
% the next S3->S4 hop starts with a large natural-drift velocity and can become
% an unsafe fly-by rather than an approach.
dv_hold_s3_lvlh = -v_s3;
if norm(dv_hold_s3_lvlh) > 1e-6
    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_hold_s3_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_hold_s3_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;
    hist_mass(end) = X_chaser(14);
    hist_p2.mass(end) = X_chaser(14);
    fprintf('   S3_hold_brake     : cycloid rel-vel canceled, dV %8.4f m/s\n', norm(dv_hold_s3_lvlh));
end

% -------------------------------------------------------------------------
% 2-2. S3 -> S4 R-bar approach.
% -------------------------------------------------------------------------
% Use the same signed R-bar side as the actual S3. This avoids commanding the
% chaser to cross through the target just because of a sign-convention mismatch.
approach_R_sign = sign(p2.S3(1));
if approach_R_sign == 0
    approach_R_sign = sign(expected_delta_R_code);
end
if approach_R_sign == 0
    approach_R_sign = 1;
end
p2.S4 = [approach_R_sign * p2.S4_R_abs; 0; 0];

if abs(p2.S3(1)) <= abs(p2.S4(1))
    warning('S3 radial distance %.2f m is already inside or near S4 %.2f m. Skipping R-bar hops and using terminal refine only.', p2.S3(1), p2.S4(1));
    rbar_hops = p2.S4(1);
else
    rbar_hops = linspace(p2.S3(1), p2.S4(1), p2.rbar_hop_count + 1);
    rbar_hops = rbar_hops(2:end);
end

for ii = 1:numel(rbar_hops)
    r_goal = [rbar_hops(ii); 0; 0];
    tof_seg = p2.tof_hop;

    [r_rel, v_rel] = rel_state_lvlh_local(X_chaser, X_target);
    n_now = target_mean_motion_local(X_target);
    dv_lvlh = cw_delta_v_to_waypoint_local(r_rel, v_rel, r_goal, tof_seg, n_now);

    [X_chaser, fuel_seg] = apply_impulse_lvlh_local(X_chaser, X_target, dv_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_lvlh);
    fuel_p2 = fuel_p2 + fuel_seg;

    fprintf('   Rbar_hop_%02d      : target [%+7.1f,%+7.1f,%+6.1f] m, TOF %5.0f s, depart dV %8.4f m/s\n', ...
            ii, r_goal(1), r_goal(2), r_goal(3), tof_seg, norm(dv_lvlh));

    [X_chaser, X_target, seg_hist] = propagate_pair_free_local(X_chaser, X_target, tof_seg, p2.dt, sys, phase2_time);
    phase2_time = seg_hist.time(end);

    hist_pos = [hist_pos, seg_hist.rel];
    hist_mass = [hist_mass, seg_hist.mass];
    hist_p2.pos = [hist_p2.pos, seg_hist.pos];
    hist_p2.time = [hist_p2.time, seg_hist.time];
    hist_p2.mass = [hist_p2.mass, seg_hist.mass];

    [r_arrive, v_arrive] = rel_state_lvlh_local(X_chaser, X_target);

    % Make each hop a real hold point by braking the residual relative velocity.
    % This is more physical/safe than only doing one final brake at S4.
    dv_brake_lvlh = -v_arrive;
    [X_chaser, fuel_brake] = apply_impulse_lvlh_local(X_chaser, X_target, dv_brake_lvlh, sys, p2);
    dV_p2 = dV_p2 + norm(dv_brake_lvlh);
    fuel_p2 = fuel_p2 + fuel_brake;
    hist_mass(end) = X_chaser(14);
    hist_p2.mass(end) = X_chaser(14);

    phase2_targets = [phase2_targets, r_goal];
    phase2_target_times = [phase2_target_times, phase2_time];
    phase2_tofs = [phase2_tofs, tof_seg];

    if ii == numel(rbar_hops)
        phase2_names(end+1) = "S4";
    else
        phase2_names(end+1) = "WP" + string(ii);
    end

    [r_hold, v_hold] = rel_state_lvlh_local(X_chaser, X_target);
    fprintf('      arrival error: pos %.3f m, brake dV %.4f m/s, post-brake rel-vel %.5f m/s\n', ...
            norm(r_hold - r_goal), norm(dv_brake_lvlh), norm(v_hold));
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
[~, v_final] = rel_state_lvlh_local(X_chaser, X_target);
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
Budget = [Budget; {"Phase 2: Cycloid + R-bar", dV_p2, fuel_p2, m_current}];

%% 4. Phase 3: Re-entry
fprintf('\n[Phase 3] 3-DOF re-entry simulation start...\n');

% Select Phase 3 mode here.
%   "HOHMANN"        : legacy direct FPA-targeted Hohmann-style descent.
%   "R_BAR_200_FPA" : lower to 200 km, wait until below the target on R-bar,
%                     then inject toward the 120 km / 4 deg FPA interface.
reentry_mode = "R_BAR_200_FPA";
reentry_mode_env = getenv('RENDEZVOUS_PHASE3_MODE');
if strlength(reentry_mode_env) > 0
    reentry_mode = string(reentry_mode_env);
end
fprintf('   selected Phase 3 mode: %s\n', reentry_mode);

% Set true to charge propellant for the final 200 km -> 120 km injection in
% R_BAR_200_FPA mode. The default keeps this injection fuel-excluded.
charge_final_reentry_fuel = false;
charge_final_reentry_fuel_env = getenv('RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL');
if strlength(charge_final_reentry_fuel_env) > 0
    charge_final_reentry_fuel = parse_bool_setting_local(charge_final_reentry_fuel_env);
end
custom_params.charge_final_reentry_fuel = charge_final_reentry_fuel;
fprintf('   charge final 200 km -> 120 km injection fuel: %s\n', bool_label_local(charge_final_reentry_fuel));

reentry_info = struct();

if reentry_mode == "HOHMANN"
    [target_reentry_r, ~] = compute_reentry_target_radius_local(norm(X_chaser(1:3)), sys);

    % Legacy behavior: feed the FPA-derived terminal radius into the phasing
    % propagator and let it perform the Hohmann-style descent.
    [X_chaser, dV_p3, fuel_p3, hist_p3, X_target] = ...
        Phasing_Propagator(sys, X_chaser, target_reentry_r, "HOHMANN", custom_params, X_target, 3);

    fprintf('   final altitude: %.2f km\n', (norm(X_chaser(1:3)) - sys.Re)/1000);

elseif reentry_mode == "R_BAR_200_FPA"
    [X_chaser, X_target, dV_p3, fuel_p3, hist_p3, reentry_info] = ...
        run_phase3_rbar_200_fpa_local(sys, X_chaser, X_target, custom_params);

    fprintf('   final altitude: %.2f km\n', (norm(X_chaser(1:3)) - sys.Re)/1000);
    if reentry_info.final_injection_fuel_charged
        fprintf('   final 200 km -> 120 km injection dV %.4f m/s and fuel %.4f kg are included in Phase 3 budget.\n', ...
                reentry_info.final_injection_dV, reentry_info.final_injection_fuel);
    else
        fprintf('   final 200 km -> 120 km injection dV %.4f m/s is included in Phase 3 dV, but its fuel is excluded.\n', ...
                reentry_info.final_injection_dV);
    end

else
    error('Unknown Phase 3 reentry_mode: %s. Use "HOHMANN" or "R_BAR_200_FPA".', reentry_mode);
end

m_current = X_chaser(14);
Budget = [Budget; {"Phase 3: "+reentry_mode, dV_p3, fuel_p3, m_current}];

legacy_phase3_block_enabled = strcmp(getenv('RENDEZVOUS_ENABLE_LEGACY_PHASE3_BLOCK'), '1');
if legacy_phase3_block_enabled
%% Legacy Phase 3 block retained for reference.
fprintf('\n[Phase 3] 3-DOF 재진입 시뮬레이션 시작 (J2 섭동 및 노이즈 역추진)...\n');

% 고추력으로 대기권 인터페이스에 하강
reentry_mode = "HOHMANN"; 
fprintf('   선택된 기동 방식: %s\n', reentry_mode);

% target_reentry_r = sys.Re + sys.h_reentry;

% 120 km 고도에서 원하는 flight_path_angle을 얻기 위한 파라미터 계산
a_tmp = norm(X_chaser(1:3));
c_tmp = sys.Re + 120000;
gamma_tmp = sys.reentry_flight_path_angle;
phi = 2 * atan(a_tmp / (a_tmp - c_tmp) * tan(gamma_tmp));
e_tmp = sin(gamma_tmp) / (sin(phi) * cos(gamma_tmp) - cos(phi) * sin(gamma_tmp));
b_tmp = (1 - e_tmp) / (1 + e_tmp) * a_tmp;
target_reentry_r = b_tmp;

% Phase 2 Berthing 직후의 X_chaser 상태를 그대로 입력하여 하강
[X_chaser, dV_p3, fuel_p3, hist_p3, X_target] = Phasing_Propagator(sys, X_chaser, target_reentry_r, reentry_mode, custom_params, X_target, 3);
m_current = X_chaser(14);

Budget = [Budget; {"Phase 3: "+reentry_mode, dV_p3, fuel_p3, m_current}];
fprintf('   도달 고도: %.2f km (목표: 200.00 km)\n', (norm(X_chaser(1:3)) - sys.Re)/1000);

% 120 km에 가장 가까운 history 지점에서 flight path angle 계산
end
target_alt_fpa = 120e3;  % [m]

alt_hist_p3 = vecnorm(hist_p3.pos, 2, 1) - sys.Re;
[~, idx_120] = min(abs(alt_hist_p3 - target_alt_fpa));

r_120 = hist_p3.pos(:, idx_120);
v_120 = hist_p3.vel(:, idx_120);

rhat_120 = r_120 / norm(r_120);
v_radial_120 = dot(v_120, rhat_120);                         % +면 상승, -면 하강
v_horizontal_120 = norm(v_120 - v_radial_120 * rhat_120);     % local horizontal 성분

fpa_120_deg = rad2deg(atan2(v_radial_120, v_horizontal_120)); % signed FPA
below_horizon_120_deg = -fpa_120_deg;                         % 하강이면 양수

fprintf('   120 km 근접 지점: 실제 고도 %.2f km, FPA %.3f deg, 수평선 아래 각도 %.3f deg\n', ...
        alt_hist_p3(idx_120)/1000, fpa_120_deg, below_horizon_120_deg);

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
plot(phase2_targets(2,:), phase2_targets(1,:), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w');
plot(0,0,'r^','MarkerSize',10,'MarkerFaceColor','r');

% Waypoint time labels
for kk = 1:size(phase2_targets, 2)
    x_wp = phase2_targets(2, kk);   % V-bar
    y_wp = phase2_targets(1, kk);   % R-bar

    label_txt = sprintf('%s, %.1f min', ...
        char(phase2_names(kk)), phase2_target_times(kk)/60);

    text(x_wp, y_wp, ['  ' label_txt], ...
        'FontSize', 8, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
end

grid on;
xlabel('V-bar (m)');
ylabel('R-bar (m)');
title('R-bar Approach Trajectory (LVLH)');
legend('Chaser Trajectory', 'Commanded Waypoints', 'Target', 'Location', 'best');
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

%% Local helper functions for Phase 3 selectable re-entry modes
function [X_chaser, X_target, dV_used, fuel_used, hist, info] = run_phase3_rbar_200_fpa_local(sys, X_chaser, X_target, custom_params)
    info = struct();
    hist = init_phase3_hist_local();
    dV_used = 0;
    fuel_used = 0;

    phase3_params = custom_params;
    phase3_params.capture_pos_tol = inf;

    target_200_r = sys.Re + sys.h_reentry;
    fprintf('   R_BAR_200_FPA leg 1: lowering to %.1f km circular parking orbit...\n', sys.h_reentry/1000);
    [X_chaser, dV_drop, fuel_drop, hist_drop, X_target] = ...
        Phasing_Propagator(sys, X_chaser, target_200_r, "HOHMANN", phase3_params, X_target, 3);

    hist = append_phase3_hist_local(hist, hist_drop);
    dV_used = dV_used + dV_drop;
    fuel_used = fuel_used + fuel_drop;

    fprintf('      200 km parking orbit reached: altitude %.2f km, dV %.4f m/s, fuel %.4f kg\n', ...
            (norm(X_chaser(1:3)) - sys.Re)/1000, dV_drop, fuel_drop);

    n_chaser = orbit_rate_from_state_local(X_chaser);
    n_target = orbit_rate_from_state_local(X_target);
    if abs(n_chaser - n_target) > 1e-12
        default_rbar_wait = 2 * 2*pi / abs(n_chaser - n_target);
    else
        default_rbar_wait = 2 * 86400;
    end

    dt_rbar = get_phase3_param_local(custom_params, 'dt_rbar_wait', 30);
    max_rbar_wait = get_phase3_param_local(custom_params, 'max_rbar_wait', default_rbar_wait);
    rbar_vbar_tol = get_phase3_param_local(custom_params, 'rbar_vbar_tol', 10);
    rbar_radial_tol = get_phase3_param_local(custom_params, 'rbar_radial_tol', 50e3);

    fprintf('   R_BAR_200_FPA leg 2: waiting for target R-bar alignment below the station...\n');
    [X_chaser, X_target, hist_wait, wait_time, rel_rbar] = ...
        propagate_until_rbar_below_local(X_chaser, X_target, sys, max_rbar_wait, dt_rbar, rbar_vbar_tol, rbar_radial_tol, hist.time_end);

    hist = append_phase3_hist_local(hist, hist_wait);
    info.rbar_wait_time = wait_time;
    info.rbar_rel_lvlh = rel_rbar;

    fprintf('      R-bar alignment after %.2f min, LVLH rel-pos [%+.1f, %+.1f, %+.1f] m\n', ...
            wait_time/60, rel_rbar(1), rel_rbar(2), rel_rbar(3));

    [target_reentry_r, fpa_calc] = compute_reentry_target_radius_local(norm(X_chaser(1:3)), sys);
    info.fpa_calc = fpa_calc;
    info.reentry_target_radius = target_reentry_r;
    charge_final_fuel = get_phase3_bool_param_local(custom_params, 'charge_final_reentry_fuel', false);

    fprintf('   R_BAR_200_FPA leg 3: injection toward 120 km / %.2f deg FPA (fuel update: %s)...\n', ...
            rad2deg(sys.reentry_flight_path_angle), bool_label_local(charge_final_fuel));
    [X_chaser, dV_entry, fuel_entry] = apply_reentry_departure_impulse_local(X_chaser, target_reentry_r, sys, charge_final_fuel);
    dV_used = dV_used + dV_entry;
    if charge_final_fuel
        fuel_used = fuel_used + fuel_entry;
    end

    info.final_injection_dV = dV_entry;
    info.final_injection_fuel = fuel_entry;
    info.final_injection_fuel_charged = charge_final_fuel;

    hist = log_phase3_state_local(hist, X_chaser, X_target, hist.time_end);

    r_start = norm(X_chaser(1:3));
    a_trans = 0.5 * (r_start + target_reentry_r);
    default_reentry_time = 1.2 * pi * sqrt(a_trans^3 / sys.mu);
    dt_reentry = get_phase3_param_local(custom_params, 'dt_reentry_coast', 2);
    max_reentry_time = get_phase3_param_local(custom_params, 'max_reentry_coast_time', default_reentry_time);

    [X_chaser, X_target, hist_entry, coast_time] = ...
        propagate_until_altitude_local(X_chaser, X_target, sys, 120e3, max_reentry_time, dt_reentry, hist.time_end);

    hist = append_phase3_hist_local(hist, hist_entry);
    info.reentry_coast_time = coast_time;

    if charge_final_fuel
        fprintf('      120 km interface reached after %.2f min; final-injection fuel charged: %.4f kg.\n', ...
                coast_time/60, fuel_entry);
    else
        fprintf('      120 km interface reached after %.2f min; fuel-equivalent %.4f kg was not charged.\n', ...
                coast_time/60, fuel_entry);
    end
end

function [target_r, details] = compute_reentry_target_radius_local(start_r, sys)
    interface_r = sys.Re + 120e3;
    gamma = sys.reentry_flight_path_angle;
    phi = 2 * atan(start_r / (start_r - interface_r) * tan(gamma));
    denom = sin(phi) * cos(gamma) - cos(phi) * sin(gamma);

    if abs(denom) < 1e-12
        error('Re-entry FPA geometry became singular. Check start radius and flight-path angle.');
    end

    ecc = sin(gamma) / denom;
    target_r = (1 - ecc) / (1 + ecc) * start_r;

    if ~isfinite(target_r) || target_r <= 0
        error('Computed invalid re-entry target radius %.6g m.', target_r);
    end

    details = struct();
    details.start_radius = start_r;
    details.interface_radius = interface_r;
    details.flight_path_angle = gamma;
    details.phi = phi;
    details.eccentricity = ecc;
    details.target_radius = target_r;
end

function [X_chaser, dV_mag, fuel_used] = apply_reentry_departure_impulse_local(X_chaser, target_r, sys, charge_fuel)
    r_current = norm(X_chaser(1:3));
    v_current = norm(X_chaser(4:6));
    a_trans = 0.5 * (r_current + target_r);
    sqrt_arg = sys.mu * (2/r_current - 1/a_trans - sys.J2 * sys.Re^2 / r_current^3 * (3 * (X_chaser(3)/r_current)^2 - 1));

    if sqrt_arg <= 0
        error('Computed invalid re-entry injection speed. sqrt argument = %.6g.', sqrt_arg);
    end

    v_trans1 = sqrt(sqrt_arg);
    dV_cmd = v_trans1 - v_current;
    dV_mag = abs(dV_cmd);

    X_chaser(4:6) = X_chaser(4:6) + dV_cmd * X_chaser(4:6) / v_current;

    m0 = X_chaser(14);
    m1 = m0 * exp(-dV_mag / (sys.Isp * sys.g0));
    fuel_used = m0 - m1;
    if charge_fuel
        X_chaser(14) = m1;
    end
end

function [X_chaser, X_target, hist, elapsed, rel_lvlh] = propagate_until_rbar_below_local(X_chaser, X_target, sys, max_wait, dt, vbar_tol, radial_tol, t0)
    hist = init_phase3_hist_local();
    elapsed = 0;
    [rel_prev, ~] = rel_state_lvlh_local(X_chaser, X_target);

    if is_below_target_rbar_local(rel_prev, X_chaser, X_target, radial_tol) && abs(rel_prev(2)) <= vbar_tol
        rel_lvlh = rel_prev;
        hist = log_phase3_state_local(hist, X_chaser, X_target, t0);
        return;
    end

    while elapsed < max_wait - 1e-12
        dt_step = min(dt, max_wait - elapsed);
        X_prev = X_chaser;
        T_prev = X_target;
        y_prev = rel_prev(2);

        [X_next, T_next] = rk4_pair_step_full_local(X_prev, T_prev, dt_step, sys, t0 + elapsed);
        [rel_next, ~] = rel_state_lvlh_local(X_next, T_next);
        y_next = rel_next(2);

        crossed_vbar = (y_prev == 0) || (y_next == 0) || (sign(y_prev) ~= sign(y_next));
        below_target_rbar = is_below_target_rbar_local(rel_next, X_next, T_next, radial_tol);
        if below_target_rbar && (abs(y_next) <= vbar_tol || crossed_vbar)
            if crossed_vbar && y_prev ~= 0 && y_next ~= 0
                [X_chaser, X_target, t_cross] = refine_vbar_crossing_local(X_prev, T_prev, y_prev, dt_step, sys, t0 + elapsed);
                elapsed = elapsed + t_cross;
            else
                X_chaser = X_next;
                X_target = T_next;
                elapsed = elapsed + dt_step;
            end

            [rel_lvlh, ~] = rel_state_lvlh_local(X_chaser, X_target);
            hist = log_phase3_state_local(hist, X_chaser, X_target, t0 + elapsed);
            return;
        end

        elapsed = elapsed + dt_step;
        X_chaser = X_next;
        X_target = T_next;
        rel_prev = rel_next;
        hist = log_phase3_state_local(hist, X_chaser, X_target, t0 + elapsed);
    end

    error('R_BAR_200_FPA did not reach target R-bar alignment within %.2f hours.', max_wait/3600);
end

function tf = is_below_target_rbar_local(rel_lvlh, X_chaser, X_target, radial_tol)
    radial_gap = norm(X_target(1:3)) - norm(X_chaser(1:3));
    tf = rel_lvlh(1) < 0 && abs(rel_lvlh(1) + radial_gap) <= radial_tol;
end

function [X_cross, T_cross, t_cross] = refine_vbar_crossing_local(X0, T0, y0, dt_window, sys, t_abs0)
    lo = 0;
    hi = dt_window;
    X_cross = X0;
    T_cross = T0;
    t_cross = 0;

    for iter = 1:40
        mid = 0.5 * (lo + hi);
        [X_mid, T_mid] = propagate_pair_state_only_full_local(X0, T0, mid, min(2, max(mid/10, 0.1)), sys, t_abs0);
        [rel_mid, ~] = rel_state_lvlh_local(X_mid, T_mid);

        if sign(rel_mid(2)) == sign(y0)
            lo = mid;
        else
            hi = mid;
            X_cross = X_mid;
            T_cross = T_mid;
            t_cross = mid;
        end
    end
end

function [X_chaser, X_target, hist, elapsed] = propagate_until_altitude_local(X_chaser, X_target, sys, target_alt, max_time, dt, t0)
    hist = init_phase3_hist_local();
    elapsed = 0;
    alt_prev = norm(X_chaser(1:3)) - sys.Re;

    if alt_prev <= target_alt
        hist = log_phase3_state_local(hist, X_chaser, X_target, t0);
        return;
    end

    while elapsed < max_time - 1e-12
        dt_step = min(dt, max_time - elapsed);
        X_prev = X_chaser;
        T_prev = X_target;

        [X_next, T_next] = rk4_pair_step_full_local(X_prev, T_prev, dt_step, sys, t0 + elapsed);
        alt_next = norm(X_next(1:3)) - sys.Re;

        if alt_next <= target_alt
            [X_chaser, X_target, t_cross] = refine_altitude_crossing_local(X_prev, T_prev, target_alt, dt_step, sys, t0 + elapsed);
            elapsed = elapsed + t_cross;
            hist = log_phase3_state_local(hist, X_chaser, X_target, t0 + elapsed);
            return;
        end

        elapsed = elapsed + dt_step;
        X_chaser = X_next;
        X_target = T_next;
        hist = log_phase3_state_local(hist, X_chaser, X_target, t0 + elapsed);
    end

    error('Re-entry coast did not reach %.1f km altitude within %.2f min.', target_alt/1000, max_time/60);
end

function [X_cross, T_cross, t_cross] = refine_altitude_crossing_local(X0, T0, target_alt, dt_window, sys, t_abs0)
    lo = 0;
    hi = dt_window;
    X_cross = X0;
    T_cross = T0;
    t_cross = 0;

    for iter = 1:40
        mid = 0.5 * (lo + hi);
        [X_mid, T_mid] = propagate_pair_state_only_full_local(X0, T0, mid, min(1, max(mid/10, 0.05)), sys, t_abs0);
        alt_mid = norm(X_mid(1:3)) - sys.Re;

        if alt_mid > target_alt
            lo = mid;
        else
            hi = mid;
            X_cross = X_mid;
            T_cross = T_mid;
            t_cross = mid;
        end
    end
end

function [X_chaser, X_target] = propagate_pair_state_only_full_local(X_chaser, X_target, duration, dt, sys, t_abs0)
    if duration <= 0
        return;
    end

    elapsed = 0;
    while elapsed < duration - 1e-12
        dt_step = min(dt, duration - elapsed);
        [X_chaser, X_target] = rk4_pair_step_full_local(X_chaser, X_target, dt_step, sys, t_abs0 + elapsed);
        elapsed = elapsed + dt_step;
    end
end

function [X_chaser, X_target] = rk4_pair_step_full_local(X_chaser, X_target, dt_step, sys, t_abs)
    X_t_state = [X_target; zeros(7,1); sys.Target_Mass];
    k1_t = Env_EOM(t_abs,             X_t_state,                [0;0;0], [0;0;0], sys, false);
    k2_t = Env_EOM(t_abs+dt_step/2,   X_t_state+k1_t*dt_step/2, [0;0;0], [0;0;0], sys, false);
    k3_t = Env_EOM(t_abs+dt_step/2,   X_t_state+k2_t*dt_step/2, [0;0;0], [0;0;0], sys, false);
    k4_t = Env_EOM(t_abs+dt_step,     X_t_state+k3_t*dt_step,   [0;0;0], [0;0;0], sys, false);
    X_target = X_target + (dt_step/6)*(k1_t(1:6) + 2*k2_t(1:6) + 2*k3_t(1:6) + k4_t(1:6));

    k1 = Env_EOM(t_abs,             X_chaser,             [0;0;0], [0;0;0], sys, true);
    k2 = Env_EOM(t_abs+dt_step/2,   X_chaser+k1*dt_step/2,[0;0;0], [0;0;0], sys, true);
    k3 = Env_EOM(t_abs+dt_step/2,   X_chaser+k2*dt_step/2,[0;0;0], [0;0;0], sys, true);
    k4 = Env_EOM(t_abs+dt_step,     X_chaser+k3*dt_step,  [0;0;0], [0;0;0], sys, true);
    X_chaser = X_chaser + (dt_step/6)*(k1 + 2*k2 + 2*k3 + k4);

    if norm(X_chaser(7:10)) > 0
        X_chaser(7:10) = X_chaser(7:10) / norm(X_chaser(7:10));
    end
end

function n = orbit_rate_from_state_local(X)
    r = X(1:3);
    v = X(4:6);
    n = norm(cross(r, v)) / norm(r)^2;
end

function value = get_phase3_param_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
end

function value = get_phase3_bool_param_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = parse_bool_setting_local(s.(name));
    else
        value = default_value;
    end
end

function value = parse_bool_setting_local(raw_value)
    if islogical(raw_value)
        value = raw_value;
        return;
    end

    if isnumeric(raw_value) && isscalar(raw_value)
        value = raw_value ~= 0;
        return;
    end

    text_value = lower(strtrim(string(raw_value)));
    if any(text_value == ["1", "true", "yes", "y", "on"])
        value = true;
    elseif any(text_value == ["0", "false", "no", "n", "off"])
        value = false;
    else
        error('Boolean setting must be true/false, on/off, yes/no, or 1/0. Got: %s', string(raw_value));
    end
end

function label = bool_label_local(value)
    if value
        label = "on";
    else
        label = "off";
    end
end

function hist = init_phase3_hist_local()
    hist.pos = [];
    hist.vel = [];
    hist.mass = [];
    hist.time = [];
    hist.time_end = 0;
    hist.target_pos = [];
    hist.target_vel = [];
    hist.rel_pos = [];
    hist.rel_pos_lvlh = [];
    hist.rel_vel_lvlh = [];
end

function hist = log_phase3_state_local(hist, X_chaser, X_target, t)
    if ~isfield(hist, 'time')
        hist = init_phase3_hist_local();
    end

    hist.pos = [hist.pos, X_chaser(1:3)];
    hist.vel = [hist.vel, X_chaser(4:6)];
    hist.mass = [hist.mass, X_chaser(14)];
    hist.time = [hist.time, t];
    hist.time_end = t;

    if ~isempty(X_target)
        hist.target_pos = [hist.target_pos, X_target(1:3)];
        hist.target_vel = [hist.target_vel, X_target(4:6)];
        hist.rel_pos = [hist.rel_pos, X_chaser(1:3) - X_target(1:3)];
        [rel_lvlh, rel_vel_lvlh] = rel_state_lvlh_local(X_chaser, X_target);
        hist.rel_pos_lvlh = [hist.rel_pos_lvlh, rel_lvlh];
        hist.rel_vel_lvlh = [hist.rel_vel_lvlh, rel_vel_lvlh];
    end
end

function hist = append_phase3_hist_local(hist, sub_hist)
    if isempty(sub_hist) || ~isfield(sub_hist, 'time') || isempty(sub_hist.time)
        return;
    end

    if ~isfield(hist, 'time') || isempty(hist.time)
        hist = sub_hist;
        if ~isfield(hist, 'time_end')
            hist.time_end = hist.time(end);
        end
        return;
    end

    fields = {'pos','vel','mass','time','target_pos','target_vel','rel_pos','rel_pos_lvlh','rel_vel_lvlh'};
    for ii = 1:numel(fields)
        f = fields{ii};
        if ~isfield(hist, f)
            hist.(f) = [];
        end
        if isfield(sub_hist, f) && ~isempty(sub_hist.(f))
            hist.(f) = [hist.(f), sub_hist.(f)];
        end
    end
    hist.time_end = hist.time(end);
end
