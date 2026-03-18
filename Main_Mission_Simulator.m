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

%% 2. Phase 1: Phasing & Homing (300 km -> 495 km)
fprintf('[Phase 1] 3-DOF 물리 시뮬레이션 시작 (J2 섭동 및 노이즈 포함)...\n');

% 체이서 초기 상태 (300 km 극궤도)
v_insert = sqrt(sys.mu / (sys.Re + sys.h_insert));
X_chaser_init = [sys.Re+sys.h_insert; 0; 0;  0; v_insert; 0;  0;0;0;1;  0;0;0;  m_current];
target_radius = sys.Re + sys.h_wait;

% 변수 세팅
custom_params.N = 3; % MULTI-HOHMANN 시 호만 전이 횟수 설정
custom_params.TOF = 10000; % LAMBERT 시 비행 시간 설정 (예: 10000초)

% 타겟의 미래 위치를 임의로 지정 (LAMBERT 모드 타겟팅용)
% 495km 고도의 특정 y축 지점을 향해 날아간다고 가정
custom_params.target_pos = [0; sys.Re + sys.h_wait; 0]; 

% 원하는 모드로 문자열만 변경하세요: "1-HOHMANN", "MULTI-HOHMANN", "LAMBERT", "CONTINUOUS" ## 300N 고추력 또는 1N 저추력 나선상승
phasing_mode = "MULTI-HOHMANN"; 
% fprintf('   선택된 기동 방식: %s\n', phasing_mode);

[X_chaser, dV_p1, fuel_p1, hist_p1] = Phasing_Propagator(sys, X_chaser_init, target_radius, phasing_mode, custom_params);
m_current = X_chaser(14);

Budget = [Budget; {"Phase 1: "+phasing_mode, dV_p1, fuel_p1, m_current}];
fprintf('   도달 고도: %.2f km (목표: 495.00 km)\n', (norm(X_chaser(1:3)) - sys.Re)/1000);

% 타겟 위성 상태 세팅 (500 km)
v_t = sqrt(sys.mu / (sys.Re + sys.h_target));
X_target = [sys.Re+sys.h_target; 0; 0; 0; v_t; 0];

%% 3. Phase 2: R-bar Approach & Berthing (6-DOF Simulation)
fprintf('\n[Phase 2] Starting 6-DOF R-bar Approach (Glideslope Controlled)...\n');

% --- [수정 1] 완벽한 초기 진입 상태 세팅 (J2 보정 정밀 각속도 적용) ---
r_t = X_target(1:3); v_t = X_target(4:6);
h_vec = cross(r_t, v_t);
i_u = r_t/norm(r_t); k_u = h_vec/norm(h_vec); j_u = cross(k_u, i_u);
C_I2L = [i_u'; j_u'; k_u']; % ECI to LVLH 변환 행렬

offset_lvlh = [-5000; 0; 0]; % R-bar 5km 아래 (-5000m)
v_app_initial = [2.0; 0; 0]; % 초기 접근 속도(2m/s)

r_c_new = r_t + C_I2L' * offset_lvlh;

% J2 섭동이 포함된 Target의 정확한 ECI 각속도 도출 (단순 n_mean 사용 폐기)
w_lvlh_eci = h_vec / norm(r_t)^2; 

v_c_new = v_t + C_I2L'*v_app_initial + cross(w_lvlh_eci, C_I2L' * offset_lvlh);
X_chaser(1:6) = [r_c_new; v_c_new]; 

% 시뮬레이션 시간 
dt = 1; 
T_sim = 3500; 
t_vec = 0:dt:T_sim;

hist_pos = zeros(3, length(t_vec));
hist_mass = zeros(1, length(t_vec));

for k = 1:length(t_vec)
    % 1. 관측 및 로깅 (Chaser와 Target이 동일한 시간 't'에 있을 때 오차 계산)
    h_vec = cross(X_target(1:3), X_target(4:6));
    i_u = X_target(1:3)/norm(X_target(1:3)); k_u = h_vec/norm(h_vec); j_u = cross(k_u, i_u);
    C_I2L = [i_u'; j_u'; k_u'];
    
    hist_pos(:,k) = C_I2L * (X_chaser(1:3) - X_target(1:3));
    hist_mass(k) = X_chaser(14);
    
    % 종료 조건 판별
    if norm(hist_pos(:,k) - [-2;0;0]) < 0.5 && mode == "FINAL_BERTHING"
        fprintf('   -> Berthing Achieved at t = %d seconds.\n', t_vec(k));
        break;
    end
    
    % 2. GNC Control Logic (동일한 시간대에서 명령 생성)
    [F_cmd, T_cmd, mode] = GNC_Controller(X_chaser, X_target, sys);
    
    % 3. Target RK4 Propagation (이제서야 t+dt로 전파)
    X_t_state = [X_target; zeros(7,1); sys.Target_Mass]; 
    k1_t = Env_EOM(t_vec(k), X_t_state, [0;0;0], [0;0;0], sys, false);
    k2_t = Env_EOM(t_vec(k)+dt/2, X_t_state + k1_t*dt/2, [0;0;0], [0;0;0], sys, false);
    k3_t = Env_EOM(t_vec(k)+dt/2, X_t_state + k2_t*dt/2, [0;0;0], [0;0;0], sys, false);
    k4_t = Env_EOM(t_vec(k)+dt, X_t_state + k3_t*dt, [0;0;0], [0;0;0], sys, false);
    X_target = X_target + (dt/6)*(k1_t(1:6) + 2*k2_t(1:6) + 2*k3_t(1:6) + k4_t(1:6));
    
    % 4. Chaser RK4 Propagation (동시에 t+dt로 전파)
    k1 = Env_EOM(t_vec(k), X_chaser, F_cmd, T_cmd, sys, true);
    k2 = Env_EOM(t_vec(k)+dt/2, X_chaser + k1*dt/2, F_cmd, T_cmd, sys, true);
    k3 = Env_EOM(t_vec(k)+dt/2, X_chaser + k2*dt/2, F_cmd, T_cmd, sys, true);
    k4 = Env_EOM(t_vec(k)+dt, X_chaser + k3*dt, F_cmd, T_cmd, sys, true);
    X_chaser = X_chaser + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    X_chaser(7:10) = X_chaser(7:10) / norm(X_chaser(7:10));
end
hist_pos = hist_pos(:, 1:k);
hist_mass = hist_mass(1:k);

%% 4. Phase 3: Re-entry (500 km -> 200 km)
fprintf('\n[Phase 3] 3-DOF 재진입 시뮬레이션 시작 (J2 섭동 및 노이즈 역추진)...\n');

% 고추력으로 대기권 인터페이스에 하강
reentry_mode = "1-HOHMANN"; 
fprintf('   선택된 기동 방식: %s\n', reentry_mode);

target_reentry_r = sys.Re + sys.h_reentry;

% Phase 2 Berthing 직후의 X_chaser 상태를 그대로 입력하여 하강
[X_chaser, dV_p3, fuel_p3, hist_p3] = Phasing_Propagator(sys, X_chaser, target_reentry_r, reentry_mode, custom_params);
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
plot(hist_pos(2,1:k), hist_pos(1,1:k), 'b-', 'LineWidth', 2); hold on;
plot(0,0,'r^','MarkerSize',10,'MarkerFaceColor','r');
grid on; xlabel('V-bar (m)'); ylabel('R-bar (m)');
title('R-bar Forced Approach Trajectory (LVLH)');
legend('Chaser Trajectory', 'Target');
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
plot3(hist_p3.pos(1,:)/1e3, hist_p3.pos(2,:)/1e3, hist_p3.pos(3,:)/1e3, 'r-', 'LineWidth', 1.5); % Phase 3 (재진입)

% 타겟 궤도 (500km) 참조선 그리기
theta = linspace(0, 2*pi, 100);
plot3((sys.Re+sys.h_target)/1e3*cos(theta), (sys.Re+sys.h_target)/1e3*sin(theta), zeros(1,100), 'k--', 'LineWidth', 1);

title('3D Mission Trajectory (Earth Centered Inertial)');
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
legend('Earth', 'Phase 1: Phasing (Ascent)', 'Phase 3: De-orbit (Descent)', 'Target Orbit (500km)', 'Location', 'best');
view(45, 30); % 보기 좋은 각도 설정

% --- (2) Altitude vs Time Profile ---
subplot(2, 2, 2);
alt_p1 = vecnorm(hist_p1.pos) - sys.Re; % 중심으로부터의 거리 - 지구반지름 = 고도
alt_p3 = vecnorm(hist_p3.pos) - sys.Re;
time_p1 = hist_p1.time / 3600; % 시간(초)을 시간(Hour)으로 변환
time_p3 = hist_p3.time / 3600;

plot(time_p1, alt_p1/1e3, 'g-', 'LineWidth', 2); hold on;
% 재진입 시간은 Phasing 직후에 이어지도록 offset 설정
plot(time_p3 + time_p1(end), alt_p3/1e3, 'r-', 'LineWidth', 2); 
yline(sys.h_insert/1e3, 'k:', 'Insertion (300km)');
yline(sys.h_target/1e3, 'k--', 'Target (500km)');
yline(sys.h_reentry/1e3, 'k:', 'Re-entry (200km)');

title('Altitude Profile');
xlabel('Mission Time (Hours)'); ylabel('Altitude (km)');
legend('Phasing Maneuver', 'Re-entry Maneuver');
grid on;

% --- (3) Mass Depletion Profile ---
subplot(2, 2, 4);
plot(time_p1, hist_p1.mass, 'g-', 'LineWidth', 2); hold on;
plot(time_p3 + time_p1(end), hist_p3.mass, 'r-', 'LineWidth', 2);

title('Spacecraft Mass Depletion (Fuel Consumption)');
xlabel('Mission Time (Hours)'); ylabel('Mass (kg)');
legend('Fuel used in Phasing', 'Fuel used in Re-entry');
grid on;