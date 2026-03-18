function sys = Mission_Config()
    % 환경 변수
    sys.mu = 3.986004418e14;   % 지구 중력 상수 [m^3/s^2]
    sys.Re = 6378137;          % 지구 적도 반지름 [m]
    sys.J2 = 1.08263e-3;       % J2 Perturbation 계수
    sys.g0 = 9.80665;          % 표준 중력 가속도 [m/s^2]

    % 위성 제원 (Target & Chaser)
    sys.Target_Mass = 2000;    % [kg]
    sys.Chaser_Mass_Init = 2000; % 초기 질량 (연료 포함) [kg]
    sys.Inertia = diag([800, 800, 600]); % 관성 모멘트 [kg*m^2]
    
    % 궤도 초기 조건
    sys.h_target = 500e3;      % 500 km
    sys.h_insert = 300e3;      % 300 km
    sys.h_wait   = 495e3;      % 495 km (Waiting Point)
    sys.h_reentry = 200e3;     % 200 km
    sys.inc = pi/2;            % 극궤도 (90 deg)
    
    % 추진 시스템 사양
    sys.Isp = 200;             % 비추력 [s]
    sys.Thrust_Impulsive = 300;% 고추력기 [N]
    sys.Thrust_Cont = 1;       % 저추력기 [N]
    
    % 노이즈 및 불확실성 모델 (Gaussian Noise 1-sigma)
    sys.noise.pos = 2.0;       % 위치 센서 노이즈 [m]
    sys.noise.vel = 0.05;      % 속도 센서 노이즈 [m/s]
    sys.noise.att = 0.001;     % 자세(Quaternion) 노이즈
    sys.noise.gyro = 1e-4;     % 각속도 노이즈 [rad/s]
    sys.noise.thrust_err = 0.02; % 추력 오차 (2%)
    
    % 제어 게인 (PD Controller)
    sys.GNC.Kp_pos = 0.05; sys.GNC.Kd_pos = 0.8;
    sys.GNC.Kp_att = 50;   sys.GNC.Kd_att = 100;
end