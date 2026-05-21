function [F_body, T_body, mode] = GNC_Controller(X_chaser, X_target, sys)
    % 1. Sensor Measurement
    r_c = X_chaser(1:3);
    v_c = X_chaser(4:6);
    q_c = X_chaser(7:10); 
    w_c = X_chaser(11:13);
    
    r_t = X_target(1:3); v_t = X_target(4:6);
    
    % 2. LVLH Relative Navigation
    h_vec = cross(r_t, v_t);
    i_unit = r_t / norm(r_t);
    k_unit = h_vec / norm(h_vec);
    j_unit = cross(k_unit, i_unit);
    C_I2L = [i_unit'; j_unit'; k_unit'];
    
    r_rel = C_I2L * (r_c - r_t);
    
    % J2가 반영된 정밀 LVLH 각속도
    w_lvlh_eci = h_vec / norm(r_t)^2;
    w_lvlh = C_I2L * w_lvlh_eci; 
    
    v_rel = C_I2L * (v_c - v_t) - cross(w_lvlh, r_rel);
    dist = norm(r_rel);
    
    % ---------------------------------------------------------------------
    % 3. [핵심] EXACT 비선형 환경 외란 계산 (J2 각가속도 반영)
    % ---------------------------------------------------------------------
    r_t_norm = norm(r_t);
    g_t = -sys.mu / r_t_norm^3 * r_t;
    z2_t = (r_t(3)/r_t_norm)^2;
    factor_t = 1.5 * sys.J2 * (sys.mu/r_t_norm^2) * (sys.Re/r_t_norm)^2;
    j2_t = factor_t * [ (r_t(1)/r_t_norm)*(5*z2_t - 1); (r_t(2)/r_t_norm)*(5*z2_t - 1); (r_t(3)/r_t_norm)*(5*z2_t - 3) ];
    a_t_eci = g_t + j2_t;
    
    % LVLH 각가속도 정밀 계산
    w_dot_eci = cross(r_t, a_t_eci)/r_t_norm^2 - 2*dot(r_t, v_t)/r_t_norm^4 * h_vec;
    w_dot_lvlh = C_I2L * w_dot_eci;
    
    r_c_norm = norm(r_c);
    g_c = -sys.mu / r_c_norm^3 * r_c;
    z2_c = (r_c(3)/r_c_norm)^2;
    factor_c = 1.5 * sys.J2 * (sys.mu/r_c_norm^2) * (sys.Re/r_c_norm)^2;
    j2_c = factor_c * [ (r_c(1)/r_c_norm)*(5*z2_c - 1); (r_c(2)/r_c_norm)*(5*z2_c - 1); (r_c(3)/r_c_norm)*(5*z2_c - 3) ];
    a_c_free_eci = g_c + j2_c;
    
    % 완벽한 물리 모델 기반의 LVLH 상대 환경 가속도 (모든 외력 항 포함)
    A_env = C_I2L * (a_c_free_eci - a_t_eci) - cross(w_dot_lvlh, r_rel) - 2*cross(w_lvlh, v_rel) - cross(w_lvlh, cross(w_lvlh, r_rel));
    % ---------------------------------------------------------------------
       
    % 4. Guidance & Control (Decoupled Glideslope)
    r_cmd = [-2; 0; 0]; % 최종 도킹 목표 (-2m)
    
    % 위치 제어 게인 (이탈 방향인 Y, Z를 강력하게 억제)
    Kp_outer = [0.002; 0.1; 0.1]; 
    v_cmd_raw = Kp_outer .* (r_cmd - r_rel);
    
    if dist > 10
        mode = "R-BAR_CLOSING";
        v_max_x = 2.0; 
        max_thrust = sys.Thrust_Impulsive; 
    else
        mode = "FINAL_BERTHING";
        v_max_x = 0.2; 
        max_thrust = sys.Thrust_Cont; 
    end
    
    % 속도 명령 독립 포화 제어 (Vector Saturation 방지)
    v_cmd = v_cmd_raw;
    if abs(v_cmd(1)) > v_max_x, v_cmd(1) = sign(v_cmd(1)) * v_max_x; end
    
    v_max_yz = 1.0; 
    if abs(v_cmd(2)) > v_max_yz, v_cmd(2) = sign(v_cmd(2)) * v_max_yz; end
    if abs(v_cmd(3)) > v_max_yz, v_cmd(3) = sign(v_cmd(3)) * v_max_yz; end
    
    % 속도 제어 게인 (Y, Z축 제동력 압도적 강화)
    Kd_inner = [0.05; 0.5; 0.5]; 
    a_pd = Kd_inner .* (v_cmd - v_rel);
    
    % 피드포워드(환경 완벽 상쇄) + 피드백(목표 추종)
    F_env = X_chaser(14) * (-A_env); 
    F_pd  = X_chaser(14) * a_pd;     
    
    % 추력 포화 방지 로직
    avail_thrust_for_pd = max_thrust - norm(F_env); 
    if norm(F_pd) > avail_thrust_for_pd && avail_thrust_for_pd > 0
        F_pd = (F_pd / norm(F_pd)) * avail_thrust_for_pd;
    elseif avail_thrust_for_pd <= 0
        F_pd = [0;0;0];
    end
    F_lvlh = F_env + F_pd;
    
    % 5. Force/Torque Mapping
    C_B2I = quat2rotm_custom(q_c)';
    F_eci = C_I2L' * F_lvlh;
    F_body = C_B2I' * F_eci;
    
    % 자세 제어
    T_body = -sys.GNC.Kd_att * w_c; 
end

function R = quat2rotm_custom(q)
    q0=q(4); q1=q(1); q2=q(2); q3=q(3);
    R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end