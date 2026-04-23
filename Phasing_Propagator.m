function [X_final, dV_used, fuel_used, hist] = Phasing_Propagator(sys, X0, target_r, mode, custom_params)
    % 상태 벡터: X = [x; y; z; vx; vy; vz; mass]
    m0 = X0(14);
    r0 = X0(1:3); v0 = X0(4:6);
    X_state = [r0; v0; m0];
    
    dV_used = 0;
    
    % 데이터 기록용 변수 초기화
    hist.pos = [];
    hist.mass = [];
    
    % --- 1. HOHMANN 모드 ---
    if mode == "HOHMANN"
        [X_state, dV, sub_hist] = execute_hohmann(sys, X_state, target_r);
        dV_used = dV_used + dV;
        hist.pos = [hist.pos, sub_hist.pos];
        hist.mass = [hist.mass, sub_hist.mass];
        
    % --- 2. LAMBERT 모드 (특정 시간 TOF 만에 타겟팅) ---
    elseif mode == "LAMBERT"
        TOF = custom_params.TOF; 
        target_pos = custom_params.target_pos; 
        fprintf('   * Lambert Problem 해석 (TOF: %.1f 초)...\n', TOF);
        
        [v1_req, ~] = solve_lambert(X_state(1:3), target_pos, TOF, sys.mu);
        
        dV1_vec = v1_req - X_state(4:6);
        dV1_mag = norm(dV1_vec) * (1 + randn * sys.noise.thrust_err);
        X_state(4:6) = X_state(4:6) + dV1_vec/norm(dV1_vec) * dV1_mag;
        X_state(7) = X_state(7) - X_state(7) * (1 - exp(-dV1_mag / (sys.Isp * sys.g0)));
        dV_used = dV_used + dV1_mag;
        
        dt = 10;
        t_span = 0:dt:TOF;
        for k = 1:length(t_span)-1
            k1 = orbit_dynamics(X_state, sys, 0, 1);
            k2 = orbit_dynamics(X_state + k1*dt/2, sys, 0, 1);
            k3 = orbit_dynamics(X_state + k2*dt/2, sys, 0, 1);
            k4 = orbit_dynamics(X_state + k3*dt, sys, 0, 1);
            X_state = X_state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
            
            % 상태 기록
            hist.pos = [hist.pos, X_state(1:3)];
            hist.mass = [hist.mass, X_state(7)];
        end
        
        v_current = X_state(4:6);
        v_circular = sqrt(sys.mu / target_r);
        h_vec = cross(X_state(1:3), v_current);
        v_target_dir = cross(h_vec, X_state(1:3));
        v_target_req = v_target_dir / norm(v_target_dir) * v_circular;
        
        dV2_vec = v_target_req - v_current;
        dV2_mag = norm(dV2_vec) * (1 + randn * sys.noise.thrust_err);
        X_state(4:6) = X_state(4:6) + dV2_vec/norm(dV2_vec) * dV2_mag;
        X_state(7) = X_state(7) - X_state(7) * (1 - exp(-dV2_mag / (sys.Isp * sys.g0)));
        dV_used = dV_used + dV2_mag;
    end
    
    % 시간 배열 생성 (모든 루프의 dt는 10초로 고정되어 있음)
    hist.time = (1:size(hist.pos, 2)) * 10; 
    
    fuel_used = m0 - X_state(7);
    X_final = zeros(14,1);
    X_final(1:6) = X_state(1:6);
    X_final(7:10) = [0;0;0;1]; 
    X_final(14) = X_state(7);
end

%% --- Helper: Hohmann Transfer Execution ---
function [X_st, dV_tot, sub_hist] = execute_hohmann(sys, X_st, target_r)
    sub_hist.pos = []; sub_hist.mass = [];
    r_current = norm(X_st(1:3));
    is_ascending = target_r > r_current;
    
    v1 = norm(X_st(4:6));
    a_trans = (r_current + target_r) / 2;
    v_trans1 = sqrt(sys.mu * (2/r_current - 1/a_trans));
    
    dV1_mag = abs(v_trans1 - v1) * (1 + randn * sys.noise.thrust_err);
    if is_ascending, X_st(4:6) = X_st(4:6)/v1 * (v1 + dV1_mag); 
    else, X_st(4:6) = X_st(4:6)/v1 * (v1 - dV1_mag); end
    X_st(7) = X_st(7) - X_st(7)*(1 - exp(-dV1_mag / (sys.Isp * sys.g0)));
    
    TOF = pi * sqrt(a_trans^3 / sys.mu);
    dt = 10;
    for k = 1:length(0:dt:TOF)-1
        k1 = orbit_dynamics(X_st, sys, 0, 1);
        k2 = orbit_dynamics(X_st + k1*dt/2, sys, 0, 1);
        k3 = orbit_dynamics(X_st + k2*dt/2, sys, 0, 1);
        k4 = orbit_dynamics(X_st + k3*dt, sys, 0, 1);
        X_st = X_st + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        
        % 상태 기록
        sub_hist.pos = [sub_hist.pos, X_st(1:3)];
        sub_hist.mass = [sub_hist.mass, X_st(7)];
    end
    
    v_curr = norm(X_st(4:6));
    v_circ = sqrt(sys.mu / norm(X_st(1:3)));
    
    dV2_mag = abs(v_circ - v_curr) * (1 + randn * sys.noise.thrust_err);
    X_st(4:6) = X_st(4:6) / v_curr * v_circ;
    X_st(7) = X_st(7) - X_st(7)*(1 - exp(-dV2_mag / (sys.Isp * sys.g0)));
    dV_tot = dV1_mag + dV2_mag;
end

%% --- Helper: Universal Variable Lambert Solver ---
function [v1, v2] = solve_lambert(r1, r2, TOF, mu)
    r1_mag = norm(r1); r2_mag = norm(r2);
    cos_dnu = dot(r1, r2) / (r1_mag * r2_mag);
    A = sqrt(r1_mag * r2_mag * (1 + cos_dnu));
    if A == 0, error('궤도가 180도입니다. Lambert 솔버 예외 처리 필요.'); end
    z = 0; C = 1/2; S = 1/6; 
    for i = 1:100
        y = r1_mag + r2_mag - A * (1 - z*S) / sqrt(C);
        x = sqrt(y / C);
        t_calc = (x^3 * S + A * sqrt(y)) / sqrt(mu);
        if abs(t_calc - TOF) < 1e-5, break; end
        z = z + (TOF - t_calc) * 0.1; 
        if z > 0
            S = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3; C = (1 - cos(sqrt(z))) / z;
        elseif z < 0
            S = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3; C = (cosh(sqrt(-z)) - 1) / (-z);
        else
            S = 1/6; C = 1/2;
        end
    end
    f = 1 - y/r1_mag; g = A * sqrt(y/mu); g_dot = 1 - y/r2_mag;
    v1 = (r2 - f*r1) / g; v2 = (g_dot*r2 - r1) / g;
end

%% --- Helper: Physics Engine ---
function dX = orbit_dynamics(X, sys, thrust_mag, dir)
    r = X(1:3); v = X(4:6); m = X(7);
    r_norm = norm(r); v_norm = norm(v);
    a_g = -sys.mu / r_norm^3 * r;
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1); (r(2)/r_norm)*(5*z2 - 1); (r(3)/r_norm)*(5*z2 - 3) ];
    a_thrust = [0;0;0]; dm = 0;
    if thrust_mag > 0
        actual_thrust = thrust_mag;
        a_thrust = (v / v_norm) * dir * (actual_thrust / m);
        dm = -actual_thrust / (sys.Isp * sys.g0);
    end
    dX = [v; a_g + a_j2 + a_thrust; dm];
end

function out = ifelse(cond, a, b)
    if cond, out = a; else, out = b; end
end