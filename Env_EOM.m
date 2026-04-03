function dX = Env_EOM(t, X, F_cmd_body, T_cmd_body, sys, is_6DOF)
    % X state: [r(3); v(3); q(4); w(3); mass(1)] 
    r = X(1:3); v = X(4:6); mass = X(14);
    r_norm = norm(r);
    
    % 1. Gravity & J2 Perturbation
    a_g = -sys.mu / r_norm^3 * r;
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1);
                      (r(2)/r_norm)*(5*z2 - 1);
                      (r(3)/r_norm)*(5*z2 - 3) ];
                  
    % 2. Dynamics Mode Selection
    if is_6DOF
        q = X(7:10); w = X(11:13);
        
        % 추력 발생
        F_actual = F_cmd_body;
        T_actual = T_cmd_body;
        
        % Body -> ECI 변환
        C_B2I = quat2rotm_custom(q)'; 
        F_eci = C_B2I * F_actual;
        a_thrust = F_eci / mass;
        
        % Rotational Dynamics
        Omega = [0 -w(3) w(2) w(1); w(3) 0 -w(1) w(2); -w(2) w(1) 0 w(3); -w(1) -w(2) -w(3) 0];
        dq = 0.5 * Omega * q;
        dw = sys.Inertia \ (T_actual - cross(w, sys.Inertia * w));
        
        % Mass Flow Rate (m_dot = -F / (Isp * g0))
        dm = -norm(F_actual) / (sys.Isp * sys.g0);
    else
        % 3-DOF Mode (Attitude is ignored, Impulsive assumed outside)
        dq = zeros(4,1); dw = zeros(3,1);
        a_thrust = zeros(3,1); dm = 0;
    end
    
    dX = [v; a_g + a_j2 + a_thrust; dq; dw; dm];
end

function R = quat2rotm_custom(q)
    q0=q(4); q1=q(1); q2=q(2); q3=q(3);
    R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end