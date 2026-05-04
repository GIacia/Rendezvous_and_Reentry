function [X_final, dV_used, fuel_used, hist, X_target_final] = Phasing_Propagator(sys, X0, target_r, mode, custom_params, X_target0, pmode)
    % Phasing_Propagator
    % 3-DOF impulsive phasing propagator with optional target co-propagation.
    %
    % HOHMANN mode is now target-aware if X_target0 is supplied:
    %   1) wait on the insertion orbit until a J2-aware phase search says the
    %      Hohmann arrival will be close to the requested LVLH capture point,
    %   2) execute the Hohmann transfer,
    %   3) return chaser and target states at the same final epoch.
    %
    % X0         : Chaser 14-state [r; v; q; w; mass]
    % X_target0 : Optional target 6-state [r; v]

    if nargin < 5 || isempty(custom_params)
        custom_params = struct();
    end

    m0 = X0(14);
    X_state = [X0(1:3); X0(4:6); m0];   % [r; v; mass]

    propagate_target = (nargin >= 6) && ~isempty(X_target0);
    if propagate_target
        X_target_state = X_target0(1:6);
    else
        X_target_state = [];
    end

    dV_used = 0;
    hist = init_hist();

    if mode == "HOHMANN"
        [X_state, X_target_state, dV, sub_hist] = execute_hohmann(sys, X_state, target_r, X_target_state, custom_params, pmode);
        dV_used = dV_used + dV;
        hist = append_hist(hist, sub_hist);

    elseif mode == "LAMBERT"
        TOF = custom_params.TOF;
        target_pos = custom_params.target_pos;
        fprintf('   * Lambert Problem 해석 (TOF: %.1f 초)...\n', TOF);

        [v1_req, ~] = solve_lambert(X_state(1:3), target_pos, TOF, sys.mu);

        dV1_vec = v1_req - X_state(4:6);
        dV1_mag_cmd = norm(dV1_vec);
        dV1_mag_actual = dV1_mag_cmd * (1 + randn * sys.noise.thrust_err);
        if dV1_mag_cmd > 0
            X_state(4:6) = X_state(4:6) + dV1_vec/dV1_mag_cmd * dV1_mag_actual;
        end
        X_state(7) = X_state(7) * exp(-dV1_mag_actual / (sys.Isp * sys.g0));
        dV_used = dV_used + dV1_mag_actual;

        dt = get_param(custom_params, 'dt_transfer', 10);
        [X_state, X_target_state, sub_hist] = propagate_for_duration(X_state, X_target_state, sys, TOF, dt, 0);
        hist = append_hist(hist, sub_hist);

        v_current = X_state(4:6);
        v_circular = sqrt(sys.mu / target_r);
        h_vec = cross(X_state(1:3), v_current);
        v_target_dir = cross(h_vec, X_state(1:3));
        v_target_req = v_target_dir / norm(v_target_dir) * v_circular;

        dV2_vec = v_target_req - v_current;
        dV2_mag_cmd = norm(dV2_vec);
        dV2_mag_actual = dV2_mag_cmd * (1 + randn * sys.noise.thrust_err);
        if dV2_mag_cmd > 0
            X_state(4:6) = X_state(4:6) + dV2_vec/dV2_mag_cmd * dV2_mag_actual;
        end
        X_state(7) = X_state(7) * exp(-dV2_mag_actual / (sys.Isp * sys.g0));
        dV_used = dV_used + dV2_mag_actual;
    else
        error('Unknown phasing mode: %s', mode);
    end

    fuel_used = m0 - X_state(7);
    X_final = zeros(14,1);
    X_final(1:6) = X_state(1:6);
    X_final(7:10) = [0;0;0;1];
    X_final(11:13) = [0;0;0];
    X_final(14) = X_state(7);

    if propagate_target
        X_target_final = X_target_state;
    else
        X_target_final = [];
    end
end

%% --- Helper: Hohmann Transfer Execution ---
function [X_st, X_target_st, dV_tot, sub_hist] = execute_hohmann(sys, X_st, target_r, X_target_st, custom_params, pmode)
    sub_hist = init_hist();
    dV_tot = 0;

    dt_wait     = get_param(custom_params, 'dt_wait', 30);      % wait propagation step [s]
    dt_transfer = get_param(custom_params, 'dt_transfer', 10);  % transfer propagation step [s]
    dt_scan     = get_param(custom_params, 'dt_scan', 60);      % candidate spacing [s]
    max_wait    = get_param(custom_params, 'max_wait', 2*86400);% search horizon [s]
    refine_span = get_param(custom_params, 'refine_span', 180); % local refinement half-width [s]
    refine_step = get_param(custom_params, 'refine_step', 10);  % local refinement spacing [s]
    pos_tol     = get_param(custom_params, 'capture_pos_tol', 1000); % warning threshold [m]

    if isfield(sys, 'capture') && isfield(sys.capture, 'r_rel0')
        desired_rel_lvlh = sys.capture.r_rel0(:);
    else
        desired_rel_lvlh = [-5000; 0; 0];
    end

    r_current = norm(X_st(1:3));
    a_trans = (r_current + target_r) / 2;
    TOF = pi * sqrt(a_trans^3 / sys.mu); % Keplerian first guess; actual transfer is numerically propagated with J2

    if ~isempty(X_target_st)
        fprintf('   * J2-aware Hohmann phase search 시작...\n');
        best_wait = find_best_wait_time(sys, X_st, X_target_st, target_r, TOF, desired_rel_lvlh, ...
                                        dt_scan, max_wait, refine_span, refine_step, dt_transfer, pmode);
        fprintf('     - selected wait time: %.1f s (%.2f orbits at insertion altitude)\n', ...
                best_wait, best_wait / (2*pi*sqrt(r_current^3/sys.mu)));

        [X_st, X_target_st, wait_hist] = propagate_for_duration(X_st, X_target_st, sys, best_wait, dt_wait, 0);
        sub_hist = append_hist(sub_hist, wait_hist);
    else
        best_wait = 0;
    end

    % First Hohmann impulse at the numerically selected epoch.
    [X_st, dV1] = apply_hohmann_departure_impulse(X_st, target_r, X_target_st, sys, false);
    dV_tot = dV_tot + dV1;

    [X_st, X_target_st, transfer_hist] = propagate_for_duration(X_st, X_target_st, sys, TOF, dt_transfer, sub_hist.time_end);
    sub_hist = append_hist(sub_hist, transfer_hist);

    % Circularization impulse at arrival radius.
    [X_st, dV2] = apply_circularization_impulse(X_st, target_r, X_target_st, sys, false);
    dV_tot = dV_tot + dV2;

    if ~isempty(X_target_st)
        [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_st(1:6), X_target_st);
        miss = rel_lvlh - desired_rel_lvlh;
        fprintf('     - capture LVLH position: [%+.1f, %+.1f, %+.1f] m\n', rel_lvlh(1), rel_lvlh(2), rel_lvlh(3));
        fprintf('     - capture error from desired: %.1f m, rel-speed: %.4f m/s\n', norm(miss), norm(rel_vel_lvlh));
        if norm(miss) > pos_tol
            fprintf('     - warning: capture error is larger than %.0f m. Increase max_wait/refinement or use Lambert targeting.\n', pos_tol);
        end
    end
end

%% --- Helper: Find a wait time that best targets the LVLH capture point under J2 ---
function best_wait = find_best_wait_time(sys, X_ch0, X_t0, target_r, TOF, desired_rel_lvlh, dt_scan, max_wait, refine_span, refine_step, dt_transfer, pmode)
    candidate_times = 0:dt_scan:max_wait;
    [best_wait, ~] = scan_wait_candidates(sys, X_ch0, X_t0, target_r, TOF, desired_rel_lvlh, candidate_times, dt_transfer, pmode);

    t1 = max(0, best_wait - refine_span);
    t2 = min(max_wait, best_wait + refine_span);
    refined_times = t1:refine_step:t2;
    [best_wait, ~] = scan_wait_candidates(sys, X_ch0, X_t0, target_r, TOF, desired_rel_lvlh, refined_times, dt_transfer, pmode);
end

function [best_wait, best_metric] = scan_wait_candidates(sys, X_ch0, X_t0, target_r, TOF, desired_rel_lvlh, candidate_times, dt_transfer, pmode)
    best_wait = candidate_times(1);
    best_metric = inf;

    X_ch_wait = X_ch0;
    X_t_wait = X_t0;
    t_prev = 0;

    for idx = 1:numel(candidate_times)
        t_now = candidate_times(idx);
        dt_to_next = t_now - t_prev;
        if dt_to_next > 0
            [X_ch_wait, X_t_wait] = propagate_state_only(X_ch_wait, X_t_wait, sys, dt_to_next, min(60, max(1, dt_to_next)));
        end
        t_prev = t_now;
        if abs(acos(dot(X_ch_wait(1:3), X_t_wait(1:3))/(norm(X_ch_wait(1:3))*norm(X_t_wait(1:3)))) - sys.phase) > sys.phase && pmode == 1
            continue
        end

        [rel_lvlh, rel_vel_lvlh] = predict_hohmann_capture(sys, X_ch_wait, X_t_wait, target_r, TOF, dt_transfer);
        pos_error = norm(rel_lvlh - desired_rel_lvlh);
        % Position targeting is the main requirement. Relative speed is a soft penalty.
        metric = pos_error + 1000 * norm(rel_vel_lvlh);

        if metric < best_metric
            best_metric = metric;
            best_wait = t_now;
        end
    end
end

function [rel_lvlh, rel_vel_lvlh] = predict_hohmann_capture(sys, X_ch_wait, X_t_wait, target_r, TOF, dt_transfer)
    X_ch = X_ch_wait;
    X_t = X_t_wait;
    [X_ch, ~] = apply_hohmann_departure_impulse(X_ch, target_r, X_t, sys, false);
    [X_ch, X_t] = propagate_state_only(X_ch, X_t, sys, TOF, dt_transfer);
    [X_ch, ~] = apply_circularization_impulse(X_ch, target_r, X_t, sys, false);
    [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
end

%% --- Helper: Impulses ---
function [X_st, dV_mag] = apply_hohmann_departure_impulse(X_st, target_r, ~, sys, use_noise)
    if nargin < 4, use_noise = true; end
    r_current = norm(X_st(1:3));
    v_current = norm(X_st(4:6));
    a_trans = (r_current + target_r) / 2;
%    v_trans1 = sqrt(sys.mu * (2/r_current - 1/a_trans));
    v_trans1 = sqrt(sys.mu * ( -1/a_trans + 2/r_current - sys.J2 * sys.Re^2 / r_current^3 * (3 * (X_st(3)/r_current)^2 - 1)));

    dV_cmd = v_trans1 - v_current;
    dV_mag = abs(dV_cmd);
    if use_noise
        dV_mag = dV_mag * (1 + randn * sys.noise.thrust_err);
    end
    dv_vec = dV_cmd * X_st(4:6) / v_current;
    X_st(4:6) = X_st(4:6) + dv_vec;
    X_st(7) = X_st(7) * exp(-dV_mag / (sys.Isp * sys.g0));
end

function [X_st, dV_mag] = apply_circularization_impulse(X_st, target_r, X_target_st, sys, use_noise)
    if nargin < 3, use_noise = true; end
    r = X_st(1:3);
    v = X_st(4:6);
    r_target_current = norm(X_target_st(1:3));
    v_target_current = norm(X_target_st(4:6));
    r_norm = norm(r);
    h_vec = cross(r, v);
    tangential_dir = cross(h_vec, r);
    tangential_dir = tangential_dir / norm(tangential_dir);
    v_circ_req = tangential_dir * sqrt(sys.mu / r_norm);
    dv_vec_cmd = v_circ_req - v;
    dV_cmd = norm(dv_vec_cmd);
    dV_mag = dV_cmd;
    if use_noise
        dV_mag = dV_mag * (1 + randn * sys.noise.thrust_err);
    end
    if dV_cmd > 0
        X_st(4:6) = X_st(4:6) + dv_vec_cmd/dV_cmd * dV_mag;
    end
    X_st(7) = X_st(7) * exp(-dV_mag / (sys.Isp * sys.g0));
end

%% --- Helper: Propagation ---
function [X_chaser, X_target, hist] = propagate_for_duration(X_chaser, X_target, sys, duration, dt_nominal, t_offset)
    if nargin < 6, t_offset = 0; end
    hist = init_hist();
    elapsed = 0;
    while elapsed < duration - 1e-9
        dt = min(dt_nominal, duration - elapsed);
        [X_chaser, X_target] = rk4_step_chaser_target(X_chaser, X_target, sys, dt);
        elapsed = elapsed + dt;
        hist = log_full_state(hist, X_chaser, X_target, t_offset + elapsed);
    end
end

function [X_chaser, X_target] = propagate_state_only(X_chaser, X_target, sys, duration, dt_nominal)
    elapsed = 0;
    while elapsed < duration - 1e-9
        dt = min(dt_nominal, duration - elapsed);
        [X_chaser, X_target] = rk4_step_chaser_target(X_chaser, X_target, sys, dt);
        elapsed = elapsed + dt;
    end
end

function [X_chaser_next, X_target_next] = rk4_step_chaser_target(X_chaser, X_target, sys, dt)
    k1 = orbit_dynamics(X_chaser, sys, 0, 1);
    k2 = orbit_dynamics(X_chaser + k1*dt/2, sys, 0, 1);
    k3 = orbit_dynamics(X_chaser + k2*dt/2, sys, 0, 1);
    k4 = orbit_dynamics(X_chaser + k3*dt, sys, 0, 1);
    X_chaser_next = X_chaser + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    if ~isempty(X_target)
        k1t = orbit_dynamics_target(X_target, sys);
        k2t = orbit_dynamics_target(X_target + k1t*dt/2, sys);
        k3t = orbit_dynamics_target(X_target + k2t*dt/2, sys);
        k4t = orbit_dynamics_target(X_target + k3t*dt, sys);
        X_target_next = X_target + (dt/6)*(k1t + 2*k2t + 2*k3t + k4t);
    else
        X_target_next = [];
    end
end

function dX = orbit_dynamics_target(X, sys)
    r = X(1:3);
    v = X(4:6);
    r_norm = norm(r);

    a_g = -sys.mu / r_norm^3 * r;
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1);
                      (r(2)/r_norm)*(5*z2 - 1);
                      (r(3)/r_norm)*(5*z2 - 3) ];

    dX = [v; a_g + a_j2];
end

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

%% --- Helper: LVLH relative state ---
function [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_chaser6, X_target6)
    r_t = X_target6(1:3);
    v_t = X_target6(4:6);
    r_c = X_chaser6(1:3);
    v_c = X_chaser6(4:6);

    h_vec = cross(r_t, v_t);
    i_u = r_t / norm(r_t);
    k_u = h_vec / norm(h_vec);
    j_u = cross(k_u, i_u);
    C_I2L = [i_u'; j_u'; k_u'];

    rho_eci = r_c - r_t;
    rel_lvlh = C_I2L * rho_eci;

    omega_lvlh_eci = h_vec / norm(r_t)^2;
    rho_dot_eci = v_c - v_t - cross(omega_lvlh_eci, rho_eci);
    rel_vel_lvlh = C_I2L * rho_dot_eci;
end

%% --- Helper: History logging ---
function hist = init_hist()
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

function hist = log_full_state(hist, X_chaser, X_target, t)
    hist.pos = [hist.pos, X_chaser(1:3)];
    hist.vel = [hist.vel, X_chaser(4:6)];
    hist.mass = [hist.mass, X_chaser(7)];
    hist.time = [hist.time, t];
    hist.time_end = t;
    if ~isempty(X_target)
        hist.target_pos = [hist.target_pos, X_target(1:3)];
        hist.target_vel = [hist.target_vel, X_target(4:6)];
        hist.rel_pos = [hist.rel_pos, X_chaser(1:3) - X_target(1:3)];
        [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_chaser(1:6), X_target);
        hist.rel_pos_lvlh = [hist.rel_pos_lvlh, rel_lvlh];
        hist.rel_vel_lvlh = [hist.rel_vel_lvlh, rel_vel_lvlh];
    end
end

function hist = append_hist(hist, sub_hist)
    if isempty(sub_hist.time)
        return;
    end
    fields = {'pos','vel','mass','time','target_pos','target_vel','rel_pos','rel_pos_lvlh','rel_vel_lvlh'};
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(sub_hist, f) && isfield(hist, f)
            hist.(f) = [hist.(f), sub_hist.(f)];
        end
    end
    hist.time_end = hist.time(end);
end

function value = get_param(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
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
