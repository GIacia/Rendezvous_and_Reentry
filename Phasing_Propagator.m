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

    elseif mode == "MULTI_HOHMANN"
        [X_state, X_target_state, dV, sub_hist] = execute_multi_hohmann(sys, X_state, target_r, X_target_state, custom_params, pmode);
        dV_used = dV_used + dV;
        hist = append_hist(hist, sub_hist);

    elseif mode == "CUSTOM_IMPULSE"
        [X_state, X_target_state, dV, sub_hist] = execute_custom_impulse(sys, X_state, X_target_state, custom_params);
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

%% --- Helper: Custom phase-triggered impulsive boost ---
function [X_st, X_target_st, dV_mag, sub_hist] = execute_custom_impulse(sys, X_st, X_target_st, custom_params)
    % User-defined phasing logic:
    %   1) propagate chaser/target until the signed phase angle reaches the input value,
    %   2) apply an impulsive delta-V tilted by gamma from local tangential direction,
    %   3) continue propagation until the requested LVLH relative position is reached
    %      or until the closest-approach point is detected.

    if isempty(X_target_st)
        error('CUSTOM_IMPULSE mode requires X_target0 = [r_target; v_target].');
    end

    phase_angle = get_angle_param_rad(custom_params, 'phase_angle', 'phase_angle_unit');
    delta_v     = get_required_scalar(custom_params, 'delta_v');
    gamma       = get_angle_param_rad(custom_params, 'gamma', 'gamma_unit');
    maneuver_opts = get_maneuver_options(sys, custom_params);

    desired_rel_lvlh = get_vector_param(custom_params, 'desired_rel_lvlh', [0; -5; 0]);

    time_step = get_param(custom_params, 'time_step', 10);
    dt_phase  = get_param(custom_params, 'dt_phase',  get_param(custom_params, 'dt_wait',     time_step));
    dt_cap    = get_param(custom_params, 'dt_capture',get_param(custom_params, 'dt_transfer', time_step));

    max_wait = get_param(custom_params, 'max_wait', 2*86400);
    default_max_capture = 2 * 2*pi * sqrt(norm(X_st(1:3))^3 / sys.mu);
    max_capture_time = get_param(custom_params, 'max_capture_time', default_max_capture);

    phase_tol = get_param(custom_params, 'phase_tol', 1e-7);        % [rad]
    event_time_tol = get_param(custom_params, 'event_time_tol', 1e-3); % [s]
    capture_pos_tol = get_param(custom_params, 'capture_pos_tol', 0.5); % [m]
    min_capture_time = get_param(custom_params, 'min_capture_time', 0); % [s]

    if delta_v < 0
        error('custom_params.delta_v must be non-negative [m/s].');
    end

    sub_hist = init_hist();

    fprintf('   * Custom phased maneuver mode started.\n');
    fprintf('     - burn model: %s\n', char(maneuver_opts.burn_model));
    fprintf('     - target signed phase angle: %.8f rad (%.6f deg)\n', phase_angle, rad2deg(phase_angle));
    fprintf('     - commanded delta-V: %.6f m/s, gamma: %.8f rad (%.6f deg)\n', delta_v, gamma, rad2deg(gamma));
    fprintf('     - desired LVLH rel-pos: [%+.3f, %+.3f, %+.3f] m\n', desired_rel_lvlh(1), desired_rel_lvlh(2), desired_rel_lvlh(3));

    [X_st, X_target_st, wait_hist, wait_time, final_phase_error] = ...
        propagate_to_phase_angle(sys, X_st, X_target_st, phase_angle, dt_phase, max_wait, phase_tol, event_time_tol);
    sub_hist = append_hist(sub_hist, wait_hist);

    fprintf('     - departure wait time: %.3f s, final phase error: %.3e rad\n', wait_time, final_phase_error);

    [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
        execute_custom_maneuver(sys, X_st, X_target_st, delta_v, gamma, maneuver_opts, wait_time);
    sub_hist = append_hist(sub_hist, burn_hist);

    [X_st, X_target_st, cap_hist, capture_time, miss, rel_lvlh, rel_vel_lvlh, reached_tol] = ...
        propagate_until_capture(sys, X_st, X_target_st, desired_rel_lvlh, dt_cap, max_capture_time, ...
                                capture_pos_tol, event_time_tol, min_capture_time, wait_time + burn_duration);
    sub_hist = append_hist(sub_hist, cap_hist);

    if reached_tol
        fprintf('     - capture tolerance reached.\n');
    else
        fprintf('     - closest-approach capture used. Check miss distance.\n');
    end
    fprintf('     - maneuver duration: %.3f s, delivered delta-V budget: %.6f m/s\n', burn_duration, dV_mag);
    fprintf('     - capture propagation time after maneuver: %.3f s\n', capture_time);
    fprintf('     - total custom phase duration: %.3f s\n', wait_time + burn_duration + capture_time);
    fprintf('     - final LVLH rel-pos: [%+.3f, %+.3f, %+.3f] m\n', rel_lvlh(1), rel_lvlh(2), rel_lvlh(3));
    fprintf('     - desired-position miss: %.6f m, rel-speed: %.6f m/s\n', miss, norm(rel_vel_lvlh));
    if miss > capture_pos_tol
        fprintf('     - warning: miss is larger than capture_pos_tol = %.3f m. Re-check phase_angle/delta_v/gamma.\n', capture_pos_tol);
    end
end

function [X_ch, X_t, hist, elapsed, final_err] = propagate_to_phase_angle(sys, X_ch, X_t, phase_angle, dt, max_wait, phase_tol, event_time_tol)
    hist = init_hist();
    elapsed = 0;
    final_err = phase_error_to_target(X_ch, X_t, phase_angle);

    if abs(final_err) <= phase_tol
        hist = log_full_state(hist, X_ch, X_t, elapsed);
        return;
    end

    X_prev = X_ch;
    T_prev = X_t;
    err_prev = final_err;
    t_prev = elapsed;

    while elapsed < max_wait - 1e-12
        dt_eff = min(dt, max_wait - elapsed);
        [X_next, T_next] = rk4_step_chaser_target(X_prev, T_prev, sys, dt_eff);
        err_next = phase_error_to_target(X_next, T_next, phase_angle);

        if is_phase_crossing(err_prev, err_next) || abs(err_next) <= phase_tol
            [X_ch, X_t, local_t, final_err] = refine_phase_crossing(sys, X_prev, T_prev, phase_angle, dt_eff, phase_tol, event_time_tol);
            elapsed = t_prev + local_t;
            hist = log_full_state(hist, X_ch, X_t, elapsed);
            return;
        end

        elapsed = elapsed + dt_eff;
        X_prev = X_next;
        T_prev = T_next;
        err_prev = err_next;
        t_prev = elapsed;
        hist = log_full_state(hist, X_prev, T_prev, elapsed);
    end

    error('Phase angle %.8f rad was not reached within max_wait = %.1f s. Last phase error = %.3e rad.', ...
          phase_angle, max_wait, err_prev);
end

function tf = is_phase_crossing(err_a, err_b)
    % Assumes dt is small enough that the wrapped phase error does not jump by pi.
    tf = (err_a == 0) || (err_b == 0) || (sign(err_a) ~= sign(err_b) && abs(err_a - err_b) < pi);
end

function [X_best, T_best, t_best, err_best] = refine_phase_crossing(sys, X0, T0, phase_angle, dt_window, phase_tol, event_time_tol)
    err0 = phase_error_to_target(X0, T0, phase_angle);
    [X1, T1] = propagate_state_only(X0, T0, sys, dt_window, min(1, max(dt_window/20, 1e-3)));
    err1 = phase_error_to_target(X1, T1, phase_angle);

    if abs(err0) <= phase_tol
        X_best = X0; T_best = T0; t_best = 0; err_best = err0; return;
    elseif abs(err1) <= phase_tol
        X_best = X1; T_best = T1; t_best = dt_window; err_best = err1; return;
    end

    if ~is_phase_crossing(err0, err1)
        [X_best, T_best, t_best, err_best] = refine_phase_min_abs(sys, X0, T0, phase_angle, dt_window, event_time_tol);
        return;
    end

    lo = 0; hi = dt_window;
    err_lo = err0;
    X_best = X1; T_best = T1; t_best = hi; err_best = err1;

    while (hi - lo) > event_time_tol
        mid = 0.5 * (lo + hi);
        [X_mid, T_mid] = propagate_state_only(X0, T0, sys, mid, min(1, max(mid/20, 1e-3)));
        err_mid = phase_error_to_target(X_mid, T_mid, phase_angle);

        if abs(err_mid) < abs(err_best)
            X_best = X_mid; T_best = T_mid; t_best = mid; err_best = err_mid;
        end
        if abs(err_mid) <= phase_tol
            return;
        end

        if is_phase_crossing(err_lo, err_mid)
            hi = mid;
        else
            lo = mid;
            err_lo = err_mid;
        end
    end

    [X_best, T_best] = propagate_state_only(X0, T0, sys, t_best, min(1, max(t_best/20, 1e-3)));
    err_best = phase_error_to_target(X_best, T_best, phase_angle);
end

function [X_best, T_best, t_best, err_best] = refine_phase_min_abs(sys, X0, T0, phase_angle, dt_window, event_time_tol)
    phi = (sqrt(5) - 1) / 2;
    a = 0; b = dt_window;
    c = b - phi*(b-a);
    d = a + phi*(b-a);
    fc = phase_abs_error_at_time(sys, X0, T0, phase_angle, c);
    fd = phase_abs_error_at_time(sys, X0, T0, phase_angle, d);
    while (b - a) > event_time_tol
        if fc > fd
            a = c; c = d; fc = fd;
            d = a + phi*(b-a);
            fd = phase_abs_error_at_time(sys, X0, T0, phase_angle, d);
        else
            b = d; d = c; fd = fc;
            c = b - phi*(b-a);
            fc = phase_abs_error_at_time(sys, X0, T0, phase_angle, c);
        end
    end
    t_best = 0.5*(a+b);
    [X_best, T_best] = propagate_state_only(X0, T0, sys, t_best, min(1, max(t_best/20, 1e-3)));
    err_best = phase_error_to_target(X_best, T_best, phase_angle);
end

function val = phase_abs_error_at_time(sys, X0, T0, phase_angle, t_eval)
    [X_eval, T_eval] = propagate_state_only(X0, T0, sys, t_eval, min(1, max(t_eval/20, 1e-3)));
    val = abs(phase_error_to_target(X_eval, T_eval, phase_angle));
end

function [X_ch, X_t, hist, capture_time, miss, rel_lvlh, rel_vel_lvlh, reached_tol] = ...
    propagate_until_capture(sys, X_ch, X_t, desired_rel_lvlh, dt, max_time, pos_tol, event_time_tol, min_capture_time, t_offset)

    hist = init_hist();
    elapsed = 0;
    reached_tol = false;

    [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
    d_now = norm(rel_lvlh - desired_rel_lvlh);

    best_X = X_ch; best_T = X_t; best_t = 0; best_d = d_now;

    X_prev2 = []; T_prev2 = []; d_prev2 = inf; t_prev2 = 0;
    X_prev = X_ch; T_prev = X_t; d_prev = d_now; t_prev = 0;

    if d_now <= pos_tol
        hist = log_full_state(hist, X_ch, X_t, t_offset);
        capture_time = 0; miss = d_now; reached_tol = true;
        return;
    end

    while elapsed < max_time - 1e-12
        dt_eff = min(dt, max_time - elapsed);
        [X_next, T_next] = rk4_step_chaser_target(X_prev, T_prev, sys, dt_eff);
        t_next = t_prev + dt_eff;
        [rel_next, ~] = relative_state_lvlh(X_next(1:6), T_next);
        d_next = norm(rel_next - desired_rel_lvlh);

        hist = log_full_state(hist, X_next, T_next, t_offset + t_next);

        if d_next < best_d
            best_X = X_next; best_T = T_next; best_t = t_next; best_d = d_next;
        end

        if d_next <= pos_tol
            [X_ref, T_ref, local_t, d_ref] = refine_distance_minimum(sys, X_prev, T_prev, desired_rel_lvlh, dt_eff, event_time_tol);
            X_ch = X_ref; X_t = T_ref; capture_time = t_prev + local_t; miss = d_ref;
            hist = trim_hist_after(hist, t_offset + capture_time);
            hist = log_full_state(hist, X_ch, X_t, t_offset + capture_time);
            [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
            reached_tol = miss <= pos_tol;
            return;
        end

        if ~isempty(X_prev2) && (t_next >= min_capture_time) && (d_prev <= d_prev2) && (d_next > d_prev)
            % A local minimum was bracketed.  Do NOT stop here unless the
            % requested tolerance is actually reached.  The external Python
            % optimizer may be targeting a later encounter, so returning at
            % the first local minimum can create km-level disagreement.
            [X_ref, T_ref, local_t, d_ref] = refine_distance_minimum(sys, X_prev2, T_prev2, desired_rel_lvlh, t_next - t_prev2, event_time_tol);
            t_ref = t_prev2 + local_t;

            if d_ref < best_d
                best_X = X_ref; best_T = T_ref; best_t = t_ref; best_d = d_ref;
            end

            if d_ref <= pos_tol
                X_ch = X_ref; X_t = T_ref; capture_time = t_ref; miss = d_ref;
                hist = trim_hist_after(hist, t_offset + capture_time);
                hist = log_full_state(hist, X_ch, X_t, t_offset + capture_time);
                [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
                reached_tol = true;
                return;
            end
            % Otherwise continue scanning until max_time and return the best
            % global sampled/refined closest approach at the end.
        end

        X_prev2 = X_prev; T_prev2 = T_prev; d_prev2 = d_prev; t_prev2 = t_prev;
        X_prev = X_next; T_prev = T_next; d_prev = d_next; t_prev = t_next;
        elapsed = t_next;
    end

    % If no local minimum was bracketed, return the best sampled state.
    X_ch = best_X; X_t = best_T; capture_time = best_t; miss = best_d;
    hist = trim_hist_after(hist, t_offset + capture_time);
    hist = log_full_state(hist, X_ch, X_t, t_offset + capture_time);
    [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
    reached_tol = miss <= pos_tol;
end

function [X_best, T_best, t_best, d_best] = refine_distance_minimum(sys, X0, T0, desired_rel_lvlh, dt_window, event_time_tol)
    phi = (sqrt(5) - 1) / 2;
    a = 0; b = dt_window;
    c = b - phi*(b-a);
    d = a + phi*(b-a);
    fc = distance_to_desired_at_time(sys, X0, T0, desired_rel_lvlh, c);
    fd = distance_to_desired_at_time(sys, X0, T0, desired_rel_lvlh, d);

    while (b - a) > event_time_tol
        if fc > fd
            a = c; c = d; fc = fd;
            d = a + phi*(b-a);
            fd = distance_to_desired_at_time(sys, X0, T0, desired_rel_lvlh, d);
        else
            b = d; d = c; fd = fc;
            c = b - phi*(b-a);
            fc = distance_to_desired_at_time(sys, X0, T0, desired_rel_lvlh, c);
        end
    end

    t_best = 0.5*(a+b);
    [X_best, T_best] = propagate_state_only(X0, T0, sys, t_best, min(1, max(t_best/20, 1e-3)));
    d_best = distance_to_desired(X_best, T_best, desired_rel_lvlh);
end

function d = distance_to_desired_at_time(sys, X0, T0, desired_rel_lvlh, t_eval)
    [X_eval, T_eval] = propagate_state_only(X0, T0, sys, t_eval, min(1, max(t_eval/20, 1e-3)));
    d = distance_to_desired(X_eval, T_eval, desired_rel_lvlh);
end

function d = distance_to_desired(X_ch, X_t, desired_rel_lvlh)
    [rel_lvlh, ~] = relative_state_lvlh(X_ch(1:6), X_t);
    d = norm(rel_lvlh - desired_rel_lvlh);
end

function [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
    execute_custom_maneuver(sys, X_st, X_target_st, dV_cmd, gamma, maneuver_opts, t_offset)

    if maneuver_opts.burn_model == "IMPULSIVE"
        [X_st, dV_mag] = apply_custom_impulse(X_st, dV_cmd, gamma, sys, maneuver_opts.direction_mode);
        burn_duration = 0;
        burn_hist = init_hist();
        burn_hist = log_full_state(burn_hist, X_st, X_target_st, t_offset);
        burn_hist = log_maneuver_stat(burn_hist, "custom_impulse", dV_mag, burn_duration);

    elseif maneuver_opts.burn_model == "FINITE_BURN"
        [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
            apply_custom_finite_burn(sys, X_st, X_target_st, dV_cmd, gamma, maneuver_opts, t_offset);

    else
        error('Unsupported custom burn model: %s', char(maneuver_opts.burn_model));
    end
end

function [X_st, dV_mag] = apply_custom_impulse(X_st, dV_cmd, gamma, sys, direction_mode)
    if nargin < 5 || isempty(direction_mode)
        direction_mode = "LOCAL_TANGENTIAL_RADIAL";
    end

    burn_dir = custom_burn_direction(X_st, gamma, direction_mode);

    X_st(4:6) = X_st(4:6) + dV_cmd * burn_dir;
    X_st(7) = X_st(7) * exp(-dV_cmd / (sys.Isp * sys.g0));
    dV_mag = dV_cmd;
end

function [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
    apply_custom_finite_burn(sys, X_st, X_target_st, dV_cmd, gamma, maneuver_opts, t_offset)

    burn_dir_eci = custom_burn_direction(X_st, gamma, maneuver_opts.direction_mode);
    [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
        apply_finite_burn_fixed_direction(sys, X_st, X_target_st, dV_cmd, burn_dir_eci, ...
                                          maneuver_opts, t_offset, "custom_finite_burn");

    if dV_cmd > 0
        fprintf('     - finite burn thrust: %.6f N, Isp %.3f s\n', ...
                maneuver_opts.finite_burn_thrust, maneuver_opts.finite_burn_isp);
    end
end

function [X_chaser_next, X_target_next] = rk4_step_chaser_target_fixed_thrust( ...
    X_chaser, X_target, sys, dt, thrust_N, Isp, burn_dir_eci)

    k1 = orbit_dynamics_fixed_thrust(X_chaser, sys, thrust_N, Isp, burn_dir_eci);
    k2 = orbit_dynamics_fixed_thrust(X_chaser + k1*dt/2, sys, thrust_N, Isp, burn_dir_eci);
    k3 = orbit_dynamics_fixed_thrust(X_chaser + k2*dt/2, sys, thrust_N, Isp, burn_dir_eci);
    k4 = orbit_dynamics_fixed_thrust(X_chaser + k3*dt, sys, thrust_N, Isp, burn_dir_eci);
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

function dX = orbit_dynamics_fixed_thrust(X, sys, thrust_N, Isp, burn_dir_eci)
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    r_norm = norm(r);

    a_g = -sys.mu / r_norm^3 * r;
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1);
                      (r(2)/r_norm)*(5*z2 - 1);
                      (r(3)/r_norm)*(5*z2 - 3) ];

    a_thrust = burn_dir_eci(:) / norm(burn_dir_eci) * (thrust_N / m);
    dm = -thrust_N / (Isp * sys.g0);

    dX = [v; a_g + a_j2 + a_thrust; dm];
end

function burn_dir = custom_burn_direction(X_st, gamma, direction_mode)
    mode = upper(string(direction_mode));

    if mode == "LOCAL_TANGENTIAL_RADIAL" || mode == "LOCAL_TANGENT_RADIAL"
        r = X_st(1:3);
        v = X_st(4:6);
        r_hat = r / norm(r);
        h_hat = cross(r, v);
        h_hat = h_hat / norm(h_hat);
        t_hat = cross(h_hat, r_hat);
        t_hat = t_hat / norm(t_hat);

        % Convention: gamma > 0 tilts tangential thrust toward radial-outward.
        burn_dir = cos(gamma) * t_hat + sin(gamma) * r_hat;
        burn_dir = burn_dir / norm(burn_dir);
    else
        error('Unsupported burn_direction_mode: %s', char(mode));
    end
end

function opts = get_maneuver_options(sys, custom_params)
    default_model = get_sys_maneuver_field(sys, 'default_burn_model', "IMPULSIVE");
    opts.burn_model = normalize_burn_model(get_param(custom_params, 'burn_model', default_model));

    default_direction = get_sys_maneuver_field(sys, 'direction_mode', "LOCAL_TANGENTIAL_RADIAL");
    opts.direction_mode = upper(string(get_param(custom_params, 'burn_direction_mode', default_direction)));

    default_thrust = get_sys_maneuver_field(sys, 'finite_burn_thrust', get_sys_field(sys, 'Thrust_Impulsive', 300.0));
    opts.finite_burn_thrust = get_param(custom_params, 'finite_burn_thrust', get_param(custom_params, 'thrust', default_thrust));

    default_isp = get_sys_maneuver_field(sys, 'finite_burn_isp', get_sys_field(sys, 'Isp', 200.0));
    opts.finite_burn_isp = get_param(custom_params, 'finite_burn_isp', get_param(custom_params, 'isp', default_isp));

    default_dt = get_sys_maneuver_field(sys, 'finite_burn_dt', 0.1);
    opts.dt_burn = get_param(custom_params, 'dt_burn', default_dt);

    opts.max_single_burn_duration = get_param(custom_params, 'max_single_burn_duration', ...
        get_sys_maneuver_field(sys, 'max_single_burn_duration', inf));
    opts.max_single_burn_delta_v = get_param(custom_params, 'max_single_burn_delta_v', ...
        get_sys_maneuver_field(sys, 'max_single_burn_delta_v', inf));
    opts.use_thrust_noise = get_bool_param(custom_params, 'use_thrust_noise', false);
end

function model = normalize_burn_model(value)
    model = upper(string(value));
    if model == "FINITE" || model == "FINITE_BURN" || model == "FINITE_IMPULSE" || model == "CONTINUOUS"
        model = "FINITE_BURN";
    elseif model == "INSTANT" || model == "INSTANTANEOUS" || model == "CUSTOM_IMPULSE"
        model = "IMPULSIVE";
    end

    if ~(model == "IMPULSIVE" || model == "FINITE_BURN")
        error('Unknown burn_model: %s. Use "IMPULSIVE" or "FINITE_BURN".', char(model));
    end
end

function value = get_sys_field(sys, name, default_value)
    if isstruct(sys) && isfield(sys, name) && ~isempty(sys.(name))
        value = sys.(name);
    else
        value = default_value;
    end
end

function value = get_sys_maneuver_field(sys, name, default_value)
    if isstruct(sys) && isfield(sys, 'maneuver') && isstruct(sys.maneuver) && ...
            isfield(sys.maneuver, name) && ~isempty(sys.maneuver.(name))
        value = sys.maneuver.(name);
    else
        value = default_value;
    end
end

function value = get_bool_param(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        raw = s.(name);
    else
        raw = default_value;
    end

    if islogical(raw)
        value = raw;
    elseif isnumeric(raw)
        value = raw ~= 0;
    else
        text = lower(strtrim(string(raw)));
        value = text == "true" || text == "on" || text == "yes" || text == "1";
    end
end

function err = phase_error_to_target(X_ch, X_t, phase_angle)
    ph = signed_phase_angle(X_ch(1:3), X_t(1:3), X_t(4:6));
    err = wrap_to_pi_custom(ph - phase_angle);
end

function ph = signed_phase_angle(r_chaser, r_target, v_target)
    % Signed angle from chaser radius vector to target radius vector, measured
    % about the target orbit angular-momentum direction. Positive means target
    % is ahead of chaser in the target orbital plane.
    h_hat = cross(r_target, v_target);
    h_hat = h_hat / norm(h_hat);
    cross_rt = cross(r_chaser, r_target);
    ph = atan2(dot(h_hat, cross_rt), dot(r_chaser, r_target));
end

function angle = get_angle_param_rad(s, name, unit_name)
    angle = get_required_scalar(s, name);
    unit = "rad";
    if isstruct(s) && isfield(s, unit_name) && ~isempty(s.(unit_name))
        unit = string(s.(unit_name));
    end
    unit = lower(unit);
    if unit == "deg" || unit == "degree" || unit == "degrees"
        angle = deg2rad(angle);
    elseif ~(unit == "rad" || unit == "radian" || unit == "radians")
        error('Unknown angle unit for %s: %s. Use "rad" or "deg".', name, unit);
    end
    angle = wrap_to_pi_custom(angle);
end

function value = get_required_scalar(s, name)
    if ~(isstruct(s) && isfield(s, name)) || isempty(s.(name))
        error('custom_params.%s must be provided before running this mode.', name);
    end
    value = s.(name);
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
        error('custom_params.%s must be a finite scalar.', name);
    end
end

function vec = get_vector_param(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        tmp = s.(name);
        vec = tmp(:);
    else
        vec = default_value(:);
    end
    if numel(vec) ~= 3 || ~isnumeric(vec) || any(~isfinite(vec))
        error('custom_params.%s must be a finite 3x1 vector.', name);
    end
end

function y = wrap_to_pi_custom(x)
    y = mod(x + pi, 2*pi) - pi;
end

function h = trim_hist_after(h, t_end)
    if isempty(h.time)
        return;
    end
    keep = h.time <= t_end + 1e-9;
    fields = {'pos','vel','mass','time','target_pos','target_vel','rel_pos','rel_pos_lvlh','rel_vel_lvlh'};
    for i = 1:numel(fields)
        f = fields{i};
        if isfield(h, f) && ~isempty(h.(f))
            h.(f) = h.(f)(:, keep);
        end
    end
    if isempty(h.time)
        h.time_end = 0;
    else
        h.time_end = h.time(end);
    end
end

function [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
    execute_delta_v_maneuver(sys, X_st, X_target_st, dv_vec_cmd, maneuver_opts, t_offset, maneuver_name)

    burn_hist = init_hist();
    dV_cmd = norm(dv_vec_cmd);

    if dV_cmd == 0
        dV_mag = 0;
        burn_duration = 0;
        burn_hist = log_full_state(burn_hist, X_st, X_target_st, t_offset);
        burn_hist = log_maneuver_stat(burn_hist, maneuver_name, dV_mag, burn_duration);
        return;
    end

    if maneuver_opts.burn_model == "IMPULSIVE"
        X_st(4:6) = X_st(4:6) + dv_vec_cmd;
        X_st(7) = X_st(7) * exp(-dV_cmd / (sys.Isp * sys.g0));
        dV_mag = dV_cmd;
        burn_duration = 0;
        burn_hist = log_full_state(burn_hist, X_st, X_target_st, t_offset);
        burn_hist = log_maneuver_stat(burn_hist, maneuver_name, dV_mag, burn_duration);

    elseif maneuver_opts.burn_model == "FINITE_BURN"
        burn_dir_eci = dv_vec_cmd / dV_cmd;
        [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
            apply_finite_burn_fixed_direction(sys, X_st, X_target_st, dV_cmd, burn_dir_eci, maneuver_opts, t_offset, maneuver_name);

    else
        error('Unsupported burn model for delta-V maneuver: %s', char(maneuver_opts.burn_model));
    end
end

function [X_st, X_target_st, dV_mag, burn_hist, burn_duration] = ...
    apply_finite_burn_fixed_direction(sys, X_st, X_target_st, dV_cmd, burn_dir_eci, maneuver_opts, t_offset, maneuver_name)

    burn_hist = init_hist();

    if dV_cmd == 0
        dV_mag = 0;
        burn_duration = 0;
        burn_hist = log_full_state(burn_hist, X_st, X_target_st, t_offset);
        burn_hist = log_maneuver_stat(burn_hist, maneuver_name, dV_mag, burn_duration);
        return;
    end

    thrust_nominal = maneuver_opts.finite_burn_thrust;
    Isp = maneuver_opts.finite_burn_isp;
    dt_burn = maneuver_opts.dt_burn;

    if thrust_nominal <= 0
        error('Finite burn requires positive thrust.');
    end
    if Isp <= 0
        error('Finite burn requires positive Isp.');
    end
    if dt_burn <= 0
        error('custom_params.dt_burn must be positive for finite burns.');
    end

    thrust_actual = thrust_nominal;
    if maneuver_opts.use_thrust_noise
        thrust_actual = thrust_nominal * (1 + randn * sys.noise.thrust_err);
    end
    if thrust_actual <= 0
        error('Finite burn actual thrust became non-positive.');
    end

    m0 = X_st(7);
    exhaust_velocity = Isp * sys.g0;
    mf_commanded = m0 * exp(-dV_cmd / exhaust_velocity);
    burn_duration = (m0 - mf_commanded) * exhaust_velocity / thrust_actual;

    elapsed = 0;
    while elapsed < burn_duration - 1e-12
        dt_eff = min(dt_burn, burn_duration - elapsed);
        [X_st, X_target_st] = rk4_step_chaser_target_fixed_thrust( ...
            X_st, X_target_st, sys, dt_eff, thrust_actual, Isp, burn_dir_eci);
        elapsed = elapsed + dt_eff;
        burn_hist = log_full_state(burn_hist, X_st, X_target_st, t_offset + elapsed);

        if X_st(7) <= 0
            error('Chaser mass depleted during finite burn.');
        end
    end

    dV_mag = exhaust_velocity * log(m0 / X_st(7));
    burn_hist = log_maneuver_stat(burn_hist, maneuver_name, dV_mag, burn_duration);
end

%% --- Helper: Hohmann Transfer Execution ---
function [X_st, X_target_st, dV_tot, sub_hist] = execute_hohmann(sys, X_st, target_r, X_target_st, custom_params, pmode)
    sub_hist = init_hist();
    dV_tot = 0;
    maneuver_opts = get_maneuver_options(sys, custom_params);

    dt_wait     = get_param(custom_params, 'dt_wait', 30);      % wait propagation step [s]
    dt_transfer = get_param(custom_params, 'dt_transfer', 10);  % transfer propagation step [s]
    dt_scan     = get_param(custom_params, 'dt_scan', 60);      % candidate spacing [s]
    max_wait    = get_param(custom_params, 'max_wait', 2*86400);% search horizon [s]
    refine_span = get_param(custom_params, 'refine_span', 180); % local refinement half-width [s]
    refine_step = get_param(custom_params, 'refine_step', 2);  % local refinement spacing [s]
    pos_tol     = get_param(custom_params, 'capture_pos_tol', 1000); % warning threshold [m]

    if pmode == 3
        max_wait = 360;
    else
        max_wait = max_wait;
    end



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

    % First Hohmann maneuver at the numerically selected epoch.
    dv1_vec = hohmann_departure_delta_v_vec(X_st, target_r, sys);
    [X_st, X_target_st, dV1, burn_hist, ~] = execute_delta_v_maneuver( ...
        sys, X_st, X_target_st, dv1_vec, maneuver_opts, sub_hist.time_end, "hohmann_departure");
    dV_tot = dV_tot + dV1;
    sub_hist = append_hist(sub_hist, burn_hist);

    [X_st, X_target_st, transfer_hist] = propagate_for_duration(X_st, X_target_st, sys, TOF, dt_transfer, sub_hist.time_end);
    sub_hist = append_hist(sub_hist, transfer_hist);

    % Circularization maneuver at arrival radius.
    dv2_vec = circularization_delta_v_vec(X_st, X_target_st, sys);
    [X_st, X_target_st, dV2, burn_hist, ~] = execute_delta_v_maneuver( ...
        sys, X_st, X_target_st, dv2_vec, maneuver_opts, sub_hist.time_end, "hohmann_circularization");
    dV_tot = dV_tot + dV2;
    sub_hist = append_hist(sub_hist, burn_hist);

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

%% --- Helper: Multi-Hohmann Transfer Execution ---
function [X_st, X_target_st, dV_tot, sub_hist] = execute_multi_hohmann(sys, X_st, target_r, X_target_st, custom_params, pmode)
    sub_hist = init_hist();
    dV_tot = 0;
    maneuver_opts = get_maneuver_options(sys, custom_params);

    dt_wait     = get_param(custom_params, 'dt_wait', 30);
    dt_transfer = get_param(custom_params, 'dt_transfer', 10);
    dt_scan     = get_param(custom_params, 'dt_scan', 60);
    max_wait    = get_param(custom_params, 'max_wait', 2*86400);
    refine_span = get_param(custom_params, 'refine_span', 180);
    refine_step = get_param(custom_params, 'refine_step', 2);
    pos_tol     = get_param(custom_params, 'capture_pos_tol', 1000);

    if pmode == 3
        max_wait = 360;
    end

    if isfield(sys, 'capture') && isfield(sys.capture, 'r_rel0')
        desired_rel_lvlh = sys.capture.r_rel0(:);
    else
        desired_rel_lvlh = [-5000; 0; 0];
    end

    leg_count = get_param(custom_params, 'multi_hohmann_legs', []);
    if isempty(leg_count)
        leg_count = choose_multi_hohmann_leg_count(sys, X_st, target_r, maneuver_opts);
    end
    leg_count = max(1, round(leg_count));

    radii = linspace(norm(X_st(1:3)), target_r, leg_count + 1);
    leg_targets = radii(2:end);
    total_tof = multi_hohmann_total_tof(sys, radii);

    fprintf('   * Multi-Hohmann mode started: %d leg(s), burn model %s\n', ...
            leg_count, char(maneuver_opts.burn_model));
    fprintf('     - thermal limits: max burn %.3f s, max single delta-V %.6f m/s\n', ...
            maneuver_opts.max_single_burn_duration, maneuver_opts.max_single_burn_delta_v);

    if ~isempty(X_target_st)
        fprintf('   * J2-aware Multi-Hohmann phase search started...\n');
        best_wait = find_best_wait_time_multi(sys, X_st, X_target_st, radii, desired_rel_lvlh, ...
                                              dt_scan, max_wait, refine_span, refine_step, dt_transfer, pmode);
        fprintf('     - selected wait time: %.1f s (total transfer TOF %.1f s)\n', best_wait, total_tof);

        [X_st, X_target_st, wait_hist] = propagate_for_duration(X_st, X_target_st, sys, best_wait, dt_wait, 0);
        sub_hist = append_hist(sub_hist, wait_hist);
    end

    for ii = 1:leg_count
        next_r = leg_targets(ii);
        r_current = norm(X_st(1:3));
        tof_leg = pi * sqrt(((r_current + next_r)/2)^3 / sys.mu);

        dv_depart = hohmann_departure_delta_v_vec(X_st, next_r, sys);
        [X_st, X_target_st, dV_leg, burn_hist, burn_duration] = execute_delta_v_maneuver( ...
            sys, X_st, X_target_st, dv_depart, maneuver_opts, sub_hist.time_end, ...
            sprintf("multi_hohmann_%02d_departure", ii));
        dV_tot = dV_tot + dV_leg;
        sub_hist = append_hist(sub_hist, burn_hist);
        fprintf('     - leg %02d departure: dV %.6f m/s, burn %.3f s\n', ii, dV_leg, burn_duration);

        [X_st, X_target_st, transfer_hist] = propagate_for_duration(X_st, X_target_st, sys, tof_leg, dt_transfer, sub_hist.time_end);
        sub_hist = append_hist(sub_hist, transfer_hist);

        dv_circ = circularization_delta_v_vec(X_st, X_target_st, sys);
        [X_st, X_target_st, dV_leg, burn_hist, burn_duration] = execute_delta_v_maneuver( ...
            sys, X_st, X_target_st, dv_circ, maneuver_opts, sub_hist.time_end, ...
            sprintf("multi_hohmann_%02d_circularization", ii));
        dV_tot = dV_tot + dV_leg;
        sub_hist = append_hist(sub_hist, burn_hist);
        fprintf('     - leg %02d circularization: dV %.6f m/s, burn %.3f s\n', ii, dV_leg, burn_duration);
    end

    if ~isempty(X_target_st)
        [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_st(1:6), X_target_st);
        miss = rel_lvlh - desired_rel_lvlh;
        fprintf('     - capture LVLH position: [%+.1f, %+.1f, %+.1f] m\n', rel_lvlh(1), rel_lvlh(2), rel_lvlh(3));
        fprintf('     - capture error from desired: %.1f m, rel-speed: %.4f m/s\n', norm(miss), norm(rel_vel_lvlh));
        if norm(miss) > pos_tol
            fprintf('     - warning: capture error is larger than %.0f m. Adjust leg count or search settings.\n', pos_tol);
        end
    end
end

function leg_count = choose_multi_hohmann_leg_count(sys, X_st, target_r, maneuver_opts)
    max_legs = 50;
    leg_count = 1;

    for n = 1:max_legs
        radii = linspace(norm(X_st(1:3)), target_r, n + 1);
        [max_dv, max_duration] = estimate_multi_hohmann_max_burn(sys, X_st(7), radii, maneuver_opts);

        dv_ok = isinf(maneuver_opts.max_single_burn_delta_v) || max_dv <= maneuver_opts.max_single_burn_delta_v;
        duration_ok = isinf(maneuver_opts.max_single_burn_duration) || max_duration <= maneuver_opts.max_single_burn_duration;

        if dv_ok && duration_ok
            leg_count = n;
            return;
        end
    end

    leg_count = max_legs;
end

function [max_dv, max_duration] = estimate_multi_hohmann_max_burn(sys, mass0, radii, maneuver_opts)
    max_dv = 0;
    max_duration = 0;
    mass = mass0;

    for ii = 1:(numel(radii)-1)
        r1 = radii(ii);
        r2 = radii(ii+1);
        a = 0.5 * (r1 + r2);
        v1 = sqrt(sys.mu / r1);
        v2 = sqrt(sys.mu / r2);
        vt1 = sqrt(sys.mu * (2/r1 - 1/a));
        vt2 = sqrt(sys.mu * (2/r2 - 1/a));
        dvs = [abs(vt1 - v1), abs(v2 - vt2)];

        for jj = 1:numel(dvs)
            dv = dvs(jj);
            duration = estimate_finite_burn_duration(sys, mass, dv, maneuver_opts);
            max_dv = max(max_dv, dv);
            max_duration = max(max_duration, duration);
            mass = mass * exp(-dv / (maneuver_opts.finite_burn_isp * sys.g0));
        end
    end
end

function duration = estimate_finite_burn_duration(sys, mass0, delta_v, maneuver_opts)
    if maneuver_opts.burn_model == "IMPULSIVE" || delta_v == 0
        duration = 0;
        return;
    end

    ve = maneuver_opts.finite_burn_isp * sys.g0;
    mf = mass0 * exp(-delta_v / ve);
    duration = (mass0 - mf) * ve / maneuver_opts.finite_burn_thrust;
end

function total_tof = multi_hohmann_total_tof(sys, radii)
    total_tof = 0;
    for ii = 1:(numel(radii)-1)
        a = 0.5 * (radii(ii) + radii(ii+1));
        total_tof = total_tof + pi * sqrt(a^3 / sys.mu);
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
        if abs(acos(dot(X_ch_wait(1:3), X_t_wait(1:3))/(norm(X_ch_wait(1:3))*norm(X_t_wait(1:3)))) - sys.phase) > 0.5 * sys.phase && pmode == 1
            continue
        end

        [rel_lvlh, rel_vel_lvlh] = predict_hohmann_capture(sys, X_ch_wait, X_t_wait, target_r, TOF, dt_transfer);
        pos_error = norm(rel_lvlh - desired_rel_lvlh);
        metric = pos_error;

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

function best_wait = find_best_wait_time_multi(sys, X_ch0, X_t0, radii, desired_rel_lvlh, dt_scan, max_wait, refine_span, refine_step, dt_transfer, pmode)
    candidate_times = 0:dt_scan:max_wait;
    [best_wait, ~] = scan_wait_candidates_multi(sys, X_ch0, X_t0, radii, desired_rel_lvlh, candidate_times, dt_transfer, pmode);

    t1 = max(0, best_wait - refine_span);
    t2 = min(max_wait, best_wait + refine_span);
    refined_times = t1:refine_step:t2;
    [best_wait, ~] = scan_wait_candidates_multi(sys, X_ch0, X_t0, radii, desired_rel_lvlh, refined_times, dt_transfer, pmode);
end

function [best_wait, best_metric] = scan_wait_candidates_multi(sys, X_ch0, X_t0, radii, desired_rel_lvlh, candidate_times, dt_transfer, pmode)
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

        if abs(acos(dot(X_ch_wait(1:3), X_t_wait(1:3))/(norm(X_ch_wait(1:3))*norm(X_t_wait(1:3)))) - sys.phase) > 0.5 * sys.phase && pmode == 1
            continue
        end

        [rel_lvlh, ~] = predict_multi_hohmann_capture(sys, X_ch_wait, X_t_wait, radii, dt_transfer);
        metric = norm(rel_lvlh - desired_rel_lvlh);

        if metric < best_metric
            best_metric = metric;
            best_wait = t_now;
        end
    end
end

function [rel_lvlh, rel_vel_lvlh] = predict_multi_hohmann_capture(sys, X_ch_wait, X_t_wait, radii, dt_transfer)
    X_ch = X_ch_wait;
    X_t = X_t_wait;

    for ii = 2:numel(radii)
        target_r = radii(ii);
        r_current = norm(X_ch(1:3));
        tof_leg = pi * sqrt(((r_current + target_r)/2)^3 / sys.mu);
        [X_ch, ~] = apply_hohmann_departure_impulse(X_ch, target_r, X_t, sys, false);
        [X_ch, X_t] = propagate_state_only(X_ch, X_t, sys, tof_leg, dt_transfer);
        [X_ch, ~] = apply_circularization_impulse(X_ch, target_r, X_t, sys, false);
    end

    [rel_lvlh, rel_vel_lvlh] = relative_state_lvlh(X_ch(1:6), X_t);
end

%% --- Helper: Impulses ---
function dv_vec = hohmann_departure_delta_v_vec(X_st, target_r, sys)
    r_current = norm(X_st(1:3));
    v_current = norm(X_st(4:6));
    a_trans = (r_current + target_r) / 2;
%    v_trans1 = sqrt(sys.mu * (2/r_current - 1/a_trans));
    v_trans1 = sqrt(sys.mu * (2/r_current - 1/a_trans - sys.J2 * sys.Re^2 / r_current^3 * (3 * (X_st(3)/r_current)^2 - 1)));

    dV_cmd = v_trans1 - v_current;
    dv_vec = dV_cmd * X_st(4:6) / v_current;
end

function dv_vec_cmd = circularization_delta_v_vec(X_st, X_target_st, sys)
    r = X_st(1:3);
    v = X_st(4:6);
    r_norm = norm(r);
    h_vec = cross(r, v);
    tangential_dir = cross(h_vec, r);
    tangential_dir = tangential_dir / norm(tangential_dir);
    v_circ_req = tangential_dir * sqrt(sys.mu / r_norm);
    dv_vec_cmd = v_circ_req - v;
end

function [X_st, dV_mag] = apply_hohmann_departure_impulse(X_st, target_r, ~, sys, use_noise)
    if nargin < 4, use_noise = true; end
    dv_vec = hohmann_departure_delta_v_vec(X_st, target_r, sys);
    dV_cmd = norm(dv_vec);
    dV_mag = abs(dV_cmd);
    if use_noise
        dV_mag = dV_mag * (1 + randn * sys.noise.thrust_err);
    end
    if dV_cmd > 0
        X_st(4:6) = X_st(4:6) + dv_vec/dV_cmd*dV_mag;
    end
    X_st(7) = X_st(7) * exp(-dV_mag / (sys.Isp * sys.g0));
end

function [X_st, dV_mag] = apply_circularization_impulse(X_st, target_r, X_target_st, sys, use_noise)
    if nargin < 3, use_noise = true; end
    dv_vec_cmd = circularization_delta_v_vec(X_st, X_target_st, sys);
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
function [X_chaser, X_target, hist] = propagate_for_duration(X_chaser, X_target, sys, duration, dt, t_offset)
    if nargin < 6, t_offset = 0; end
    hist = init_hist();
    elapsed = 0;
    while elapsed < duration - 1e-12
        dt_eff = min(dt, duration - elapsed);
        [X_chaser, X_target] = rk4_step_chaser_target(X_chaser, X_target, sys, dt_eff);
        elapsed = elapsed + dt_eff;
        hist = log_full_state(hist, X_chaser, X_target, t_offset + elapsed);
    end
end

function [X_chaser, X_target] = propagate_state_only(X_chaser, X_target, sys, duration, dt)
    if duration <= 0
        return;
    end
    elapsed = 0;
    while elapsed < duration - 1e-12
        dt_eff = min(dt, duration - elapsed);
        [X_chaser, X_target] = rk4_step_chaser_target(X_chaser, X_target, sys, dt_eff);
        elapsed = elapsed + dt_eff;
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
    hist.maneuver_delta_v = [];
    hist.maneuver_duration = [];
    hist.maneuver_name = strings(1,0);
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

function hist = log_maneuver_stat(hist, name, delta_v, duration)
    if ~isfield(hist, 'maneuver_delta_v')
        hist.maneuver_delta_v = [];
    end
    if ~isfield(hist, 'maneuver_duration')
        hist.maneuver_duration = [];
    end
    if ~isfield(hist, 'maneuver_name')
        hist.maneuver_name = strings(1,0);
    end

    hist.maneuver_delta_v = [hist.maneuver_delta_v, delta_v];
    hist.maneuver_duration = [hist.maneuver_duration, duration];
    hist.maneuver_name(end+1) = string(name);
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
    stat_fields = {'maneuver_delta_v','maneuver_duration','maneuver_name'};
    for i = 1:numel(stat_fields)
        f = stat_fields{i};
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
