function results = Validate_Reentry_Core_Equivalence()
    % Regression checks for the extracted reusable re-entry physics kernel.
    %
    % The local reference equations below freeze the pre-extraction
    % Reentry_Propagator implementation. The end-to-end constants were also
    % recorded before the propagator was rewired to +reentry_core. These are
    % numerical-equivalence checks, not validation of the assumed physics.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(project_root);

    sys = Mission_Config();
    % Force the repository-owned atmosphere implementation so equivalence
    % does not depend on Aerospace Toolbox availability or release details.
    sys.environment.atmospheric_drag.use_matlab_atmosisa = false;
    heat_k = sys.reentry_vehicle.sutton_graves_k;

    X14 = synthetic_entry_state_local(sys);
    X7 = [X14(1:6); X14(14)];
    spaceplane = spaceplane_shape_local(sys);

    aux_core = reentry_core.evaluate_state( ...
        X7, sys, spaceplane, 20, 25, true, heat_k);
    aux_reference = legacy_aux_reference_local( ...
        X7, sys, spaceplane, 20, 25, true, heat_k);
    assert_aux_equal_local(aux_core, aux_reference, 'SPACEPLANE state evaluation');

    dX_core = reentry_core.dynamics( ...
        X7, sys, spaceplane, 20, 25, true, heat_k);
    dX_reference = legacy_dynamics_reference_local( ...
        X7, sys, spaceplane, 20, 25, true, heat_k);
    assert_close_local(dX_core, dX_reference, 1e-13, 1e-13, ...
        'SPACEPLANE state derivative');

    X_step_core = reentry_core.rk4_step( ...
        X7, sys, spaceplane, 20, 25, true, heat_k, 0.5);
    X_step_reference = legacy_rk4_reference_local( ...
        X7, sys, spaceplane, 20, 25, true, heat_k, 0.5);
    assert_close_local(X_step_core, X_step_reference, 1e-12, 1e-13, ...
        'SPACEPLANE RK4 step');

    j2_sys = sys;
    j2_sys.environment.atmospheric_drag.co_rotate_atmosphere = false;
    j2_shape = spaceplane;
    j2_shape.gravity_model = "J2";
    aux_j2_core = reentry_core.evaluate_state( ...
        X7, j2_sys, j2_shape, 23, -40, false, heat_k);
    aux_j2_reference = legacy_aux_reference_local( ...
        X7, j2_sys, j2_shape, 23, -40, false, heat_k);
    assert_aux_equal_local(aux_j2_core, aux_j2_reference, ...
        'J2/nonrotating-atmosphere state evaluation');

    capsule = capsule_shape_local(sys);
    X_capsule = X7;
    X_capsule(7) = capsule.mass_kg;
    aux_capsule_core = reentry_core.evaluate_state( ...
        X_capsule, sys, capsule, capsule.default_aoa_deg, 0, true, heat_k);
    aux_capsule_reference = legacy_aux_reference_local( ...
        X_capsule, sys, capsule, capsule.default_aoa_deg, 0, true, heat_k);
    assert_aux_equal_local(aux_capsule_core, aux_capsule_reference, ...
        'CAPSULE state evaluation');

    target_altitude = 100e3;
    [X_cross_core, t_cross_core] = reentry_core.refine_altitude_crossing( ...
        X7, sys, spaceplane, 20, 0, false, heat_k, target_altitude, 30);
    [X_cross_reference, t_cross_reference] = ...
        legacy_refine_altitude_reference_local( ...
            X7, sys, spaceplane, 20, 0, false, heat_k, ...
            target_altitude, 30);
    assert_close_local(X_cross_core, X_cross_reference, 1e-11, 1e-13, ...
        'Altitude-event crossing state');
    assert_close_local(t_cross_core, t_cross_reference, 1e-12, 0, ...
        'Altitude-event crossing time');

    [X_at_initial_altitude, t_at_initial_altitude] = ...
        reentry_core.refine_altitude_crossing( ...
            X7, sys, spaceplane, 20, 0, false, heat_k, ...
            norm(X7(1:3))-sys.Re, 1);
    assert_close_local(X_at_initial_altitude, X7, 0, 0, ...
        'Initial altitude event state');
    assert(t_at_initial_altitude == 0, ...
        'An event at the initial altitude must return t=0.');
    assert_error_id_local(@() reentry_core.refine_altitude_crossing( ...
        X7, sys, spaceplane, 20, 0, false, heat_k, -1e6, 1), ...
        'reentry_core:refineAltitude:UnbracketedCrossing');
    assert_error_id_local(@() reentry_core.refine_speed_crossing( ...
        X7, sys, spaceplane, 20, 0, false, heat_k, 0, 1), ...
        'reentry_core:refineSpeed:UnbracketedCrossing');

    relay = synthetic_relay_state_local(sys);
    [clear_core, clearance_core, elevation_core] = ...
        reentry_core.line_of_sight_geometry( ...
            relay(1:3), X14(1:3), sys, 0);
    [clear_reference, clearance_reference, elevation_reference] = ...
        legacy_los_reference_local(relay(1:3), X14(1:3), sys, 0);
    assert(clear_core == clear_reference, 'LOS Boolean changed during extraction.');
    assert_close_local([clearance_core, elevation_core], ...
        [clearance_reference, elevation_reference], 1e-12, 1e-13, ...
        'LOS geometry');
    raap_core = reentry_core.raap_deg( ...
        X14(1:3), aux_core.v_rel, relay(1:3), ...
        aux_core.aoa_deg, aux_core.bank_angle_deg, [-1;0;0]);
    raap_reference = legacy_raap_reference_local( ...
        X14(1:3), aux_reference.v_rel, relay(1:3), ...
        aux_reference.aoa_deg, aux_reference.bank_angle_deg, [-1;0;0]);
    assert_close_local(raap_core, raap_reference, 1e-12, 1e-13, ...
        'RAAP geometry');

    % Frozen end-to-end values captured from the pre-extraction propagator.
    base = struct('shape_name', "COMPROMISE", 'max_time', 300, ...
        'terminal_altitude', 100e3, 'lift_enabled', false, 'dt', 2.0);
    [X_space_final, ~, hist_space, summary_space] = ...
        Reentry_Propagator(sys, X14, [], 0, base);
    expected_space_state = [ ...
        6475356.70756091; 189774.867421921; 0; ...
        -1046.94527424721; 7753.68099883433; 0; 1000];
    expected_space_metrics = [ ...
        14; 24.467853307724; -2.39228829741478e-05; ...
        15.1552932003681; 234502.560490548; 2747553.99771515; ...
        0.00294056181276533; 40; 0.961602319629914; 0; 0];
    actual_space_metrics = [ ...
        numel(hist_space.time); summary_space.duration_s; ...
        summary_space.terminal_altitude_error_m; ...
        summary_space.max_dynamic_pressure_Pa; ...
        summary_space.max_heat_flux_W_m2; ...
        summary_space.total_heat_load_J_m2; summary_space.max_g_load; ...
        hist_space.aoa_deg(1); hist_space.cd(1); ...
        hist_space.cl(1); hist_space.ld(1)];
    assert_close_local(X_space_final([1:6,14]), expected_space_state, ...
        1e-8, 1e-12, 'Frozen SPACEPLANE terminal state');
    assert_close_local(actual_space_metrics, expected_space_metrics, ...
        1e-8, 1e-12, 'Frozen SPACEPLANE metrics');

    capsule_case = struct('vehicle_mode', "CAPSULE", 'max_time', 20, ...
        'terminal_altitude', 100e3, 'lift_enabled', true, 'dt', 1.0);
    [X_capsule_final, ~, hist_capsule, summary_capsule] = ...
        Reentry_Propagator(sys, X14, [], 0, capsule_case);
    expected_capsule_state = [ ...
        6479940.98556156; 155126.57295919; 0; ...
        -1004.33087734998; 7754.38402485689; 0; 60];
    expected_capsule_metrics = [ ...
        21; 20; 103660.546205724; 7.96622277745646; ...
        60877.1606226645; 662709.921393009; ...
        0.0100507416050252; 0; 1.3; 0.325; 0.25];
    actual_capsule_metrics = [ ...
        numel(hist_capsule.time); summary_capsule.duration_s; ...
        summary_capsule.final_altitude_m; ...
        summary_capsule.max_dynamic_pressure_Pa; ...
        summary_capsule.max_heat_flux_W_m2; ...
        summary_capsule.total_heat_load_J_m2; summary_capsule.max_g_load; ...
        hist_capsule.aoa_deg(1); hist_capsule.cd(1); ...
        hist_capsule.cl(1); hist_capsule.ld(1)];
    assert_close_local(X_capsule_final([1:6,14]), expected_capsule_state, ...
        1e-8, 1e-12, 'Frozen CAPSULE terminal state');
    assert_close_local(actual_capsule_metrics, expected_capsule_metrics, ...
        1e-8, 1e-12, 'Frozen CAPSULE metrics');

    results = struct();
    results.passed = true;
    results.spaceplane_terminal_state_max_abs_error = ...
        max(abs(X_space_final([1:6,14]) - expected_space_state));
    results.capsule_terminal_state_max_abs_error = ...
        max(abs(X_capsule_final([1:6,14]) - expected_capsule_state));
    results.altitude_event_time_error_s = ...
        abs(t_cross_core - t_cross_reference);
    results.raap_error_deg = abs(raap_core - raap_reference);

    fprintf('Reentry core equivalence validation passed.\n');
    fprintf('  SPACEPLANE frozen-state max error: %.6g\n', ...
        results.spaceplane_terminal_state_max_abs_error);
    fprintf('  CAPSULE frozen-state max error   : %.6g\n', ...
        results.capsule_terminal_state_max_abs_error);
    fprintf('  altitude-event time error        : %.6g s\n', ...
        results.altitude_event_time_error_s);
    fprintf('  RAAP reference error             : %.6g deg\n', ...
        results.raap_error_deg);
end

function shape = spaceplane_shape_local(sys)
    shape = sys.reentry_vehicle.shapes.COMPROMISE;
    cfg = sys.reentry_vehicle.spaceplane;
    shape.vehicle_mode = "SPACEPLANE";
    shape.gravity_model = "CENTRAL_SPHERICAL";
    shape.aero_model = cfg.aero_model;
    shape.aoa_profile_mode = cfg.aoa_profile_mode;
    shape.aoa_speed_grid_m_s = cfg.aoa_speed_grid_m_s;
    shape.aoa_values_deg = cfg.aoa_values_deg;
    shape.cl_polynomial = cfg.cl_polynomial;
    shape.cd_from_cl_polynomial = cfg.cd_from_cl_polynomial;
    shape.density_scale = sys.reentry_vehicle.uncertainty.density_scale;
    shape.cd_scale = sys.reentry_vehicle.uncertainty.cd_scale;
    shape.ld_scale = sys.reentry_vehicle.uncertainty.ld_scale;
end

function shape = capsule_shape_local(sys)
    shape = sys.reentry_vehicle.shapes.CAPSULE;
    cfg = sys.reentry_vehicle.capsule;
    shape.vehicle_mode = "CAPSULE";
    shape.gravity_model = "CENTRAL_SPHERICAL";
    shape.aero_model = cfg.aero_model;
    shape.aoa_profile_mode = "CONSTANT_TRIM";
    shape.default_aoa_deg = cfg.trim_aoa_deg;
    shape.nominal_ld = cfg.nominal_ld;
    shape.mass_kg = cfg.mass_kg;
    shape.guidance_activation_drag_g = cfg.guidance_activation_drag_g;
    shape.guidance_cutoff_mach = cfg.guidance_cutoff_mach;
    shape.density_scale = sys.reentry_vehicle.uncertainty.density_scale;
    shape.cd_scale = sys.reentry_vehicle.uncertainty.cd_scale;
    shape.ld_scale = sys.reentry_vehicle.uncertainty.ld_scale;
end

function aux = legacy_aux_reference_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    r = X(1:3);
    v = X(4:6);
    mass = X(7);
    altitude = norm(r) - sys.Re;

    atmosphere_opts = sys.environment.atmospheric_drag;
    atmosphere_opts.enabled = true;
    atmosphere_opts.model = "ISA76";
    [rho, ~, ~, a_sound] = ...
        Standard_Atmosphere_Density(altitude, atmosphere_opts);
    rho = rho * field_or_local(shape, 'density_scale', 1.0);

    if bool_or_local(atmosphere_opts, 'co_rotate_atmosphere', true)
        omega = field_or_local(atmosphere_opts, ...
            'earth_rotation_rad_s', 7.2921159e-5);
        if isscalar(omega)
            omega_vec = [0;0;omega];
        else
            omega_vec = omega(:);
        end
        v_atm = cross(omega_vec, r);
    else
        v_atm = zeros(3,1);
    end
    v_rel = v - v_atm;
    speed_rel = norm(v_rel);
    mach = speed_rel / max(a_sound, eps);
    [aoa_actual_deg, bank_actual_deg] = ...
        legacy_resolve_commands_local( ...
            shape, aoa_deg, bank_angle_deg, speed_rel, mach);

    a_gravity = legacy_gravity_local(r, sys, shape);
    a_aero = zeros(3,1);
    a_drag = zeros(3,1);
    a_lift = zeros(3,1);
    ld = 0;
    cd = NaN;
    cl = 0;
    dynamic_pressure = 0;
    heat_flux = 0;

    if speed_rel > 0 && rho > 0 && mass > 0
        v_rel_hat = v_rel / speed_rel;
        dynamic_pressure = 0.5 * rho * speed_rel^2;
        [cd, cl, ld] = legacy_aero_local(shape, aoa_actual_deg, mach);
        area = field_or_local(shape, 'reference_area_m2', 1.0);
        drag_accel_mag = dynamic_pressure * cd * area / mass;
        a_drag = -drag_accel_mag * v_rel_hat;

        if lift_enabled
            lift_dir = legacy_lift_direction_local( ...
                r, v_rel, bank_actual_deg);
            a_lift = ld * drag_accel_mag * lift_dir;
        else
            cl = 0;
            ld = 0;
        end

        nose_radius = max(field_or_local( ...
            shape, 'nose_radius_m', 0.05), 1e-4);
        heat_flux = heat_k * sqrt(rho / nose_radius) * speed_rel^3;
        a_aero = a_drag + a_lift;
    end

    aux = struct();
    aux.altitude = altitude;
    aux.rho = rho;
    aux.v_rel = v_rel;
    aux.speed_rel = speed_rel;
    aux.fpa_deg = legacy_fpa_local(r, v);
    aux.fpa_rel_deg = legacy_fpa_local(r, v_rel);
    aux.mach = mach;
    aux.dynamic_pressure = dynamic_pressure;
    aux.heat_flux = heat_flux;
    aux.g_load = norm(a_aero) / sys.g0;
    aux.drag_g = norm(a_drag) / sys.g0;
    aux.aoa_deg = aoa_actual_deg;
    aux.bank_angle_deg = bank_actual_deg;
    aux.cd = cd;
    aux.cl = cl;
    aux.ld = ld;
    guidance_threshold = field_or_local( ...
        shape, 'guidance_activation_drag_g', inf);
    aux.capsule_guidance_active = aux.drag_g >= guidance_threshold;
    guidance_cutoff_mach = field_or_local( ...
        shape, 'guidance_cutoff_mach', -inf);
    aux.capsule_guidance_window_active = ...
        aux.capsule_guidance_active && aux.mach >= guidance_cutoff_mach;
    aux.a_gravity = a_gravity;
    aux.a_aero = a_aero;
    aux.a_drag = a_drag;
    aux.a_lift = a_lift;
end

function dX = legacy_dynamics_reference_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    aux = legacy_aux_reference_local( ...
        X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    dX = [X(4:6); aux.a_gravity + aux.a_aero; 0];
end

function X_next = legacy_rk4_reference_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, dt)
    f = @(X_now) legacy_dynamics_reference_local( ...
        X_now, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    k1 = f(X);
    k2 = f(X + 0.5 * dt * k1);
    k3 = f(X + 0.5 * dt * k2);
    k4 = f(X + dt * k3);
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function [X_cross, t_cross] = legacy_refine_altitude_reference_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, target_altitude, dt_window)
    lo = 0;
    hi = dt_window;
    X_cross = legacy_rk4_reference_local( ...
        X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, hi);

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        X_mid = legacy_rk4_reference_local( ...
            X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, mid);
        altitude_mid = norm(X_mid(1:3)) - sys.Re;
        if altitude_mid > target_altitude
            lo = mid;
        else
            hi = mid;
            X_cross = X_mid;
        end
        if hi - lo <= 1e-7
            break;
        end
    end
    t_cross = hi;
end

function [aoa_deg, bank_deg] = legacy_resolve_commands_local(shape, aoa_fallback_deg, bank_fallback_deg, speed_rel, mach)
    profile_mode = upper(string(field_or_local( ...
        shape, 'aoa_profile_mode', "CONSTANT")));
    bank_deg = bank_fallback_deg;
    if profile_mode == "SPEED_SCHEDULE"
        speed_grid = field_or_local(shape, ...
            'aoa_speed_grid_m_s', [0 8000]);
        aoa_grid = field_or_local(shape, ...
            'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        speed_query = min(max(speed_rel, speed_grid(1)), speed_grid(end));
        aoa_deg = interp1(speed_grid, aoa_grid, speed_query, 'linear');
    elseif profile_mode == "MACH_SCHEDULE"
        mach_grid = field_or_local(shape, 'aoa_mach_grid', [0 30]);
        aoa_grid = field_or_local(shape, ...
            'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        mach_query = min(max(mach, mach_grid(1)), mach_grid(end));
        aoa_deg = interp1(mach_grid, aoa_grid, mach_query, 'linear');
    else
        aoa_deg = aoa_fallback_deg;
    end
end

function [cd, cl, ld] = legacy_aero_local(shape, aoa_deg, ~)
    aero_model = upper(string(field_or_local( ...
        shape, 'aero_model', "LEGACY_SHAPE_LD")));
    if aero_model == "PAPER_RLV_POLYNOMIAL"
        cl_coeff = field_or_local(shape, 'cl_polynomial', ...
            [-0.041065, 0.016292, 0.0002602]);
        cd_coeff = field_or_local(shape, 'cd_from_cl_polynomial', ...
            [0.080505, -0.03026, 0.86495]);
        cl = cl_coeff(1) + cl_coeff(2)*aoa_deg + cl_coeff(3)*aoa_deg^2;
        cd = cd_coeff(1) + cd_coeff(2)*cl + cd_coeff(3)*cl^2;
        cd = max(cd, 1e-6);
        ld = cl / cd;
    elseif aero_model == "PAPER_CAPSULE_REDUCED"
        cd = field_or_local(shape, 'cd', 1.3);
        ld = field_or_local(shape, 'nominal_ld', 0.25);
    else
        cd = field_or_local(shape, 'cd', 1.2);
        aoa_grid = field_or_local(shape, 'ld_aoa_deg', [0 20 80]);
        ld_grid = field_or_local(shape, 'ld_values', [0 1 0]);
        aoa = min(max(aoa_deg, min(aoa_grid)), max(aoa_grid));
        ld = max(0, interp1(aoa_grid, ld_grid, aoa, 'pchip'));
    end
    cd = cd * field_or_local(shape, 'cd_scale', 1.0);
    ld = ld * field_or_local(shape, 'ld_scale', 1.0);
    cl = cd * ld;
end

function a = legacy_gravity_local(r, sys, shape)
    r_norm = norm(r);
    a_g = -sys.mu / r_norm^3 * r;
    gravity_model = upper(string(field_or_local( ...
        shape, 'gravity_model', "CENTRAL_SPHERICAL")));
    if gravity_model == "CENTRAL_SPHERICAL" || gravity_model == "CENTRAL"
        a = a_g;
        return;
    end
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ ...
        (r(1)/r_norm)*(5*z2 - 1); ...
        (r(2)/r_norm)*(5*z2 - 1); ...
        (r(3)/r_norm)*(5*z2 - 3) ];
    a = a_g + a_j2;
end

function dir = legacy_lift_direction_local(r, v_rel, bank_angle_deg)
    dir = zeros(3,1);
    h = cross(r, v_rel);
    if norm(h) <= eps || norm(v_rel) <= eps
        return;
    end
    vhat = v_rel / norm(v_rel);
    hhat = h / norm(h);
    dir = cross(vhat, hhat);
    if norm(dir) <= eps
        dir = zeros(3,1);
        return;
    end
    dir = dir / norm(dir);
    bank = deg2rad(bank_angle_deg);
    dir = dir*cos(bank) + cross(vhat, dir)*sin(bank) + ...
        vhat*dot(vhat, dir)*(1-cos(bank));
    dir = dir / norm(dir);
end

function fpa_deg = legacy_fpa_local(r, v)
    r_hat = r / norm(r);
    v_radial = dot(v, r_hat);
    v_horizontal = norm(v - v_radial * r_hat);
    fpa_deg = rad2deg(atan2(v_radial, v_horizontal));
end

function [clear, clearance, elevation_deg] = legacy_los_reference_local(r_relay, r_rv, sys, margin_altitude)
    los = r_relay - r_rv;
    los_norm = norm(los);
    if los_norm <= eps
        clearance = norm(r_rv) - sys.Re;
        clear = clearance > margin_altitude;
        elevation_deg = 90;
        return;
    end
    u = r_relay - r_rv;
    tau = -dot(r_rv, u) / dot(u, u);
    tau = min(1, max(0, tau));
    closest = r_rv + tau * u;
    clearance = norm(closest) - sys.Re;
    clear = clearance > margin_altitude;
    up = r_rv / norm(r_rv);
    elevation_deg = rad2deg(asin(dot(los / los_norm, up)));
end

function raap_deg = legacy_raap_reference_local(r, v_rel, r_relay, aoa_deg, bank_deg, boresight_body)
    raap_deg = NaN;
    if norm(v_rel) <= eps || norm(r_relay - r) <= eps
        return;
    end
    x_t = v_rel / norm(v_rel);
    r_hat = r / norm(r);
    y_t = r_hat - dot(r_hat, x_t) * x_t;
    if norm(y_t) <= eps
        return;
    end
    y_t = y_t / norm(y_t);
    z_t = cross(x_t, y_t);
    z_t = z_t / norm(z_t);
    w_t = [x_t, y_t, z_t]' * (r_relay - r);
    alpha = deg2rad(aoa_deg);
    sigma = deg2rad(bank_deg);
    g_alpha = [cos(alpha), sin(alpha), 0; ...
              -sin(alpha), cos(alpha), 0; 0, 0, 1];
    g_bank = [1, 0, 0; 0, cos(sigma), sin(sigma); ...
              0, -sin(sigma), cos(sigma)];
    w_body = g_alpha * g_bank * w_t;
    boresight = boresight_body(:);
    cosine = dot(boresight, w_body) / ...
        (norm(boresight) * norm(w_body));
    raap_deg = rad2deg(acos(min(1, max(-1, cosine))));
end

function assert_aux_equal_local(actual, expected, label)
    assert(isequal(fieldnames(actual), fieldnames(expected)), ...
        '%s field contract changed.', label);
    names = fieldnames(expected);
    for ii = 1:numel(names)
        name = names{ii};
        assert_close_local(actual.(name), expected.(name), ...
            1e-13, 1e-13, sprintf('%s.%s', label, name));
    end
end

function assert_close_local(actual, expected, abs_tol, rel_tol, label)
    assert(isequal(size(actual), size(expected)), ...
        '%s size changed.', label);
    assert(isequal(isnan(actual), isnan(expected)), ...
        '%s NaN pattern changed.', label);
    assert(isequal(isinf(actual), isinf(expected)), ...
        '%s Inf pattern changed.', label);
    finite = isfinite(expected);
    if ~any(finite(:))
        return;
    end
    delta = abs(actual(finite) - expected(finite));
    tolerance = abs_tol + rel_tol * abs(expected(finite));
    assert(all(delta <= tolerance), ...
        '%s changed (maximum absolute error %.6g).', label, max(delta));
end

function assert_error_id_local(action, expected_id)
    thrown = false;
    try
        action();
    catch exception
        thrown = true;
        assert(strcmp(exception.identifier, expected_id), ...
            'Expected error %s but received %s.', ...
            expected_id, exception.identifier);
    end
    assert(thrown, 'Expected error %s was not thrown.', expected_id);
end

function value = field_or_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
end

function value = bool_or_local(s, name, default_value)
    value = field_or_local(s, name, default_value);
    if islogical(value)
        return;
    end
    if isnumeric(value)
        value = value ~= 0;
        return;
    end
    txt = lower(strtrim(string(value)));
    value = txt == "true" || txt == "1" || txt == "yes" || txt == "on";
end

function X = synthetic_entry_state_local(sys)
    radius = sys.Re + 120e3;
    speed = 7800;
    fpa = deg2rad(-6);
    X = zeros(14,1);
    X(1:3) = [radius; 0; 0];
    X(4:6) = [speed*sin(fpa); speed*cos(fpa); 0];
    X(7:10) = [0;0;0;1];
    X(14) = 1000;
end

function X = synthetic_relay_state_local(sys)
    radius = sys.Re + 500e3;
    X = zeros(14,1);
    X(1:3) = [radius; 0; 0];
    X(4:6) = [0; sqrt(sys.mu/radius); 0];
    X(7:10) = [0;0;0;1];
    X(14) = sys.Target_Mass;
end
