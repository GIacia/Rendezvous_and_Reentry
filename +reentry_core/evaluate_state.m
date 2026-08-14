function aux = evaluate_state(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    % Evaluate the shared seven-state atmospheric-entry physics model.
    % X = [r_ECI_m; v_ECI_m_s; mass_kg].

    if ~isnumeric(X) || ~isreal(X) || numel(X) ~= 7 || ...
            any(~isfinite(X(:))) || norm(X(1:3)) <= 0 || X(7) <= 0
        error('reentry_core:evaluateState:InvalidState', ...
              ['X must be a finite real seven-element state with nonzero ' ...
               'position and positive mass.']);
    end
    validateattributes(heat_k, {'numeric'}, ...
        {'real','scalar','finite','nonnegative'}, mfilename, 'heat_k');
    X = X(:);

    r = X(1:3);
    v = X(4:6);
    mass = X(7);
    altitude = norm(r) - sys.Re;

    atmosphere_opts = atmosphere_options(sys);
    [rho, ~, ~, a_sound] = Standard_Atmosphere_Density(altitude, atmosphere_opts);
    rho = rho * get_field(shape, 'density_scale', 1.0);

    if get_bool(atmosphere_opts, 'co_rotate_atmosphere', true)
        omega = get_field(atmosphere_opts, 'earth_rotation_rad_s', 7.2921159e-5);
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
    [aoa_actual_deg, bank_actual_deg] = reentry_core.resolve_commands( ...
        shape, aoa_deg, bank_angle_deg, speed_rel, mach);

    a_gravity = reentry_core.gravity_acceleration(r, sys, shape);
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
        [cd, cl, ld] = reentry_core.aerodynamic_coefficients( ...
            shape, aoa_actual_deg, mach);
        area = get_field(shape, 'reference_area_m2', 1.0);
        drag_accel_mag = dynamic_pressure * cd * area / mass;
        a_drag = -drag_accel_mag * v_rel_hat;

        if lift_enabled
            lift_dir = reentry_core.lift_direction(r, v_rel, bank_actual_deg);
            a_lift = ld * drag_accel_mag * lift_dir;
        else
            cl = 0;
            ld = 0;
        end

        nose_radius = max(get_field(shape, 'nose_radius_m', 0.05), 1e-4);
        heat_flux = heat_k * sqrt(rho / nose_radius) * speed_rel^3;
        a_aero = a_drag + a_lift;
    end

    aux = struct();
    aux.altitude = altitude;
    aux.rho = rho;
    aux.v_rel = v_rel;
    aux.speed_rel = speed_rel;
    aux.fpa_deg = reentry_core.flight_path_angle_deg(r, v);
    aux.fpa_rel_deg = reentry_core.flight_path_angle_deg(r, v_rel);
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
    guidance_threshold = get_field(shape, 'guidance_activation_drag_g', inf);
    aux.capsule_guidance_active = aux.drag_g >= guidance_threshold;
    guidance_cutoff_mach = get_field(shape, 'guidance_cutoff_mach', -inf);
    aux.capsule_guidance_window_active = ...
        aux.capsule_guidance_active && aux.mach >= guidance_cutoff_mach;
    aux.a_gravity = a_gravity;
    aux.a_aero = a_aero;
    aux.a_drag = a_drag;
    aux.a_lift = a_lift;
end
