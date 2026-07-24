function [X_rv_final, X_chaser_final, hist, summary] = Reentry_Propagator(sys, X_entry0, X_chaser_orbit0, phase3_elapsed_s, params)
    % Atmospheric re-entry propagation from the configured entry interface.
    %
    % The translational state remains ECI, but aerodynamic velocity is
    % computed relative to a co-rotating atmosphere. This gives the ECEF
    % transition needed for drag without changing the rest of the simulator's
    % inertial-state convention.

    if nargin < 5 || isempty(params)
        params = struct();
    end
    if nargin < 4 || isempty(phase3_elapsed_s)
        phase3_elapsed_s = 0;
    end

    shape = get_reentry_shape_local(sys, params);
    dt = get_param_local(params, 'dt', get_reentry_vehicle_field_local(sys, 'dt', 0.5));
    max_time = get_param_local(params, 'max_time', get_reentry_vehicle_field_local(sys, 'max_time', 2500));
    terminal_altitude = get_param_local(params, 'terminal_altitude', get_reentry_vehicle_field_local(sys, 'terminal_altitude', 20e3));
    aoa_deg = get_param_local(params, 'aoa_deg', get_reentry_vehicle_field_local(sys, 'aoa_deg', get_shape_field_local(shape, 'default_aoa_deg', 20)));
    bank_angle_deg = get_param_local(params, 'bank_angle_deg', get_reentry_vehicle_field_local(sys, 'bank_angle_deg', 0));
    lift_enabled = get_bool_param_local(params, 'lift_enabled', get_reentry_vehicle_field_local(sys, 'lift_enabled', true));
    los_margin_altitude = get_param_local(params, 'los_margin_altitude', get_reentry_vehicle_field_local(sys, 'los_margin_altitude', 0));
    chaser_dt = get_param_local(params, 'chaser_dt', max(1, min(20, dt * 20)));
    heat_k = get_param_local(params, 'sutton_graves_k', get_reentry_vehicle_field_local(sys, 'sutton_graves_k', 1.83e-4));

    if dt <= 0
        error('Reentry_Propagator requires positive dt.');
    end
    if max_time <= 0
        error('Reentry_Propagator requires positive max_time.');
    end

    X_chaser_final = propagate_orbit_chaser_local(X_chaser_orbit0, sys, phase3_elapsed_s, chaser_dt, 0);
    X_chaser = X_chaser_final;
    X_rv = [X_entry0(1:6); X_entry0(14)];

    hist = init_reentry_hist_local();
    elapsed = 0;
    heat_load = 0;
    [hist, aux_prev] = log_reentry_state_local(hist, X_rv, X_chaser, elapsed, heat_load, sys, shape, ...
                                               aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k);

    while elapsed < max_time - 1e-12
        altitude = norm(X_rv(1:3)) - sys.Re;
        if altitude <= terminal_altitude
            break;
        end

        dt_eff = min(dt, max_time - elapsed);
        X_rv = rk4_reentry_step_local(X_rv, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, dt_eff);
        X_chaser = rk4_chaser_step_local(X_chaser, sys, dt_eff, phase3_elapsed_s + elapsed);
        elapsed = elapsed + dt_eff;

        aux_now = reentry_aux_local(X_rv, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
        heat_load = heat_load + 0.5 * (aux_prev.heat_flux + aux_now.heat_flux) * dt_eff;
        [hist, aux_prev] = log_reentry_state_local(hist, X_rv, X_chaser, elapsed, heat_load, sys, shape, ...
                                                   aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k);
    end

    X_rv_final = zeros(14,1);
    X_rv_final(1:6) = X_rv(1:6);
    X_rv_final(7:10) = [0;0;0;1];
    X_rv_final(11:13) = [0;0;0];
    X_rv_final(14) = X_rv(7);
    X_chaser_final = X_chaser;

    summary = summarize_reentry_local(hist, shape, terminal_altitude);
end

function hist = init_reentry_hist_local()
    hist.time = [];
    hist.rv_pos = [];
    hist.rv_vel = [];
    hist.chaser_pos = [];
    hist.chaser_vel = [];
    hist.mass = [];
    hist.altitude = [];
    hist.speed = [];
    hist.speed_rel = [];
    hist.fpa_deg = [];
    hist.rho = [];
    hist.mach = [];
    hist.dynamic_pressure = [];
    hist.heat_flux = [];
    hist.heat_load = [];
    hist.g_load = [];
    hist.aoa_deg = [];
    hist.bank_angle_deg = [];
    hist.ld = [];
    hist.los_clear = [];
    hist.los_clearance = [];
    hist.los_elevation_deg = [];
end

function [hist, aux] = log_reentry_state_local(hist, X_rv, X_chaser, t, heat_load, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k)
    aux = reentry_aux_local(X_rv, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    [los_clear, los_clearance, los_elevation_deg] = line_of_sight_geometry_local(X_chaser(1:3), X_rv(1:3), sys, los_margin_altitude);

    hist.time = [hist.time, t];
    hist.rv_pos = [hist.rv_pos, X_rv(1:3)];
    hist.rv_vel = [hist.rv_vel, X_rv(4:6)];
    hist.chaser_pos = [hist.chaser_pos, X_chaser(1:3)];
    hist.chaser_vel = [hist.chaser_vel, X_chaser(4:6)];
    hist.mass = [hist.mass, X_rv(7)];
    hist.altitude = [hist.altitude, norm(X_rv(1:3)) - sys.Re];
    hist.speed = [hist.speed, norm(X_rv(4:6))];
    hist.speed_rel = [hist.speed_rel, aux.speed_rel];
    hist.fpa_deg = [hist.fpa_deg, aux.fpa_deg];
    hist.rho = [hist.rho, aux.rho];
    hist.mach = [hist.mach, aux.mach];
    hist.dynamic_pressure = [hist.dynamic_pressure, aux.dynamic_pressure];
    hist.heat_flux = [hist.heat_flux, aux.heat_flux];
    hist.heat_load = [hist.heat_load, heat_load];
    hist.g_load = [hist.g_load, aux.g_load];
    hist.aoa_deg = [hist.aoa_deg, aoa_deg];
    hist.bank_angle_deg = [hist.bank_angle_deg, bank_angle_deg];
    hist.ld = [hist.ld, aux.ld];
    hist.los_clear = [hist.los_clear, los_clear];
    hist.los_clearance = [hist.los_clearance, los_clearance];
    hist.los_elevation_deg = [hist.los_elevation_deg, los_elevation_deg];
end

function summary = summarize_reentry_local(hist, shape, terminal_altitude)
    summary = struct();
    summary.shape_name = string(shape.name);
    summary.duration_s = hist.time(end);
    summary.terminal_altitude_m = terminal_altitude;
    summary.final_altitude_m = hist.altitude(end);
    summary.initial_fpa_deg = hist.fpa_deg(1);
    summary.final_fpa_deg = hist.fpa_deg(end);
    summary.aoa_deg = hist.aoa_deg(1);
    summary.bank_angle_deg = hist.bank_angle_deg(1);
    summary.max_heat_flux_W_m2 = max(hist.heat_flux);
    summary.total_heat_load_J_m2 = hist.heat_load(end);
    summary.max_dynamic_pressure_Pa = max(hist.dynamic_pressure);
    summary.max_g_load = max(hist.g_load);
    summary.los_maintained = all(hist.los_clear);
    summary.min_los_clearance_m = min(hist.los_clearance);
    summary.min_los_elevation_deg = min(hist.los_elevation_deg);

    loss_idx = find(~hist.los_clear, 1, 'first');
    if isempty(loss_idx)
        summary.first_los_loss_time_s = NaN;
    else
        summary.first_los_loss_time_s = hist.time(loss_idx);
    end
end

function X_next = rk4_reentry_step_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, dt)
    f = @(X_now) reentry_dynamics_local(X_now, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    k1 = f(X);
    k2 = f(X + 0.5 * dt * k1);
    k3 = f(X + 0.5 * dt * k2);
    k4 = f(X + dt * k3);
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function dX = reentry_dynamics_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    aux = reentry_aux_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    dX = [X(4:6); aux.a_gravity + aux.a_aero; 0];
end

function aux = reentry_aux_local(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    r = X(1:3);
    v = X(4:6);
    mass = X(7);
    altitude = norm(r) - sys.Re;

    atmosphere_opts = get_atmosphere_opts_local(sys);
    [rho, ~, ~, a_sound] = Standard_Atmosphere_Density(altitude, atmosphere_opts);

    omega = get_atmos_field_local(atmosphere_opts, 'earth_rotation_rad_s', 7.2921159e-5);
    v_atm = cross([0;0;omega], r);
    v_rel = v - v_atm;
    speed_rel = norm(v_rel);

    a_gravity = gravity_j2_local(r, sys);
    a_aero = zeros(3,1);
    a_drag = zeros(3,1);
    a_lift = zeros(3,1);
    ld = 0;
    dynamic_pressure = 0;
    heat_flux = 0;

    if speed_rel > 0 && rho > 0 && mass > 0
        v_rel_hat = v_rel / speed_rel;
        dynamic_pressure = 0.5 * rho * speed_rel^2;
        cd = get_shape_field_local(shape, 'cd', 1.2);
        area = get_shape_field_local(shape, 'reference_area_m2', 1.0);
        drag_accel_mag = dynamic_pressure * cd * area / mass;
        a_drag = -drag_accel_mag * v_rel_hat;

        if lift_enabled
            mach = speed_rel / max(a_sound, eps);
            ld = lookup_ld_local(shape, aoa_deg, mach);
            lift_dir = lift_direction_local(r, v_rel, bank_angle_deg);
            a_lift = ld * drag_accel_mag * lift_dir;
        end

        nose_radius = max(get_shape_field_local(shape, 'nose_radius_m', 0.05), 1e-4);
        heat_flux = heat_k * sqrt(rho / nose_radius) * speed_rel^3;
        a_aero = a_drag + a_lift;
    end

    aux = struct();
    aux.rho = rho;
    aux.speed_rel = speed_rel;
    aux.fpa_deg = flight_path_angle_deg_local(r, v);
    aux.mach = speed_rel / max(a_sound, eps);
    aux.dynamic_pressure = dynamic_pressure;
    aux.heat_flux = heat_flux;
    aux.g_load = norm(a_aero) / sys.g0;
    aux.ld = ld;
    aux.a_gravity = a_gravity;
    aux.a_aero = a_aero;
    aux.a_drag = a_drag;
    aux.a_lift = a_lift;
end

function fpa_deg = flight_path_angle_deg_local(r, v)
    r_hat = r / norm(r);
    v_radial = dot(v, r_hat);
    v_horizontal = norm(v - v_radial * r_hat);
    fpa_deg = rad2deg(atan2(v_radial, v_horizontal));
end

function a = gravity_j2_local(r, sys)
    r_norm = norm(r);
    a_g = -sys.mu / r_norm^3 * r;
    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1);
                      (r(2)/r_norm)*(5*z2 - 1);
                      (r(3)/r_norm)*(5*z2 - 3) ];
    a = a_g + a_j2;
end

function dir = lift_direction_local(r, v_rel, bank_angle_deg)
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
    dir = dir*cos(bank) + cross(vhat, dir)*sin(bank) + vhat*dot(vhat, dir)*(1-cos(bank));
    dir = dir / norm(dir);
end

function [clear, clearance, elevation_deg] = line_of_sight_geometry_local(r_chaser, r_rv, sys, margin_altitude)
    los = r_chaser - r_rv;
    los_norm = norm(los);
    if los_norm <= eps
        clearance = norm(r_rv) - sys.Re;
        clear = clearance > margin_altitude;
        elevation_deg = 90;
        return;
    end

    u = r_chaser - r_rv;
    tau = -dot(r_rv, u) / dot(u, u);
    tau = min(1, max(0, tau));
    closest = r_rv + tau * u;
    clearance = norm(closest) - sys.Re;
    clear = clearance > margin_altitude;

    up = r_rv / norm(r_rv);
    elevation_deg = rad2deg(asin(dot(los / los_norm, up)));
end

function X = propagate_orbit_chaser_local(X, sys, duration, dt, t0)
    if isempty(X) || duration <= 0
        return;
    end

    elapsed = 0;
    while elapsed < duration - 1e-12
        dt_eff = min(dt, duration - elapsed);
        X = rk4_chaser_step_local(X, sys, dt_eff, t0 + elapsed);
        elapsed = elapsed + dt_eff;
    end
end

function X_next = rk4_chaser_step_local(X, sys, dt, t_abs)
    k1 = Env_EOM(t_abs,          X,             [0;0;0], [0;0;0], sys, true, "chaser");
    k2 = Env_EOM(t_abs + dt/2,   X + k1*dt/2,   [0;0;0], [0;0;0], sys, true, "chaser");
    k3 = Env_EOM(t_abs + dt/2,   X + k2*dt/2,   [0;0;0], [0;0;0], sys, true, "chaser");
    k4 = Env_EOM(t_abs + dt,     X + k3*dt,     [0;0;0], [0;0;0], sys, true, "chaser");
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    if norm(X_next(7:10)) > 0
        X_next(7:10) = X_next(7:10) / norm(X_next(7:10));
    end
end

function shape = get_reentry_shape_local(sys, params)
    default_shape = get_reentry_vehicle_field_local(sys, 'selected_shape', "COMPROMISE");
    name = get_string_param_local(params, 'shape_name', default_shape);
    key = normalize_shape_key_local(name);

    if ~isfield(sys, 'reentry_vehicle') || ~isfield(sys.reentry_vehicle, 'shapes') || ...
            ~isfield(sys.reentry_vehicle.shapes, key)
        error('Unknown reentry vehicle shape: %s.', char(name));
    end

    shape = sys.reentry_vehicle.shapes.(key);
    shape.name = string(get_shape_field_local(shape, 'name', key));
end

function key = normalize_shape_key_local(name)
    txt = upper(strtrim(string(name)));
    txt = replace(txt, "-", "_");
    txt = replace(txt, " ", "_");

    if txt == "HEATLOADMIN" || txt == "HEAT_LOAD_MIN"
        txt = "HEATLOAD_MIN";
    elseif txt == "PAYLOADMAX" || txt == "PAY_LOAD_MAX"
        txt = "PAYLOAD_MAX";
    elseif txt == "TPSMIN"
        txt = "TPS_MIN";
    end

    key = char(txt);
end

function ld = lookup_ld_local(shape, aoa_deg, ~)
    aoa_grid = get_shape_field_local(shape, 'ld_aoa_deg', [0 20 80]);
    ld_grid = get_shape_field_local(shape, 'ld_values', [0 1 0]);
    aoa = min(max(aoa_deg, min(aoa_grid)), max(aoa_grid));
    ld = interp1(aoa_grid, ld_grid, aoa, 'pchip');
    ld = max(0, ld);
end

function atmosphere_opts = get_atmosphere_opts_local(sys)
    atmosphere_opts = struct();
    if isfield(sys, 'environment') && isfield(sys.environment, 'atmospheric_drag')
        atmosphere_opts = sys.environment.atmospheric_drag;
    end
    atmosphere_opts.enabled = true;
    atmosphere_opts.model = "ISA76";
end

function value = get_reentry_vehicle_field_local(sys, name, default_value)
    if isfield(sys, 'reentry_vehicle') && isfield(sys.reentry_vehicle, name) && ~isempty(sys.reentry_vehicle.(name))
        value = sys.reentry_vehicle.(name);
    else
        value = default_value;
    end
end

function value = get_shape_field_local(shape, name, default_value)
    if isstruct(shape) && isfield(shape, name) && ~isempty(shape.(name))
        value = shape.(name);
    else
        value = default_value;
    end
end

function value = get_atmos_field_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
end

function value = get_param_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
end

function value = get_string_param_local(s, name, default_value)
    value = get_param_local(s, name, default_value);
    value = string(value);
end

function value = get_bool_param_local(s, name, default_value)
    value = get_param_local(s, name, default_value);
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
