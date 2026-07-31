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
    vehicle_mode = string(get_shape_field_local(shape, 'vehicle_mode', "SPACEPLANE"));
    dt = get_param_local(params, 'dt', get_reentry_vehicle_field_local(sys, 'dt', 0.5));
    max_time = get_param_local(params, 'max_time', get_reentry_vehicle_field_local(sys, 'max_time', 2500));
    terminal_altitude = get_param_local(params, 'terminal_altitude', get_reentry_vehicle_field_local(sys, 'terminal_altitude', 20e3));
    aoa_deg = get_param_local(params, 'aoa_deg', get_reentry_vehicle_field_local(sys, 'aoa_deg', get_shape_field_local(shape, 'default_aoa_deg', 20)));
    bank_angle_deg = get_param_local(params, 'bank_angle_deg', get_reentry_vehicle_field_local(sys, 'bank_angle_deg', 0));
    lift_enabled = get_bool_param_local(params, 'lift_enabled', get_reentry_vehicle_field_local(sys, 'lift_enabled', true));
    los_margin_altitude = get_param_local(params, 'los_margin_altitude', get_reentry_vehicle_field_local(sys, 'los_margin_altitude', 0));
    chaser_dt = get_param_local(params, 'chaser_dt', max(1, min(20, dt * 20)));
    heat_k = get_param_local(params, 'sutton_graves_k', get_reentry_vehicle_field_local(sys, 'sutton_graves_k', 1.83e-4));
    terminal_speed = get_shape_field_local(shape, 'terminal_speed_m_s', NaN);
    altitude_termination_enabled = get_bool_param_local(shape, 'altitude_termination_enabled', true);
    safety_floor_altitude = get_shape_field_local(shape, 'safety_floor_altitude_m', 0);
    relay_mode = upper(string(get_shape_field_local(shape, 'relay_mode', "MISSION_TARGET_DYNAMIC_ORBIT")));

    if dt <= 0
        error('Reentry_Propagator requires positive dt.');
    end
    if max_time <= 0
        error('Reentry_Propagator requires positive max_time.');
    end
    if ~isscalar(safety_floor_altitude) || ~isfinite(safety_floor_altitude)
        error('Reentry safety_floor_altitude_m must be a finite scalar.');
    end

    if numel(X_chaser_orbit0) == 6
        relay_mass = get_struct_field_local(sys, 'Target_Mass', 1);
        X_chaser_orbit0 = [X_chaser_orbit0(1:6); 0;0;0;1; 0;0;0; relay_mass];
    end
    if relay_mode == "PAPER_TDRS_STATIC_EARTH_FIXED"
        X_chaser_final = paper_tdrs_state_local(sys, shape, phase3_elapsed_s);
    else
        X_chaser_final = propagate_orbit_chaser_local(X_chaser_orbit0, sys, phase3_elapsed_s, chaser_dt, 0);
    end
    X_chaser = X_chaser_final;
    X_rv = [X_entry0(1:6); X_entry0(14)];
    if vehicle_mode == "CAPSULE"
        separation_mode = upper(string(get_param_local(params, 'separation_mode', ...
            get_shape_field_local(shape, 'separation_mode', "ENTRY_INTERFACE"))));
        if separation_mode == "ENTRY_INTERFACE"
            capsule_mass = get_param_local(params, 'capsule_mass_kg', ...
                get_shape_field_local(shape, 'mass_kg', X_rv(7)));
            if ~isscalar(capsule_mass) || ~isfinite(capsule_mass) || capsule_mass <= 0
                error('Capsule separation mass must be a finite positive scalar.');
            end
            if capsule_mass > X_rv(7)
                error('Capsule separation mass %.3f kg exceeds entry stack mass %.3f kg.', ...
                      capsule_mass, X_rv(7));
            end
            X_rv(7) = capsule_mass;
        elseif separation_mode ~= "ATTACHED"
            error('Unknown capsule separation_mode: %s. Use ENTRY_INTERFACE or ATTACHED.', ...
                  char(separation_mode));
        end
    end

    hist = init_reentry_hist_local();
    elapsed = 0;
    heat_load = 0;
    integration_steps = 0;
    termination_reason = "MAX_TIME";
    [hist, aux_prev] = log_reentry_state_local(hist, X_rv, X_chaser, elapsed, heat_load, sys, shape, ...
                                               aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k);

    while elapsed < max_time - 1e-12
        altitude = norm(X_rv(1:3)) - sys.Re;
        if altitude_termination_enabled && altitude <= terminal_altitude
            termination_reason = "TERMINAL_ALTITUDE";
            break;
        elseif ~altitude_termination_enabled && altitude <= safety_floor_altitude
            termination_reason = "GROUND_IMPACT";
            break;
        end
        if isfinite(terminal_speed) && aux_prev.speed_rel <= terminal_speed
            termination_reason = "PARACHUTE_SPEED";
            break;
        end

        dt_eff = min(dt, max_time - elapsed);
        X_rv_prev = X_rv;
        X_rv_trial = rk4_reentry_step_local(X_rv_prev, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, dt_eff);
        altitude_trial = norm(X_rv_trial(1:3)) - sys.Re;
        aux_trial = reentry_aux_local(X_rv_trial, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);

        if altitude_termination_enabled
            altitude_event_target = terminal_altitude;
        else
            altitude_event_target = safety_floor_altitude;
        end
        altitude_crossed = altitude_trial <= altitude_event_target;
        speed_crossed = isfinite(terminal_speed) && aux_trial.speed_rel <= terminal_speed;

        if altitude_crossed
            [X_altitude_cross, t_altitude_cross] = refine_reentry_altitude_crossing_local( ...
                X_rv_prev, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, ...
                altitude_event_target, dt_eff);
        else
            X_altitude_cross = [];
            t_altitude_cross = inf;
        end

        if speed_crossed
            [X_speed_cross, t_speed_cross] = refine_reentry_speed_crossing_local( ...
                X_rv_prev, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, ...
                terminal_speed, dt_eff);
        else
            X_speed_cross = [];
            t_speed_cross = inf;
        end

        if isfinite(t_altitude_cross) && t_altitude_cross <= t_speed_cross
            X_rv = X_altitude_cross;
            dt_advance = t_altitude_cross;
            if altitude_termination_enabled
                termination_reason = "TERMINAL_ALTITUDE";
            else
                termination_reason = "GROUND_IMPACT";
            end
        elseif isfinite(t_speed_cross) && t_speed_cross < t_altitude_cross
            X_rv = X_speed_cross;
            dt_advance = t_speed_cross;
            termination_reason = "PARACHUTE_SPEED";
        else
            X_rv = X_rv_trial;
            dt_advance = dt_eff;
        end

        if relay_mode == "PAPER_TDRS_STATIC_EARTH_FIXED"
            X_chaser = paper_tdrs_state_local( ...
                sys, shape, phase3_elapsed_s + elapsed + dt_advance);
        else
            X_chaser = rk4_chaser_step_local( ...
                X_chaser, sys, dt_advance, phase3_elapsed_s + elapsed);
        end
        elapsed = elapsed + dt_advance;
        integration_steps = integration_steps + 1;

        aux_now = reentry_aux_local(X_rv, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
        heat_load = heat_load + 0.5 * (aux_prev.heat_flux + aux_now.heat_flux) * dt_advance;
        [hist, aux_prev] = log_reentry_state_local(hist, X_rv, X_chaser, elapsed, heat_load, sys, shape, ...
                                                   aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k);

        if termination_reason == "TERMINAL_ALTITUDE" || ...
                termination_reason == "PARACHUTE_SPEED" || ...
                termination_reason == "GROUND_IMPACT"
            break;
        end
    end

    X_rv_final = zeros(14,1);
    X_rv_final(1:6) = X_rv(1:6);
    X_rv_final(7:10) = [0;0;0;1];
    X_rv_final(11:13) = [0;0;0];
    X_rv_final(14) = X_rv(7);
    X_chaser_final = X_chaser;

    summary = summarize_reentry_local(hist, shape, terminal_altitude, termination_reason, integration_steps);
end

function hist = init_reentry_hist_local()
    hist.time = [];
    hist.rv_pos = [];
    hist.rv_vel = [];
    hist.relay_pos = [];
    hist.relay_vel = [];
    % Deprecated aliases retained for existing plotting code.
    hist.chaser_pos = [];
    hist.chaser_vel = [];
    hist.mass = [];
    hist.altitude = [];
    hist.speed = [];
    hist.speed_rel = [];
    hist.fpa_deg = [];
    hist.fpa_rel_deg = [];
    hist.rho = [];
    hist.mach = [];
    hist.dynamic_pressure = [];
    hist.heat_flux = [];
    hist.heat_load = [];
    hist.g_load = [];
    hist.aoa_deg = [];
    hist.bank_angle_deg = [];
    hist.cd = [];
    hist.cl = [];
    hist.ld = [];
    hist.drag_g = [];
    hist.path_constraint_satisfied = [];
    hist.los_clear = [];
    hist.los_clearance = [];
    hist.los_elevation_deg = [];
    hist.raap_deg = [];
    hist.relay_range_m = [];
    hist.antenna_constraint_active = [];
    hist.blackout_zone_active = [];
    hist.blackout_zone_phase = [];
    hist.antenna_in_beam = [];
    hist.antenna_within_cone = [];
    hist.antenna_within_range = [];
    hist.communication_geometry_available = [];
    hist.best_feasible_raap_deg = [];
    hist.best_feasible_bank_deg = [];
    hist.capsule_guidance_active = [];
    hist.capsule_guidance_window_active = [];
    hist.capsule_guidance_threshold_met = [];
end

function [hist, aux] = log_reentry_state_local(hist, X_rv, X_chaser, t, heat_load, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, los_margin_altitude, heat_k)
    aux = reentry_aux_local(X_rv, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    if isempty(X_chaser) || numel(X_chaser) < 6
        los_clear = NaN;
        los_clearance = NaN;
        los_elevation_deg = NaN;
        chaser_pos = nan(3,1);
        chaser_vel = nan(3,1);
    else
        [los_clear, los_clearance, los_elevation_deg] = line_of_sight_geometry_local(X_chaser(1:3), X_rv(1:3), sys, los_margin_altitude);
        chaser_pos = X_chaser(1:3);
        chaser_vel = X_chaser(4:6);
    end
    [blackout_active, antenna_constraint_active, blackout_phase] = ...
        antenna_constraint_activity_local(hist, aux, shape);
    antenna = antenna_geometry_local(X_rv, X_chaser, aux, shape, los_clear);
    antenna.blackout_active = blackout_active;
    antenna.constraint_active = antenna_constraint_active;
    path_constraint_satisfied = path_constraints_satisfied_local(aux, shape);

    hist.time = [hist.time, t];
    hist.rv_pos = [hist.rv_pos, X_rv(1:3)];
    hist.rv_vel = [hist.rv_vel, X_rv(4:6)];
    hist.relay_pos = [hist.relay_pos, chaser_pos];
    hist.relay_vel = [hist.relay_vel, chaser_vel];
    hist.chaser_pos = [hist.chaser_pos, chaser_pos];
    hist.chaser_vel = [hist.chaser_vel, chaser_vel];
    hist.mass = [hist.mass, X_rv(7)];
    hist.altitude = [hist.altitude, norm(X_rv(1:3)) - sys.Re];
    hist.speed = [hist.speed, norm(X_rv(4:6))];
    hist.speed_rel = [hist.speed_rel, aux.speed_rel];
    hist.fpa_deg = [hist.fpa_deg, aux.fpa_deg];
    hist.fpa_rel_deg = [hist.fpa_rel_deg, aux.fpa_rel_deg];
    hist.rho = [hist.rho, aux.rho];
    hist.mach = [hist.mach, aux.mach];
    hist.dynamic_pressure = [hist.dynamic_pressure, aux.dynamic_pressure];
    hist.heat_flux = [hist.heat_flux, aux.heat_flux];
    hist.heat_load = [hist.heat_load, heat_load];
    hist.g_load = [hist.g_load, aux.g_load];
    hist.aoa_deg = [hist.aoa_deg, aux.aoa_deg];
    hist.bank_angle_deg = [hist.bank_angle_deg, aux.bank_angle_deg];
    hist.cd = [hist.cd, aux.cd];
    hist.cl = [hist.cl, aux.cl];
    hist.ld = [hist.ld, aux.ld];
    hist.drag_g = [hist.drag_g, aux.drag_g];
    hist.path_constraint_satisfied = [hist.path_constraint_satisfied, path_constraint_satisfied];
    hist.los_clear = [hist.los_clear, los_clear];
    hist.los_clearance = [hist.los_clearance, los_clearance];
    hist.los_elevation_deg = [hist.los_elevation_deg, los_elevation_deg];
    hist.raap_deg = [hist.raap_deg, antenna.raap_deg];
    hist.relay_range_m = [hist.relay_range_m, antenna.range_m];
    hist.antenna_constraint_active = [hist.antenna_constraint_active, antenna.constraint_active];
    hist.blackout_zone_active = [hist.blackout_zone_active, antenna.blackout_active];
    hist.blackout_zone_phase = [hist.blackout_zone_phase, blackout_phase];
    hist.antenna_in_beam = [hist.antenna_in_beam, antenna.in_beam];
    hist.antenna_within_cone = [hist.antenna_within_cone, antenna.within_cone];
    hist.antenna_within_range = [hist.antenna_within_range, antenna.within_range];
    hist.communication_geometry_available = [hist.communication_geometry_available, antenna.communication_available];
    hist.best_feasible_raap_deg = [hist.best_feasible_raap_deg, antenna.best_raap_deg];
    hist.best_feasible_bank_deg = [hist.best_feasible_bank_deg, antenna.best_bank_deg];
    threshold_met = aux.capsule_guidance_active;
    if isempty(hist.capsule_guidance_active)
        guidance_latched = threshold_met;
    else
        guidance_latched = hist.capsule_guidance_active(end) || threshold_met;
    end
    guidance_cutoff_mach = get_shape_field_local(shape, 'guidance_cutoff_mach', -inf);
    hist.capsule_guidance_threshold_met = [hist.capsule_guidance_threshold_met, threshold_met];
    hist.capsule_guidance_active = [hist.capsule_guidance_active, guidance_latched];
    hist.capsule_guidance_window_active = [hist.capsule_guidance_window_active, ...
        guidance_latched && aux.mach >= guidance_cutoff_mach];
end

function summary = summarize_reentry_local(hist, shape, terminal_altitude, termination_reason, integration_steps)
    summary = struct();
    summary.vehicle_mode = string(get_shape_field_local(shape, 'vehicle_mode', "SPACEPLANE"));
    summary.shape_name = string(shape.name);
    summary.separation_mode = string(get_shape_field_local(shape, 'separation_mode', "NOT_APPLICABLE"));
    summary.duration_s = hist.time(end);
    summary.termination_reason = string(termination_reason);
    summary.reached_terminal_altitude = termination_reason == "TERMINAL_ALTITUDE";
    summary.reached_parachute_speed = termination_reason == "PARACHUTE_SPEED";
    summary.reached_ground_impact = termination_reason == "GROUND_IMPACT";
    summary.completed_nominally = summary.reached_terminal_altitude || summary.reached_parachute_speed;
    summary.altitude_termination_enabled = get_bool_param_local(shape, 'altitude_termination_enabled', true);
    summary.safety_floor_altitude_m = get_shape_field_local(shape, 'safety_floor_altitude_m', 0);
    summary.integration_steps = integration_steps;
    summary.terminal_altitude_m = terminal_altitude;
    summary.final_altitude_m = hist.altitude(end);
    summary.minimum_altitude_m = min(hist.altitude);
    if summary.altitude_termination_enabled
        summary.terminal_altitude_error_m = hist.altitude(end) - terminal_altitude;
    else
        summary.terminal_altitude_error_m = NaN;
    end
    summary.ascending_at_termination = dot(hist.rv_vel(:,end), hist.rv_pos(:,end)) > 0;
    summary.initial_fpa_deg = hist.fpa_deg(1);
    summary.final_fpa_deg = hist.fpa_deg(end);
    summary.initial_fpa_rel_deg = hist.fpa_rel_deg(1);
    summary.final_fpa_rel_deg = hist.fpa_rel_deg(end);
    summary.aoa_deg = hist.aoa_deg(1);
    summary.initial_aoa_deg = hist.aoa_deg(1);
    summary.final_aoa_deg = hist.aoa_deg(end);
    summary.bank_angle_deg = hist.bank_angle_deg(1);
    summary.max_heat_flux_W_m2 = max(hist.heat_flux);
    summary.total_heat_load_J_m2 = hist.heat_load(end);
    summary.max_dynamic_pressure_Pa = max(hist.dynamic_pressure);
    summary.max_g_load = max(hist.g_load);
    summary.max_drag_g = max(hist.drag_g);
    summary.max_bank_rate_deg_s = max_rate_local(hist.bank_angle_deg, hist.time);
    constraints = get_shape_field_local(shape, 'constraints', struct());
    summary.dynamic_pressure_limit_Pa = get_struct_field_local(constraints, 'max_dynamic_pressure_Pa', inf);
    summary.g_load_limit = get_struct_field_local(constraints, 'max_g_load', inf);
    summary.heat_flux_limit_W_m2 = get_struct_field_local(constraints, 'max_heat_flux_W_m2', inf);
    summary.bank_angle_limit_deg = get_struct_field_local(constraints, 'max_bank_angle_deg', inf);
    summary.bank_rate_limit_deg_s = get_struct_field_local(constraints, 'max_bank_rate_deg_s', inf);
    summary.evaluated_path_constraints_satisfied = all(hist.path_constraint_satisfied ~= 0) && ...
        summary.max_bank_rate_deg_s <= summary.bank_rate_limit_deg_s;
    summary.heat_constraint_comparable = get_bool_param_local(shape, 'heat_constraint_comparable', true);
    summary.heat_constraint_evaluated = isfinite(summary.heat_flux_limit_W_m2) && ...
        summary.heat_constraint_comparable;
    summary.surrogate_heat_below_reference_limit = ...
        ~isfinite(summary.heat_flux_limit_W_m2) || ...
        summary.max_heat_flux_W_m2 <= summary.heat_flux_limit_W_m2;
    if summary.heat_constraint_evaluated
        summary.heat_constraint_satisfied = summary.surrogate_heat_below_reference_limit;
    else
        summary.heat_constraint_satisfied = NaN;
    end
    summary.path_constraints_fully_evaluated = ...
        ~isfinite(summary.heat_flux_limit_W_m2) || summary.heat_constraint_evaluated;
    if ~summary.evaluated_path_constraints_satisfied
        summary.path_constraints_satisfied = false;
    elseif summary.path_constraints_fully_evaluated
        summary.path_constraints_satisfied = true;
    else
        summary.path_constraints_satisfied = NaN;
    end
    valid_los = isfinite(hist.los_clear);
    summary.los_evaluated = any(valid_los);
    if summary.los_evaluated
        summary.los_maintained = all(hist.los_clear(valid_los) ~= 0);
        summary.min_los_clearance_m = min(hist.los_clearance(valid_los));
        summary.min_los_elevation_deg = min(hist.los_elevation_deg(valid_los));
    else
        summary.los_maintained = NaN;
        summary.min_los_clearance_m = NaN;
        summary.min_los_elevation_deg = NaN;
    end

    loss_idx = find(valid_los & hist.los_clear == 0, 1, 'first');
    if isempty(loss_idx)
        summary.first_los_loss_time_s = NaN;
    else
        summary.first_los_loss_time_s = hist.time(loss_idx);
    end

    cd = get_shape_field_local(shape, 'cd', 1.2);
    area = get_shape_field_local(shape, 'reference_area_m2', 1.0);
    summary.initial_mass_kg = hist.mass(1);
    if isfinite(hist.cd(1)) && hist.cd(1) > 0
        cd = hist.cd(1);
    end
    summary.ballistic_coefficient_kg_m2 = hist.mass(1) / (cd * area);
    summary.initial_cd = hist.cd(1);
    summary.initial_cl = hist.cl(1);
    summary.aerodynamic_model = string(get_shape_field_local(shape, 'aero_model', "CONSTANT_CD_WITH_AOA_ONLY_LD"));
    summary.attitude_guidance_model = string(get_shape_field_local(shape, 'aoa_profile_mode', "CONSTANT_AOA_AND_BANK"));
    summary.gravity_model = string(get_shape_field_local(shape, 'gravity_model', "CENTRAL_SPHERICAL"));
    summary.heating_model = string(get_shape_field_local(shape, 'heating_model', "SUTTON_GRAVES_SURROGATE"));
    summary.paper_heating_model = string(get_shape_field_local(shape, 'paper_heating_model', "NOT_SPECIFIED"));
    summary.relay_mode = string(get_shape_field_local(shape, 'relay_mode', "MISSION_TARGET_DYNAMIC_ORBIT"));
    summary.paper_initial_altitude_m = get_shape_field_local(shape, 'paper_initial_altitude_m', NaN);
    summary.paper_terminal_altitude_m = get_shape_field_local(shape, 'paper_terminal_altitude_m', NaN);
    summary.paper_terminal_speed_m_s = get_shape_field_local(shape, 'paper_terminal_speed_m_s', NaN);

    valid_raap = isfinite(hist.raap_deg);
    valid_antenna_geometry = valid_raap & isfinite(hist.relay_range_m) & ...
        isfinite(hist.los_clear);
    requested_antenna = hist.antenna_constraint_active ~= 0;
    requested_blackout = hist.blackout_zone_active ~= 0;
    active_antenna = requested_antenna & valid_antenna_geometry;
    summary.antenna_evaluated = any(valid_raap);
    summary.antenna_beam_half_angle_deg = get_shape_field_local(shape, 'antenna_beam_half_angle_deg', NaN);
    summary.antenna_min_range_m = get_shape_field_local(shape, 'antenna_min_range_m', NaN);
    summary.antenna_max_range_m = get_shape_field_local(shape, 'antenna_max_range_m', NaN);
    if summary.antenna_evaluated
        summary.min_raap_deg = min(hist.raap_deg(valid_raap));
        summary.max_raap_deg = max(hist.raap_deg(valid_raap));
        summary.min_relay_range_m = min(hist.relay_range_m(valid_raap));
        summary.max_relay_range_m = max(hist.relay_range_m(valid_raap));
    else
        summary.min_raap_deg = NaN;
        summary.max_raap_deg = NaN;
        summary.min_relay_range_m = NaN;
        summary.max_relay_range_m = NaN;
    end
    antenna_interval_present = any(requested_antenna);
    antenna_interval_fully_evaluated = antenna_interval_present && ...
        all(valid_antenna_geometry(requested_antenna));
    antenna_any_evaluated_violation = any(active_antenna & ...
        hist.communication_geometry_available == 0);
    summary.antenna_constraint_any_evaluated_violation = ...
        antenna_any_evaluated_violation;
    if antenna_any_evaluated_violation
        summary.antenna_tracking_maintained = false;
    elseif antenna_interval_fully_evaluated
        summary.antenna_tracking_maintained = true;
    else
        summary.antenna_tracking_maintained = NaN;
    end
    bank_scan_enabled = get_bool_param_local(shape, ...
        'evaluate_bank_feasibility', false);
    summary.antenna_bank_scan_enabled = bank_scan_enabled;
    angle_tolerance = get_shape_field_local(shape, ...
        'antenna_angle_tolerance_deg', 1e-10);
    bank_sample_feasible = bank_scan_enabled & ...
        hist.best_feasible_raap_deg <= ...
            summary.antenna_beam_half_angle_deg + angle_tolerance & ...
        hist.antenna_within_range ~= 0;
    antenna_bank_any_evaluated_violation = bank_scan_enabled && ...
        any(active_antenna & ~bank_sample_feasible);
    summary.antenna_bank_any_evaluated_violation = ...
        antenna_bank_any_evaluated_violation;
    summary.antenna_bank_scan_fully_evaluated = bank_scan_enabled && ...
        antenna_interval_fully_evaluated;
    if summary.antenna_bank_scan_fully_evaluated
        summary.max_best_feasible_raap_deg = max(hist.best_feasible_raap_deg(active_antenna));
    else
        summary.max_best_feasible_raap_deg = NaN;
    end
    if antenna_bank_any_evaluated_violation
        summary.antenna_bank_feasible = false;
    elseif summary.antenna_bank_scan_fully_evaluated
        summary.antenna_bank_feasible = true;
    else
        summary.antenna_bank_feasible = NaN;
    end
    summary.max_instantaneous_best_raap_deg = summary.max_best_feasible_raap_deg;
    summary.instantaneous_bank_geometric_reachable = summary.antenna_bank_feasible;
    blackout_interval_present = any(requested_blackout);
    blackout_interval_fully_evaluated = blackout_interval_present && ...
        all(valid_antenna_geometry(requested_blackout));
    active_blackout = requested_blackout & valid_antenna_geometry;
    bzc_any_evaluated_violation = any(active_blackout & ...
        hist.communication_geometry_available == 0);
    summary.bzc_constraint_evaluated = blackout_interval_fully_evaluated;
    summary.bzc_constraint_any_evaluated_violation = ...
        bzc_any_evaluated_violation;
    if bzc_any_evaluated_violation
        summary.bzc_constraint_maintained = false;
    elseif blackout_interval_fully_evaluated
        summary.bzc_constraint_maintained = true;
    else
        summary.bzc_constraint_maintained = NaN;
    end

    summary.antenna_constraint_required = get_bool_param_local(shape, 'antenna_enabled', false);
    summary.antenna_constraint_fully_evaluated = ...
        ~summary.antenna_constraint_required || antenna_interval_fully_evaluated;
    antenna_constraint_violated = summary.antenna_constraint_required && ...
        antenna_any_evaluated_violation;
    summary.any_evaluated_constraint_violated = ...
        ~summary.evaluated_path_constraints_satisfied || antenna_constraint_violated;
    summary.all_configured_constraints_evaluated = ...
        summary.path_constraints_fully_evaluated && ...
        summary.antenna_constraint_fully_evaluated;
    if summary.any_evaluated_constraint_violated
        summary.overall_constraints_satisfied = false;
    elseif summary.all_configured_constraints_evaluated
        summary.overall_constraints_satisfied = true;
    else
        summary.overall_constraints_satisfied = NaN;
    end

    summary.capsule_guidance_activation_reached = any(hist.capsule_guidance_active ~= 0);
    summary.capsule_guidance_window_entered = any(hist.capsule_guidance_window_active ~= 0);
    summary.capsule_guidance_threshold_first_time_s = NaN;
    threshold_idx = find(hist.capsule_guidance_threshold_met ~= 0, 1, 'first');
    if ~isempty(threshold_idx)
        summary.capsule_guidance_threshold_first_time_s = hist.time(threshold_idx);
    end
    summary.reference_entry_interface_altitude_m = get_shape_field_local(shape, 'reference_entry_interface_altitude_m', NaN);
    summary.reference_entry_speed_m_s = get_shape_field_local(shape, 'reference_entry_speed_m_s', NaN);
    summary.reference_entry_fpa_deg = get_shape_field_local(shape, 'reference_entry_fpa_deg', NaN);
    summary.paper_entry_condition_scope = string(get_shape_field_local(shape, 'paper_entry_condition_scope', "NOT_APPLICABLE"));
    summary.paper_entry_altitude_range_m = get_shape_field_local(shape, 'paper_entry_altitude_range_m', [NaN, NaN]);
    summary.paper_entry_fpa_range_deg = get_shape_field_local(shape, 'paper_entry_fpa_range_deg', [NaN, NaN]);
    summary.paper_capsule_mass_kg = get_shape_field_local(shape, 'paper_capsule_mass_kg', NaN);
    if isfinite(summary.paper_capsule_mass_kg)
        summary.paper_capsule_ballistic_coefficient_kg_m2 = ...
            summary.paper_capsule_mass_kg / (get_shape_field_local(shape, 'cd', NaN) * ...
            get_shape_field_local(shape, 'reference_area_m2', NaN));
    else
        summary.paper_capsule_ballistic_coefficient_kg_m2 = NaN;
    end
    summary.reference_total_heat_load_J_m2 = get_shape_field_local(shape, 'reference_total_heat_load_J_m2', NaN);
    summary.surrogate_exceeds_reference_total_heat_load = ...
        isfinite(summary.reference_total_heat_load_J_m2) && ...
        summary.total_heat_load_J_m2 > summary.reference_total_heat_load_J_m2;
    if summary.heat_constraint_comparable
        summary.exceeds_reference_total_heat_load = ...
            summary.surrogate_exceeds_reference_total_heat_load;
    else
        summary.exceeds_reference_total_heat_load = NaN;
    end
    summary.density_uncertainty_fraction = get_shape_field_local(shape, 'density_uncertainty_fraction', NaN);
    summary.explicit_cd_uncertainty_fraction = get_shape_field_local(shape, 'explicit_cd_uncertainty_fraction', NaN);
    summary.cd_uncertainty_fraction = get_shape_field_local(shape, 'cd_uncertainty_fraction', NaN);
    summary.ld_recession_uncertainty_fraction = get_shape_field_local(shape, 'ld_recession_uncertainty_fraction', NaN);
    summary.ld_bounds = get_shape_field_local(shape, 'ld_bounds', [NaN, NaN]);
    summary.density_scale = get_shape_field_local(shape, 'density_scale', 1.0);
    summary.cd_scale = get_shape_field_local(shape, 'cd_scale', 1.0);
    summary.ld_scale = get_shape_field_local(shape, 'ld_scale', 1.0);
end

function [X_cross, t_cross] = refine_reentry_altitude_crossing_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, target_altitude, dt_window)
    lo = 0;
    hi = dt_window;
    X_cross = rk4_reentry_step_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, hi);

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        X_mid = rk4_reentry_step_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, mid);
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

function [X_cross, t_cross] = refine_reentry_speed_crossing_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, target_speed, dt_window)
    lo = 0;
    hi = dt_window;
    X_cross = rk4_reentry_step_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, hi);

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        X_mid = rk4_reentry_step_local(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, mid);
        aux_mid = reentry_aux_local(X_mid, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);

        if aux_mid.speed_rel > target_speed
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
    rho = rho * get_shape_field_local(shape, 'density_scale', 1.0);

    if get_bool_param_local(atmosphere_opts, 'co_rotate_atmosphere', true)
        omega = get_atmos_field_local(atmosphere_opts, 'earth_rotation_rad_s', 7.2921159e-5);
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
    [aoa_actual_deg, bank_actual_deg] = resolve_entry_commands_local( ...
        shape, aoa_deg, bank_angle_deg, speed_rel, mach);

    a_gravity = gravity_acceleration_local(r, sys, shape);
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
        [cd, cl, ld] = aerodynamic_coefficients_local(shape, aoa_actual_deg, mach);
        area = get_shape_field_local(shape, 'reference_area_m2', 1.0);
        drag_accel_mag = dynamic_pressure * cd * area / mass;
        a_drag = -drag_accel_mag * v_rel_hat;

        if lift_enabled
            lift_dir = lift_direction_local(r, v_rel, bank_actual_deg);
            a_lift = ld * drag_accel_mag * lift_dir;
        else
            cl = 0;
            ld = 0;
        end

        nose_radius = max(get_shape_field_local(shape, 'nose_radius_m', 0.05), 1e-4);
        heat_flux = heat_k * sqrt(rho / nose_radius) * speed_rel^3;
        a_aero = a_drag + a_lift;
    end

    aux = struct();
    aux.altitude = altitude;
    aux.rho = rho;
    aux.v_rel = v_rel;
    aux.speed_rel = speed_rel;
    aux.fpa_deg = flight_path_angle_deg_local(r, v);
    aux.fpa_rel_deg = flight_path_angle_deg_local(r, v_rel);
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
    guidance_threshold = get_shape_field_local(shape, 'guidance_activation_drag_g', inf);
    aux.capsule_guidance_active = aux.drag_g >= guidance_threshold;
    guidance_cutoff_mach = get_shape_field_local(shape, 'guidance_cutoff_mach', -inf);
    aux.capsule_guidance_window_active = aux.capsule_guidance_active && aux.mach >= guidance_cutoff_mach;
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

function a = gravity_acceleration_local(r, sys, shape)
    r_norm = norm(r);
    a_g = -sys.mu / r_norm^3 * r;
    gravity_model = upper(string(get_shape_field_local(shape, 'gravity_model', "CENTRAL_SPHERICAL")));
    if gravity_model == "CENTRAL_SPHERICAL" || gravity_model == "CENTRAL"
        a = a_g;
        return;
    end

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

function [blackout_active, constraint_active, phase] = ...
        antenna_constraint_activity_local(hist, aux, shape)
    blackout_active = false;
    constraint_active = false;
    phase = 0;
    if ~get_bool_param_local(shape, 'antenna_enabled', false)
        return;
    end

    scope = upper(string(get_shape_field_local(shape, ...
        'antenna_tracking_scope', "CONTINUOUS")));
    lower_alt = get_shape_field_local(shape, ...
        'blackout_lower_altitude_m', 60e3);
    upper_alt = get_shape_field_local(shape, ...
        'blackout_upper_altitude_m', 80e3);

    if scope == "ALTITUDE_BAND"
        blackout_active = aux.altitude >= lower_alt && ...
            aux.altitude <= upper_alt;
        phase = double(blackout_active);
    else
        % Zhang Eq. (34) starts at the paper's 80 km initial point and
        % remains active until the first downward 60 km crossing, including
        % any intervening skip above 80 km. Mission trajectories begin at
        % 120 km, so phase 0 explicitly waits for the first downward entry
        % through 80 km instead of treating all h >= 60 km as blackout.
        initial_sample = isempty(hist.blackout_zone_phase);
        if ~initial_sample
            phase = hist.blackout_zone_phase(end);
        end
        descending = aux.fpa_rel_deg <= 0;
        if initial_sample && aux.altitude < lower_alt
            phase = 2;
        elseif initial_sample && aux.altitude <= upper_alt
            phase = 1;
        elseif phase == 0 && descending && ...
                aux.altitude <= upper_alt && aux.altitude >= lower_alt
            phase = 1;
        elseif phase == 0 && descending && aux.altitude < lower_alt
            % A coarse integration step may skip the whole sampled band.
            % No BZC sample can be reconstructed, but terminal phase 2
            % prevents a later spurious reactivation.
            phase = 2;
        elseif phase == 1 && descending && aux.altitude < lower_alt
            phase = 2;
        end
        blackout_active = phase == 1;
    end

    if scope == "BLACKOUT_ONLY" || scope == "PAPER_BZC" || ...
            scope == "ALTITUDE_BAND"
        constraint_active = blackout_active;
    else
        constraint_active = true;
    end
end

function antenna = antenna_geometry_local(X_rv, X_relay, aux, shape, los_clear)
    antenna = struct('raap_deg', NaN, 'range_m', NaN, 'constraint_active', false, ...
                     'blackout_active', false, ...
                     'in_beam', false, 'within_cone', false, ...
                     'within_range', false, 'communication_available', false, ...
                     'best_raap_deg', NaN, ...
                     'best_bank_deg', NaN);
    if ~get_bool_param_local(shape, 'antenna_enabled', false)
        return;
    end

    if isempty(X_relay) || numel(X_relay) < 3 || aux.speed_rel <= eps
        return;
    end

    boresight = get_shape_field_local(shape, 'antenna_boresight_body', [-1;0;0]);
    antenna.range_m = norm(X_relay(1:3) - X_rv(1:3));
    antenna.raap_deg = raap_deg_local(X_rv(1:3), aux.v_rel, X_relay(1:3), ...
                                      aux.aoa_deg, aux.bank_angle_deg, boresight);

    beam_half_angle = get_shape_field_local(shape, 'antenna_beam_half_angle_deg', 45);
    angle_tolerance = get_shape_field_local(shape, 'antenna_angle_tolerance_deg', 1e-10);
    min_range = get_shape_field_local(shape, 'antenna_min_range_m', 0);
    max_range = get_shape_field_local(shape, 'antenna_max_range_m', inf);
    antenna.within_cone = antenna.raap_deg <= beam_half_angle + angle_tolerance;
    antenna.in_beam = antenna.within_cone;
    antenna.within_range = antenna.range_m >= min_range && antenna.range_m <= max_range;
    antenna.communication_available = antenna.within_cone && antenna.within_range && ...
                                      isfinite(los_clear) && los_clear ~= 0;

    bank_scan_enabled = get_bool_param_local(shape, ...
        'evaluate_bank_feasibility', false);
    if bank_scan_enabled
        constraints = get_shape_field_local(shape, 'constraints', struct());
        max_bank = get_struct_field_local(constraints, 'max_bank_angle_deg', 60);
        bank_grid = linspace(-max_bank, max_bank, max(3, ceil(2*max_bank) + 1));
        raap_grid = nan(size(bank_grid));
        for ii = 1:numel(bank_grid)
            raap_grid(ii) = raap_deg_local(X_rv(1:3), aux.v_rel, X_relay(1:3), ...
                                           aux.aoa_deg, bank_grid(ii), boresight);
        end
        [antenna.best_raap_deg, idx] = min(raap_grid);
        antenna.best_bank_deg = bank_grid(idx);
    else
        antenna.best_raap_deg = NaN;
        antenna.best_bank_deg = NaN;
    end
    if bank_scan_enabled && (~isfinite(los_clear) || los_clear == 0)
        antenna.best_raap_deg = inf;
    end
end

function raap_deg = raap_deg_local(r, v_rel, r_relay, aoa_deg, bank_deg, boresight_body)
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
              -sin(alpha), cos(alpha), 0; ...
               0, 0, 1];
    g_bank = [1, 0, 0; ...
              0, cos(sigma), sin(sigma); ...
              0, -sin(sigma), cos(sigma)];
    w_body = g_alpha * g_bank * w_t;

    boresight = boresight_body(:);
    cosine = dot(boresight, w_body) / (norm(boresight) * norm(w_body));
    raap_deg = rad2deg(acos(min(1, max(-1, cosine))));
end

function satisfied = path_constraints_satisfied_local(aux, shape)
    constraints = get_shape_field_local(shape, 'constraints', struct());
    max_q = get_struct_field_local(constraints, 'max_dynamic_pressure_Pa', inf);
    max_g = get_struct_field_local(constraints, 'max_g_load', inf);
    max_heat = get_struct_field_local(constraints, 'max_heat_flux_W_m2', inf);
    max_bank = get_struct_field_local(constraints, 'max_bank_angle_deg', inf);
    heat_comparable = get_bool_param_local(shape, 'heat_constraint_comparable', true);

    satisfied = aux.dynamic_pressure <= max_q && ...
                aux.g_load <= max_g && ...
                (~heat_comparable || aux.heat_flux <= max_heat) && ...
                abs(aux.bank_angle_deg) <= max_bank;
end

function rate = max_rate_local(values, time)
    if numel(values) < 2
        rate = 0;
        return;
    end
    dt = diff(time);
    valid = dt > 0;
    if ~any(valid)
        rate = 0;
    else
        dv = abs(diff(values));
        rate = max(dv(valid) ./ dt(valid));
    end
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

function X = paper_tdrs_state_local(sys, shape, t_abs)
    radius = get_shape_field_local(shape, 'paper_tdrs_geocentric_radius_m', 42164e3);
    longitude = deg2rad(get_shape_field_local(shape, 'paper_tdrs_longitude_deg', 77));
    latitude = deg2rad(get_shape_field_local(shape, 'paper_tdrs_latitude_deg', 0));
    atmosphere_opts = get_atmosphere_opts_local(sys);
    omega = get_atmos_field_local(atmosphere_opts, 'earth_rotation_rad_s', 7.2921159e-5);
    if ~isscalar(omega)
        omega = omega(3);
    end

    inertial_longitude = longitude + omega * t_abs;
    cos_latitude = cos(latitude);
    r = radius * [ ...
        cos_latitude * cos(inertial_longitude); ...
        cos_latitude * sin(inertial_longitude); ...
        sin(latitude)];
    v = cross([0;0;omega], r);

    X = zeros(14,1);
    X(1:3) = r;
    X(4:6) = v;
    X(7:10) = [0;0;0;1];
    X(14) = get_struct_field_local(sys, 'Target_Mass', 1);
end

function X_next = rk4_chaser_step_local(X, sys, dt, t_abs)
    if isempty(X) || dt <= 0
        X_next = X;
        return;
    end

    k1 = Env_EOM(t_abs,          X,             [0;0;0], [0;0;0], sys, true, "target");
    k2 = Env_EOM(t_abs + dt/2,   X + k1*dt/2,   [0;0;0], [0;0;0], sys, true, "target");
    k3 = Env_EOM(t_abs + dt/2,   X + k2*dt/2,   [0;0;0], [0;0;0], sys, true, "target");
    k4 = Env_EOM(t_abs + dt,     X + k3*dt,     [0;0;0], [0;0;0], sys, true, "target");
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

    if norm(X_next(7:10)) > 0
        X_next(7:10) = X_next(7:10) / norm(X_next(7:10));
    end
end

function shape = get_reentry_shape_local(sys, params)
    vehicle_mode = upper(get_string_param_local(params, 'vehicle_mode', ...
        get_reentry_vehicle_field_local(sys, 'vehicle_mode', "SPACEPLANE")));
    if vehicle_mode ~= "SPACEPLANE" && vehicle_mode ~= "CAPSULE"
        error('Unknown re-entry vehicle_mode: %s. Use SPACEPLANE or CAPSULE.', char(vehicle_mode));
    end

    default_shape = get_reentry_vehicle_field_local(sys, 'selected_shape', "COMPROMISE");
    if vehicle_mode == "CAPSULE"
        name = "CAPSULE";
    else
        name = get_string_param_local(params, 'shape_name', default_shape);
    end
    key = normalize_shape_key_local(name);

    if ~isfield(sys, 'reentry_vehicle') || ~isfield(sys.reentry_vehicle, 'shapes') || ...
            ~isfield(sys.reentry_vehicle.shapes, key)
        error('Unknown reentry vehicle shape: %s.', char(name));
    end

    shape = sys.reentry_vehicle.shapes.(key);
    shape.name = string(get_shape_field_local(shape, 'name', key));
    shape.vehicle_mode = vehicle_mode;
    shape.gravity_model = upper(get_string_param_local(params, 'gravity_model', ...
        get_reentry_vehicle_field_local(sys, 'gravity_model', "CENTRAL_SPHERICAL")));
    if shape.gravity_model ~= "CENTRAL_SPHERICAL" && shape.gravity_model ~= "CENTRAL" && ...
            shape.gravity_model ~= "J2"
        error('Unknown re-entry gravity_model: %s. Use CENTRAL_SPHERICAL or J2.', ...
              char(shape.gravity_model));
    end

    if vehicle_mode == "SPACEPLANE"
        cfg = sys.reentry_vehicle.spaceplane;
        shape.aero_model = get_struct_field_local(cfg, 'aero_model', "PAPER_RLV_POLYNOMIAL");
        shape.aoa_profile_mode = get_struct_field_local(cfg, 'aoa_profile_mode', "SPEED_SCHEDULE");
        shape.aoa_speed_grid_m_s = get_struct_field_local(cfg, 'aoa_speed_grid_m_s', [0 2000 5000 8000]);
        shape.aoa_values_deg = get_struct_field_local(cfg, 'aoa_values_deg', [15 15 40 40]);
        if numel(shape.aoa_speed_grid_m_s) ~= numel(shape.aoa_values_deg) || ...
                numel(shape.aoa_speed_grid_m_s) < 2 || ...
                any(~isfinite(shape.aoa_speed_grid_m_s)) || ...
                any(~isfinite(shape.aoa_values_deg)) || ...
                any(diff(shape.aoa_speed_grid_m_s) <= 0)
            error('SPACEPLANE AoA schedule requires finite, equal-length values on a strictly increasing speed grid.');
        end
        has_param_aoa_override = isfield(params, 'aoa_deg') && ~isempty(params.aoa_deg);
        has_sys_aoa_override = isfield(sys.reentry_vehicle, 'aoa_deg') && ...
                               ~isempty(sys.reentry_vehicle.aoa_deg);
        if has_param_aoa_override || has_sys_aoa_override
            shape.aoa_profile_mode = "CONSTANT_OVERRIDE";
        end
        shape.cl_polynomial = get_struct_field_local(cfg, 'cl_polynomial', [-0.041065 0.016292 0.0002602]);
        shape.cd_from_cl_polynomial = get_struct_field_local(cfg, 'cd_from_cl_polynomial', [0.080505 -0.03026 0.86495]);
        shape.heating_model = get_struct_field_local(cfg, 'heating_model', "SUTTON_GRAVES_SURROGATE");
        shape.paper_heating_model = get_struct_field_local(cfg, 'paper_heating_model', "ZHANG_EQ32_KQ_SQRT_RHO_V_3P15");
        shape.heat_constraint_comparable = get_bool_param_local(cfg, 'heat_constraint_comparable', false);
        shape.constraints = get_struct_field_local(cfg, 'constraints', struct());
        shape.paper_initial_altitude_m = get_struct_field_local(cfg, 'paper_initial_altitude_m', 80e3);
        shape.paper_initial_scenarios = get_struct_field_local(cfg, 'paper_initial_scenarios', nan(3,5));
        shape.paper_terminal_longitude_deg = get_struct_field_local(cfg, 'paper_terminal_longitude_deg', 111.3);
        shape.paper_terminal_latitude_deg = get_struct_field_local(cfg, 'paper_terminal_latitude_deg', 42.2);
        shape.paper_terminal_altitude_m = get_struct_field_local(cfg, 'paper_terminal_altitude_m', 25e3);
        shape.paper_terminal_speed_m_s = get_struct_field_local(cfg, 'paper_terminal_speed_m_s', 800);

        communication = get_struct_field_local(cfg, 'communication', struct());
        shape.antenna_enabled = get_bool_param_local(communication, 'enabled', true);
        shape.relay_mode = upper(string(get_struct_field_local( ...
            communication, 'relay_mode', "MISSION_TARGET_DYNAMIC_ORBIT")));
        if shape.relay_mode ~= "MISSION_TARGET_DYNAMIC_ORBIT" && ...
                shape.relay_mode ~= "PAPER_TDRS_STATIC_EARTH_FIXED"
            error(['Unknown relay_mode: %s. Use MISSION_TARGET_DYNAMIC_ORBIT ' ...
                   'or PAPER_TDRS_STATIC_EARTH_FIXED.'], char(shape.relay_mode));
        end
        shape.paper_tdrs_longitude_deg = get_struct_field_local(communication, 'paper_tdrs_longitude_deg', 77);
        shape.paper_tdrs_latitude_deg = get_struct_field_local(communication, 'paper_tdrs_latitude_deg', 0);
        shape.paper_tdrs_geocentric_radius_m = get_struct_field_local( ...
            communication, 'paper_tdrs_geocentric_radius_m', 42164e3);
        if ~isfinite(shape.paper_tdrs_longitude_deg) || ...
                ~isfinite(shape.paper_tdrs_latitude_deg) || ...
                abs(shape.paper_tdrs_latitude_deg) > 90 || ...
                ~isfinite(shape.paper_tdrs_geocentric_radius_m) || ...
                shape.paper_tdrs_geocentric_radius_m <= sys.Re
            error('Paper TDRS longitude/latitude/radius reference is invalid.');
        end
        antenna_mount = upper(string(get_struct_field_local(communication, 'antenna_mount', "AFT")));
        if antenna_mount == "PAPER_TOP"
            shape.antenna_boresight_body = get_struct_field_local(communication, 'paper_top_boresight_body', [0;1;0]);
        elseif antenna_mount == "AFT"
            shape.antenna_boresight_body = get_struct_field_local(communication, 'aft_boresight_body', [-1;0;0]);
        else
            error('Unknown antenna_mount: %s. Use AFT or PAPER_TOP.', char(antenna_mount));
        end
        shape.antenna_mount = antenna_mount;
        shape.antenna_beam_half_angle_deg = get_struct_field_local(communication, 'beam_half_angle_deg', 45);
        shape.antenna_angle_tolerance_deg = get_struct_field_local(communication, 'angle_tolerance_deg', 1e-10);
        if ~isfinite(shape.antenna_beam_half_angle_deg) || ...
                shape.antenna_beam_half_angle_deg < 0 || shape.antenna_beam_half_angle_deg > 180
            error('Antenna beam_half_angle_deg must be finite and within [0, 180].');
        end
        if ~isscalar(shape.antenna_angle_tolerance_deg) || ...
                ~isfinite(shape.antenna_angle_tolerance_deg) || shape.antenna_angle_tolerance_deg < 0
            error('Antenna angle_tolerance_deg must be a finite nonnegative scalar.');
        end
        if numel(shape.antenna_boresight_body) ~= 3 || ...
                any(~isfinite(shape.antenna_boresight_body)) || norm(shape.antenna_boresight_body) <= eps
            error('Antenna boresight must be a finite nonzero three-vector.');
        end
        shape.antenna_min_range_m = get_struct_field_local(communication, 'min_range_m', 0);
        shape.antenna_max_range_m = get_struct_field_local(communication, 'max_range_m', inf);
        if ~isscalar(shape.antenna_min_range_m) || ~isscalar(shape.antenna_max_range_m) || ...
                ~isfinite(shape.antenna_min_range_m) || isnan(shape.antenna_max_range_m) || ...
                shape.antenna_min_range_m < 0 || ...
                shape.antenna_max_range_m < shape.antenna_min_range_m
            error('Antenna range limits require 0 <= min_range_m <= max_range_m.');
        end
        shape.antenna_tracking_scope = upper(string(get_struct_field_local(communication, 'tracking_scope', "CONTINUOUS")));
        if shape.antenna_tracking_scope ~= "CONTINUOUS" && ...
                shape.antenna_tracking_scope ~= "BLACKOUT_ONLY" && ...
                shape.antenna_tracking_scope ~= "PAPER_BZC" && ...
                shape.antenna_tracking_scope ~= "ALTITUDE_BAND"
            error(['Unknown antenna tracking_scope: %s. Use CONTINUOUS, ' ...
                   'PAPER_BZC, BLACKOUT_ONLY, or ALTITUDE_BAND.'], ...
                  char(shape.antenna_tracking_scope));
        end
        shape.blackout_upper_altitude_m = get_struct_field_local(communication, 'blackout_upper_altitude_m', 80e3);
        shape.blackout_lower_altitude_m = get_struct_field_local(communication, 'blackout_lower_altitude_m', 60e3);
        if shape.blackout_lower_altitude_m > shape.blackout_upper_altitude_m
            error('blackout_lower_altitude_m must not exceed blackout_upper_altitude_m.');
        end
        shape.evaluate_bank_feasibility = get_bool_param_local(communication, 'evaluate_bank_feasibility', true);
    else
        cfg = sys.reentry_vehicle.capsule;
        shape.aero_model = get_struct_field_local(cfg, 'aero_model', "PAPER_CAPSULE_REDUCED");
        shape.aoa_profile_mode = "CONSTANT_TRIM";
        shape.default_aoa_deg = get_struct_field_local(cfg, 'trim_aoa_deg', 0);
        shape.nominal_ld = get_struct_field_local(cfg, 'nominal_ld', 0.25);
        shape.mass_kg = get_struct_field_local(cfg, 'mass_kg', 60);
        shape.altitude_termination_enabled = get_bool_param_local(cfg, 'altitude_termination_enabled', false);
        shape.safety_floor_altitude_m = get_struct_field_local(cfg, 'safety_floor_altitude_m', 0);
        shape.separation_mode = upper(get_string_param_local(params, 'separation_mode', ...
            get_struct_field_local(cfg, 'separation_mode', "ENTRY_INTERFACE")));
        shape.terminal_speed_m_s = get_struct_field_local(cfg, 'parachute_deploy_speed_m_s', 240);
        shape.guidance_activation_drag_g = get_struct_field_local(cfg, 'guidance_activation_drag_g', 0.20);
        shape.guidance_cutoff_mach = get_struct_field_local(cfg, 'guidance_cutoff_mach', 3);
        shape.heating_model = get_struct_field_local(cfg, 'heating_model', "SUTTON_GRAVES_SURROGATE");
        shape.paper_heating_model = get_struct_field_local(cfg, 'paper_heating_model', "DETRA_KEMP_RIDDELL_EQ28");
        shape.heat_constraint_comparable = get_bool_param_local(cfg, 'heat_constraint_comparable', false);
        shape.constraints = get_struct_field_local(cfg, 'constraints', struct());
        shape.reference_entry_interface_altitude_m = get_struct_field_local(cfg, 'entry_interface_altitude_m', 120e3);
        shape.reference_entry_speed_m_s = get_struct_field_local(cfg, 'reference_entry_speed_m_s', 7.97e3);
        shape.reference_entry_fpa_deg = get_struct_field_local(cfg, 'reference_entry_fpa_deg', -1.16);
        shape.density_uncertainty_fraction = get_struct_field_local(cfg, 'density_uncertainty_fraction', 0.10);
        shape.explicit_cd_uncertainty_fraction = get_struct_field_local(cfg, 'explicit_cd_uncertainty_fraction', 0.10);
        shape.cd_uncertainty_fraction = get_struct_field_local(cfg, 'cd_uncertainty_fraction', 0.20);
        shape.ld_recession_uncertainty_fraction = get_struct_field_local(cfg, 'ld_recession_uncertainty_fraction', 0.10);
        shape.ld_bounds = get_struct_field_local(cfg, 'ld_bounds', [0.20, 0.30]);
        shape.paper_entry_altitude_range_m = get_struct_field_local(cfg, 'paper_entry_altitude_range_m', [121.01e3, 121.15e3]);
        shape.paper_entry_fpa_range_deg = get_struct_field_local(cfg, 'paper_entry_fpa_range_deg', [-1.17, -1.16]);
        shape.paper_entry_range_to_go_m = get_struct_field_local(cfg, 'paper_entry_range_to_go_m', [4852e3, 4882e3, 4890e3, 4920e3]);
        shape.paper_capsule_mass_kg = get_struct_field_local(cfg, 'paper_capsule_mass_kg', 150);
        shape.paper_entry_condition_scope = get_struct_field_local(cfg, 'paper_entry_condition_scope', "ALTITUDE_AND_FPA_ONLY");
        shape.reference_total_heat_load_J_m2 = get_struct_field_local(cfg, 'reference_total_heat_load_J_m2', 200e6);
        shape.antenna_enabled = false;
    end

    % Explicit study overrides support transparent ballistic-coefficient
    % sensitivity cases without modifying the source vehicle database.
    shape.reference_area_m2 = get_param_local(params, 'reference_area_m2', ...
        get_shape_field_local(shape, 'reference_area_m2', 1.0));
    shape.cd = get_param_local(params, 'cd', get_shape_field_local(shape, 'cd', 1.2));
    shape.nominal_ld = get_param_local(params, 'nominal_ld', ...
        get_shape_field_local(shape, 'nominal_ld', NaN));
    if ~isscalar(shape.reference_area_m2) || ~isfinite(shape.reference_area_m2) || ...
            shape.reference_area_m2 <= 0
        error('Re-entry reference_area_m2 must be a finite positive scalar.');
    end
    if ~isscalar(shape.cd) || ~isfinite(shape.cd) || shape.cd <= 0
        error('Re-entry Cd override must be a finite positive scalar.');
    end
    if vehicle_mode == "CAPSULE" && ...
            (~isscalar(shape.nominal_ld) || ~isfinite(shape.nominal_ld) || shape.nominal_ld < 0)
        error('CAPSULE nominal_ld must be a finite nonnegative scalar.');
    end

    uncertainty = get_struct_field_local(sys.reentry_vehicle, 'uncertainty', struct());
    shape.density_scale = get_param_local(params, 'density_scale', ...
        get_struct_field_local(uncertainty, 'density_scale', 1.0));
    shape.cd_scale = get_param_local(params, 'cd_scale', ...
        get_struct_field_local(uncertainty, 'cd_scale', 1.0));
    shape.ld_scale = get_param_local(params, 'ld_scale', ...
        get_struct_field_local(uncertainty, 'ld_scale', 1.0));
    if shape.density_scale <= 0 || shape.cd_scale <= 0 || shape.ld_scale < 0
        error('Re-entry density_scale and cd_scale must be positive; ld_scale must be nonnegative.');
    end
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

function [aoa_deg, bank_deg] = resolve_entry_commands_local(shape, aoa_fallback_deg, bank_fallback_deg, speed_rel, mach)
    profile_mode = upper(string(get_shape_field_local(shape, 'aoa_profile_mode', "CONSTANT")));
    bank_deg = bank_fallback_deg;

    if profile_mode == "SPEED_SCHEDULE"
        speed_grid = get_shape_field_local(shape, 'aoa_speed_grid_m_s', [0 8000]);
        aoa_grid = get_shape_field_local(shape, 'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        speed_query = min(max(speed_rel, speed_grid(1)), speed_grid(end));
        aoa_deg = interp1(speed_grid, aoa_grid, speed_query, 'linear');
    elseif profile_mode == "MACH_SCHEDULE"
        mach_grid = get_shape_field_local(shape, 'aoa_mach_grid', [0 30]);
        aoa_grid = get_shape_field_local(shape, 'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        mach_query = min(max(mach, mach_grid(1)), mach_grid(end));
        aoa_deg = interp1(mach_grid, aoa_grid, mach_query, 'linear');
    else
        aoa_deg = aoa_fallback_deg;
    end
end

function [cd, cl, ld] = aerodynamic_coefficients_local(shape, aoa_deg, mach)
    aero_model = upper(string(get_shape_field_local(shape, 'aero_model', "LEGACY_SHAPE_LD")));

    if aero_model == "PAPER_RLV_POLYNOMIAL"
        cl_coeff = get_shape_field_local(shape, 'cl_polynomial', [-0.041065, 0.016292, 0.0002602]);
        cd_coeff = get_shape_field_local(shape, 'cd_from_cl_polynomial', [0.080505, -0.03026, 0.86495]);
        cl = cl_coeff(1) + cl_coeff(2)*aoa_deg + cl_coeff(3)*aoa_deg^2;
        cd = cd_coeff(1) + cd_coeff(2)*cl + cd_coeff(3)*cl^2;
        cd = max(cd, 1e-6);
        ld = cl / cd;
    elseif aero_model == "PAPER_CAPSULE_REDUCED"
        cd = get_shape_field_local(shape, 'cd', 1.3);
        ld = get_shape_field_local(shape, 'nominal_ld', 0.25);
    else
        cd = get_shape_field_local(shape, 'cd', 1.2);
        ld = lookup_ld_local(shape, aoa_deg, mach);
    end

    cd = cd * get_shape_field_local(shape, 'cd_scale', 1.0);
    ld = ld * get_shape_field_local(shape, 'ld_scale', 1.0);
    cl = cd * ld;
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

function value = get_struct_field_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
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
