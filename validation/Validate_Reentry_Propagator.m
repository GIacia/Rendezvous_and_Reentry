function results = Validate_Reentry_Propagator()
    % Lightweight deterministic checks for Reentry_Propagator.
    %
    % These tests verify numerical/event-handling behavior only. They do not
    % validate the assumed aerodynamic or aerothermal data.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(project_root);

    sys = Mission_Config();
    X0 = synthetic_entry_state_local(sys);

    base = struct();
    base.shape_name = "COMPROMISE";
    base.max_time = 300;
    base.terminal_altitude = 100e3;
    base.lift_enabled = false;

    coarse = base;
    coarse.dt = 2.0;
    [~, relay_empty, hist_coarse, summary_coarse] = ...
        Reentry_Propagator(sys, X0, [], 0, coarse);

    assert(isempty(relay_empty), 'An omitted relay state must remain empty.');
    assert(~summary_coarse.los_evaluated, 'LOS must be marked unevaluated when no relay is supplied.');
    assert(summary_coarse.reached_terminal_altitude, 'The descending test case must reach the terminal altitude.');
    assert(summary_coarse.termination_reason == "TERMINAL_ALTITUDE", 'Unexpected terminal event classification.');
    assert(abs(summary_coarse.terminal_altitude_error_m) < 1e-3, 'Terminal-altitude event was not refined accurately.');
    assert(isfield(hist_coarse, 'fpa_rel_deg'), 'Atmosphere-relative FPA history is missing.');
    assert(summary_coarse.gravity_model == "CENTRAL_SPHERICAL", ...
           'Paper-based re-entry must default to central spherical gravity.');
    assert(~summary_coarse.heat_constraint_evaluated && ...
           isnan(summary_coarse.path_constraints_satisfied), ...
           'A Sutton-Graves surrogate must not claim compliance with the paper heat law.');
    assert(~summary_coarse.antenna_constraint_fully_evaluated && ...
           isnan(summary_coarse.overall_constraints_satisfied), ...
           'Missing relay geometry must leave a required antenna constraint unevaluated.');

    violated_sys = sys;
    violated_sys.reentry_vehicle.spaceplane.constraints.max_dynamic_pressure_Pa = 0;
    [~, ~, ~, summary_violated] = ...
        Reentry_Propagator(violated_sys, X0, [], 0, coarse);
    assert(~summary_violated.path_constraints_satisfied && ...
           ~summary_violated.overall_constraints_satisfied, ...
           ['A definite evaluated-path violation must report false even when ' ...
            'heat and antenna constraints remain unevaluated.']);

    fine = base;
    fine.dt = 0.5;
    [~, ~, ~, summary_fine] = Reentry_Propagator(sys, X0, [], 0, fine);
    event_time_difference_s = abs(summary_coarse.duration_s - summary_fine.duration_s);
    assert(event_time_difference_s < 0.05, 'Terminal event time did not converge between dt=2.0 s and dt=0.5 s.');

    timeout = base;
    timeout.dt = 0.5;
    timeout.max_time = 1.0;
    timeout.terminal_altitude = 20e3;
    [~, ~, ~, summary_timeout] = Reentry_Propagator(sys, X0, [], 0, timeout);
    assert(~summary_timeout.reached_terminal_altitude, 'A one-second run must not report terminal-altitude completion.');
    assert(summary_timeout.termination_reason == "MAX_TIME", 'Timeout must be reported explicitly.');

    no_rotation_sys = sys;
    no_rotation_sys.environment.atmospheric_drag.co_rotate_atmosphere = false;
    initial_only = base;
    initial_only.dt = 1.0;
    initial_only.terminal_altitude = 130e3;
    [~, ~, hist_no_rotation, ~] = Reentry_Propagator(no_rotation_sys, X0, [], 0, initial_only);
    assert(abs(hist_no_rotation.speed_rel(1) - norm(X0(4:6))) < 1e-9, ...
           'Disabled atmosphere co-rotation must use inertial velocity as aerodynamic velocity.');

    expected_beta = X0(14) / ( ...
        summary_coarse.initial_cd * ...
        sys.reentry_vehicle.shapes.COMPROMISE.reference_area_m2);
    assert(abs(summary_coarse.ballistic_coefficient_kg_m2 - expected_beta) < 1e-10, ...
           'Reported ballistic coefficient is inconsistent with mass, Cd, and reference area.');
    assert(abs(hist_coarse.aoa_deg(1) - 40) < 1e-12, ...
           'The paper RLV speed-based AoA schedule must command 40 deg at entry speed.');

    constant_aoa_case = base;
    constant_aoa_case.aoa_deg = 23;
    constant_aoa_case.terminal_altitude = 130e3;
    [~, ~, hist_constant_aoa, ~] = ...
        Reentry_Propagator(sys, X0, [], 0, constant_aoa_case);
    assert(all(abs(hist_constant_aoa.aoa_deg - 23) < 1e-12), ...
           'An explicit scalar AoA must override the paper speed schedule.');

    zero_aoa_case = constant_aoa_case;
    zero_aoa_case.aoa_deg = 0;
    zero_aoa_case.lift_enabled = true;
    [~, ~, hist_zero_aoa, ~] = Reentry_Propagator(sys, X0, [], 0, zero_aoa_case);
    assert(abs(hist_zero_aoa.cl(1) - sys.reentry_vehicle.spaceplane.cl_polynomial(1)) < 1e-12, ...
           'The paper CL polynomial must retain its negative zero-AoA value without clipping.');

    endpoint_sys = sys;
    endpoint_sys.reentry_vehicle.spaceplane.aoa_speed_grid_m_s = [0, 2000, 5000];
    endpoint_sys.reentry_vehicle.spaceplane.aoa_values_deg = [10, 40, 20];
    [~, ~, hist_endpoint, ~] = ...
        Reentry_Propagator(endpoint_sys, X0, [], 0, initial_only);
    assert(abs(hist_endpoint.aoa_deg(1) - 20) < 1e-12, ...
           'Out-of-range speed must clamp to the schedule endpoint value.');

    relay = synthetic_relay_state_local(sys);
    antenna_case = base;
    antenna_case.dt = 1;
    antenna_case.terminal_altitude = 130e3;
    [~, ~, hist_antenna, summary_antenna] = Reentry_Propagator(sys, X0, relay, 0, antenna_case);
    assert(summary_antenna.antenna_evaluated, 'SPACEPLANE antenna geometry was not evaluated.');
    assert(isfinite(hist_antenna.raap_deg(1)), 'SPACEPLANE RAAP must be finite with a relay state.');
    assert(isfinite(hist_antenna.best_feasible_raap_deg(1)), 'Bank-feasibility RAAP must be reported.');
    assert(isequaln(hist_antenna.relay_pos, hist_antenna.chaser_pos), ...
           'The relay history and deprecated chaser alias must remain compatible.');

    no_scan_sys = sys;
    no_scan_sys.reentry_vehicle.spaceplane.communication.evaluate_bank_feasibility = false;
    [~, ~, hist_no_scan, summary_no_scan] = ...
        Reentry_Propagator(no_scan_sys, X0, relay, 0, antenna_case);
    assert(~summary_no_scan.antenna_bank_scan_enabled && ...
           all(isnan(hist_no_scan.best_feasible_raap_deg)) && ...
           isnan(summary_no_scan.instantaneous_bank_geometric_reachable), ...
           'A disabled bank scan must remain explicitly unevaluated.');

    body_relay_direction = relay_direction_body_local(sys, X0, relay, 40, 0);
    aligned_sys = sys;
    aligned_sys.reentry_vehicle.spaceplane.communication.aft_boresight_body = body_relay_direction;
    [~, ~, hist_aligned, ~] = Reentry_Propagator(aligned_sys, X0, relay, 0, antenna_case);
    assert(abs(hist_aligned.raap_deg(1)) < 1e-9, ...
           'A boresight aligned with the relay must produce zero RAAP.');
    assert(hist_aligned.communication_geometry_available(1), ...
           'Aligned in-range Earth-clear antenna geometry must be available.');

    basis = [1;0;0];
    if abs(dot(basis, body_relay_direction)) > 0.9
        basis = [0;1;0];
    end
    perpendicular = basis - dot(basis, body_relay_direction) * body_relay_direction;
    perpendicular = perpendicular / norm(perpendicular);
    boundary_boresight = cosd(45) * body_relay_direction + sind(45) * perpendicular;
    boundary_sys = sys;
    boundary_sys.reentry_vehicle.spaceplane.communication.aft_boresight_body = boundary_boresight;
    [~, ~, hist_boundary, ~] = Reentry_Propagator(boundary_sys, X0, relay, 0, antenna_case);
    assert(abs(hist_boundary.raap_deg(1) - 45) < 1e-9, ...
           'Constructed antenna geometry must lie on the 45 deg cone boundary.');
    assert(hist_boundary.antenna_within_cone(1), ...
           'The exact antenna cone boundary must be accepted.');

    range_limited_sys = sys;
    range_limited_sys.reentry_vehicle.spaceplane.communication.max_range_m = 1;
    [~, ~, hist_range_limited, summary_range_limited] = ...
        Reentry_Propagator(range_limited_sys, X0, relay, 0, antenna_case);
    assert(~hist_range_limited.antenna_within_range(1), ...
           'A one-metre maximum range must reject the synthetic relay.');
    assert(~hist_range_limited.communication_geometry_available(1), ...
           'Communication geometry must include the configured range limit.');
    assert(~summary_range_limited.antenna_tracking_maintained, ...
           'Tracking must not be reported as maintained after a range failure.');

    partial_sys = sys;
    partial_sys.reentry_vehicle.spaceplane.communication.max_range_m = 1;
    partial_relay = X0;
    partial_relay(4:6) = -X0(4:6);
    partial_relay(14) = sys.Target_Mass;
    partial_case = antenna_case;
    partial_case.dt = 1;
    partial_case.max_time = 10;
    partial_case.terminal_altitude = 119e3;
    [~, ~, ~, summary_partial_antenna] = ...
        Reentry_Propagator(partial_sys, X0, partial_relay, 0, partial_case);
    assert(~summary_partial_antenna.antenna_constraint_fully_evaluated && ...
           summary_partial_antenna.antenna_constraint_any_evaluated_violation && ...
           ~summary_partial_antenna.antenna_tracking_maintained && ...
           ~summary_partial_antenna.overall_constraints_satisfied, ...
           ['A known antenna failure must remain false when another required ' ...
            'sample has unavailable geometry.']);

    occulted_relay = relay;
    occulted_relay(1:3) = -relay(1:3);
    [~, ~, hist_occulted, ~] = Reentry_Propagator(sys, X0, occulted_relay, 0, antenna_case);
    assert(~hist_occulted.los_clear(1), ...
           'An Earth-occulted relay must fail the line-of-sight test.');
    assert(~hist_occulted.communication_geometry_available(1), ...
           'Earth occultation must make communication geometry unavailable.');

    blackout_only_sys = sys;
    blackout_only_sys.reentry_vehicle.spaceplane.communication.tracking_scope = "PAPER_BZC";
    [~, ~, hist_above_upper_blackout, ~] = ...
        Reentry_Propagator(blackout_only_sys, X0, relay, 0, antenna_case);
    assert(~hist_above_upper_blackout.blackout_zone_active(1), ...
           'A mission trajectory at 120 km must wait for downward entry through 80 km.');

    skip_entry_state = X0;
    skip_entry_state(1:3) = (sys.Re + 80e3) * X0(1:3) / norm(X0(1:3));
    skip_speed = 7800;
    skip_fpa = deg2rad(-0.1);
    skip_entry_state(4:6) = [skip_speed*sin(skip_fpa); ...
                             skip_speed*cos(skip_fpa); 0];
    skip_case = antenna_case;
    skip_case.dt = 0.25;
    skip_case.terminal_altitude = 50e3;
    skip_case.max_time = 120;
    skip_case.lift_enabled = true;
    [~, ~, hist_skip_blackout, ~] = ...
        Reentry_Propagator(blackout_only_sys, skip_entry_state, relay, 0, skip_case);
    skip_above_idx = find(hist_skip_blackout.altitude > 80e3 & ...
                          hist_skip_blackout.time > 0, 1, 'first');
    assert(hist_skip_blackout.blackout_zone_active(1), ...
           'Paper BZC must activate at the descending 80 km entry point.');
    assert(~isempty(skip_above_idx) && ...
           hist_skip_blackout.blackout_zone_active(skip_above_idx), ...
           'Paper BZC must remain latched during a post-entry skip above 80 km.');

    inside_ascending_state = X0;
    inside_ascending_state(1:3) = ...
        (sys.Re + 70e3) * X0(1:3) / norm(X0(1:3));
    inside_speed = 7800;
    inside_fpa = deg2rad(5);
    inside_ascending_state(4:6) = [inside_speed*sin(inside_fpa); ...
                                   inside_speed*cos(inside_fpa); 0];
    inside_case = skip_case;
    inside_case.max_time = 1;
    [~, ~, hist_inside_ascending, ~] = ...
        Reentry_Propagator(blackout_only_sys, inside_ascending_state, relay, 0, inside_case);
    assert(hist_inside_ascending.blackout_zone_phase(1) == 1 && ...
           hist_inside_ascending.blackout_zone_active(1), ...
           'An initial state inside the paper BZC must start active despite unknown prehistory.');

    exit_entry_state = X0;
    exit_entry_state(1:3) = (sys.Re + 80e3) * X0(1:3) / norm(X0(1:3));
    exit_case = skip_case;
    exit_case.dt = 1;
    exit_case.max_time = 300;
    exit_case.lift_enabled = false;
    [~, ~, hist_blackout_exit, ~] = ...
        Reentry_Propagator(blackout_only_sys, exit_entry_state, relay, 0, exit_case);
    exited_idx = find(hist_blackout_exit.blackout_zone_phase == 2, 1, 'first');
    assert(~isempty(exited_idx), ...
           'Paper BZC must enter its terminal phase after the downward 60 km crossing.');
    assert(all(hist_blackout_exit.blackout_zone_phase(exited_idx:end) == 2) && ...
           ~any(hist_blackout_exit.blackout_zone_active(exited_idx:end)), ...
           'Paper BZC must remain inactive after its first downward 60 km exit.');

    below_blackout_state = X0;
    below_blackout_state(1:3) = (sys.Re + 55e3) * X0(1:3) / norm(X0(1:3));
    below_blackout_case = antenna_case;
    below_blackout_case.terminal_altitude = 50e3;
    below_blackout_case.max_time = 1;
    [~, ~, hist_no_blackout, summary_no_blackout] = ...
        Reentry_Propagator(blackout_only_sys, below_blackout_state, relay, 0, below_blackout_case);
    assert(hist_no_blackout.blackout_zone_phase(1) == 2 && ...
           ~hist_no_blackout.blackout_zone_active(1), ...
           'An initial state below 60 km must start in the terminal BZC phase.');
    assert(~summary_no_blackout.bzc_constraint_evaluated, ...
           'A trajectory entirely below the paper 60 km BZC boundary must be unevaluated.');
    assert(isnan(summary_no_blackout.bzc_constraint_maintained), ...
           'An unevaluated blackout constraint must not report a vacuous success.');

    paper_tdrs_sys = sys;
    paper_tdrs_sys.reentry_vehicle.spaceplane.communication.relay_mode = ...
        "PAPER_TDRS_STATIC_EARTH_FIXED";
    [~, ~, hist_paper_tdrs, summary_paper_tdrs] = ...
        Reentry_Propagator(paper_tdrs_sys, X0, [], 123, antenna_case);
    assert(summary_paper_tdrs.relay_mode == "PAPER_TDRS_STATIC_EARTH_FIXED", ...
           'Paper TDRS relay mode was not selected.');
    assert(abs(norm(hist_paper_tdrs.relay_pos(:,1)) - 42164e3) < 1e-6, ...
           'Paper TDRS must use the published 42,164 km geocentric radius.');

    capsule_case = base;
    capsule_case.vehicle_mode = "CAPSULE";
    capsule_case.dt = 1;
    capsule_case.terminal_altitude = 100e3;
    capsule_case.lift_enabled = true;
    [X_capsule, ~, hist_capsule, summary_capsule] = ...
        Reentry_Propagator(sys, X0, [], 0, capsule_case);
    expected_capsule_beta = 60 / (1.3 * 0.554);
    assert(summary_capsule.vehicle_mode == "CAPSULE", 'CAPSULE mode selection failed.');
    assert(abs(X_capsule(14) - 60) < 1e-12, 'CAPSULE entry-interface separation must set the requested 60 kg mass.');
    assert(abs(summary_capsule.ballistic_coefficient_kg_m2 - expected_capsule_beta) < 1e-10, ...
           'CAPSULE ballistic coefficient is inconsistent with the paper area and Cd.');
    assert(all(abs(hist_capsule.ld - 0.25) < 1e-12), 'CAPSULE reduced model must use nominal L/D=0.25.');
    assert(summary_capsule.reference_entry_fpa_deg == -1.16, 'CAPSULE paper reference FPA is missing.');
    assert(~summary_capsule.altitude_termination_enabled, ...
           'CAPSULE mode must prioritize the paper 240 m/s termination over the generic altitude stop.');
    assert(isnan(summary_capsule.terminal_altitude_error_m), ...
           'Disabled capsule altitude termination must not report a terminal-altitude error.');
    activation_idx = find(hist_capsule.capsule_guidance_threshold_met, 1, 'first');
    if ~isempty(activation_idx)
        assert(all(hist_capsule.capsule_guidance_active(activation_idx:end)), ...
               'The 0.20 g aeroassist activation marker must latch after first crossing.');
    end
    assert(abs(summary_capsule.ld_recession_uncertainty_fraction - 0.10) < 1e-12, ...
           'The paper recession-style L/D uncertainty reference is missing.');

    low_mass_state = X0;
    low_mass_state(14) = 50;
    mass_error_thrown = false;
    try
        Reentry_Propagator(sys, low_mass_state, [], 0, capsule_case);
    catch err
        mass_error_thrown = contains(err.message, 'exceeds entry stack mass');
    end
    assert(mass_error_thrown, ...
           'CAPSULE separation must reject a capsule heavier than the entry stack.');

    results = struct();
    results.passed = true;
    results.event_time_difference_s = event_time_difference_s;
    results.terminal_altitude_error_m = summary_coarse.terminal_altitude_error_m;
    results.ballistic_coefficient_kg_m2 = summary_coarse.ballistic_coefficient_kg_m2;
    results.capsule_ballistic_coefficient_kg_m2 = summary_capsule.ballistic_coefficient_kg_m2;
    results.initial_raap_deg = hist_antenna.raap_deg(1);

    fprintf('Reentry validation passed.\n');
    fprintf('  terminal altitude error : %.6g m\n', results.terminal_altitude_error_m);
    fprintf('  dt event-time difference: %.6g s\n', results.event_time_difference_s);
    fprintf('  ballistic coefficient   : %.3f kg/m^2\n', results.ballistic_coefficient_kg_m2);
    fprintf('  capsule beta (60 kg)    : %.3f kg/m^2\n', results.capsule_ballistic_coefficient_kg_m2);
    fprintf('  initial antenna RAAP    : %.3f deg\n', results.initial_raap_deg);
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

function direction_body = relay_direction_body_local(sys, X_rv, X_relay, aoa_deg, bank_deg)
    r = X_rv(1:3);
    v = X_rv(4:6);
    if sys.environment.atmospheric_drag.co_rotate_atmosphere
        omega = sys.environment.atmospheric_drag.earth_rotation_rad_s;
        v = v - cross([0;0;omega], r);
    end

    x_t = v / norm(v);
    r_hat = r / norm(r);
    y_t = r_hat - dot(r_hat, x_t) * x_t;
    y_t = y_t / norm(y_t);
    z_t = cross(x_t, y_t);
    z_t = z_t / norm(z_t);
    w_t = [x_t, y_t, z_t]' * (X_relay(1:3) - r);

    alpha = deg2rad(aoa_deg);
    sigma = deg2rad(bank_deg);
    g_alpha = [cos(alpha), sin(alpha), 0; ...
              -sin(alpha), cos(alpha), 0; ...
               0, 0, 1];
    g_bank = [1, 0, 0; ...
              0, cos(sigma), sin(sigma); ...
              0, -sin(sigma), cos(sigma)];
    direction_body = g_alpha * g_bank * w_t;
    direction_body = direction_body / norm(direction_body);
end
