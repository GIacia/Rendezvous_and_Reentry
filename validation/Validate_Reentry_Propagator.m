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

    occulted_relay = relay;
    occulted_relay(1:3) = -relay(1:3);
    [~, ~, hist_occulted, ~] = Reentry_Propagator(sys, X0, occulted_relay, 0, antenna_case);
    assert(~hist_occulted.los_clear(1), ...
           'An Earth-occulted relay must fail the line-of-sight test.');
    assert(~hist_occulted.communication_geometry_available(1), ...
           'Earth occultation must make communication geometry unavailable.');

    blackout_only_sys = sys;
    blackout_only_sys.reentry_vehicle.spaceplane.communication.tracking_scope = "BLACKOUT_ONLY";
    [~, ~, ~, summary_no_blackout] = ...
        Reentry_Propagator(blackout_only_sys, X0, relay, 0, antenna_case);
    assert(~summary_no_blackout.bzc_constraint_evaluated, ...
           'A trajectory outside the blackout altitude band must be marked unevaluated.');
    assert(isnan(summary_no_blackout.bzc_constraint_maintained), ...
           'An unevaluated blackout constraint must not report a vacuous success.');

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
