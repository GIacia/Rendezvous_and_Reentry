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
        sys.reentry_vehicle.shapes.COMPROMISE.cd * ...
        sys.reentry_vehicle.shapes.COMPROMISE.reference_area_m2);
    assert(abs(summary_coarse.ballistic_coefficient_kg_m2 - expected_beta) < 1e-10, ...
           'Reported ballistic coefficient is inconsistent with mass, Cd, and reference area.');

    results = struct();
    results.passed = true;
    results.event_time_difference_s = event_time_difference_s;
    results.terminal_altitude_error_m = summary_coarse.terminal_altitude_error_m;
    results.ballistic_coefficient_kg_m2 = summary_coarse.ballistic_coefficient_kg_m2;

    fprintf('Reentry validation passed.\n');
    fprintf('  terminal altitude error : %.6g m\n', results.terminal_altitude_error_m);
    fprintf('  dt event-time difference: %.6g s\n', results.event_time_difference_s);
    fprintf('  ballistic coefficient   : %.3f kg/m^2\n', results.ballistic_coefficient_kg_m2);
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
