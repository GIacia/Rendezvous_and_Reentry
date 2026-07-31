function run = Mission_Run_Config()
    % Mission_Run_Config
    % Single user-facing control panel for normal mission runs.
    %
    % Typical workflow:
    %   1) edit this file
    %   2) run Main_Mission_Simulator
    %
    % Lower-level constants and vehicle shape tables remain in Mission_Config.m.
    defaults = Mission_Config();

    %% Runtime / external config selection
    % python_config.mode:
    %   "AUTO"    : select the newest matching configs/python_runs case
    %   "NONE"    : ignore Python optimizer JSON and use MATLAB defaults plus this file
    %   "FILE"    : load python_config.file
    %   "CASE_ID" : select one archived configs/python_runs case by case_id
    %   "HASH"    : select one archived case by settings_hash
    run.python_config.mode = "AUTO";
    run.python_config.file = "configs/latest_python_solution.json";
    run.python_config.case_id = "";
    run.python_config.hash = "";

    % Keep legacy setenv(...) overrides available for batch scripts. Set false
    % when you want this file to be the only run-control surface.
    run.runtime.allow_environment_overrides = true;

    %% Scenario
    run.scenario.h_insert_km = defaults.h_insert / 1e3;
    run.scenario.h_target_km = defaults.h_target / 1e3;
    run.scenario.h_wait_km = defaults.h_wait / 1e3;
    run.scenario.initial_chaser_angle_deg = 0;
    run.scenario.initial_phase_angle_deg = 90;

    %% Maneuver execution
    % "IMPULSIVE" or "FINITE_BURN".
    run.maneuver.burn_model = defaults.maneuver.default_burn_model;
    run.maneuver.finite_burn_thrust_N = defaults.maneuver.finite_burn_thrust;
    run.maneuver.finite_burn_isp_s = defaults.maneuver.finite_burn_isp;
    run.maneuver.finite_burn_dt_s = defaults.maneuver.finite_burn_dt;
    % No validated per-burn duration limit is available yet. Set a finite
    % value only when a real thermal/engine-cycle limit is being studied.
    run.maneuver.max_single_burn_duration_s = defaults.maneuver.max_single_burn_duration;
    run.maneuver.max_single_burn_delta_v_m_s = defaults.maneuver.max_single_burn_delta_v;
    run.maneuver.use_thrust_noise = false;
    run.maneuver.multi_hohmann_legs = [];

    %% Phase 1: phasing / homing
    run.phase1.mode = "CUSTOM_IMPULSE"; % "CUSTOM_IMPULSE", "HOHMANN"; "MULTI_HOHMANN" is preliminary

    % Leave [] to use the selected Python optimizer JSON values when available.
    run.phase1.phase_angle_deg = [];
    run.phase1.delta_v_m_s = [];
    run.phase1.gamma_deg = [];

    run.phase1.desired_rel_lvlh_m = [0; -5000; 0];
    run.phase1.time_step_s = 10;
    run.phase1.dt_phase_s = 5;
    run.phase1.dt_capture_s = 5;
    run.phase1.max_wait_s = [];
    run.phase1.max_capture_time_s = 30000;
    run.phase1.phase_tol_rad = 1e-7;
    run.phase1.event_time_tol_s = 1e-3;
    run.phase1.capture_pos_tol_m = 0.5;
    run.phase1.min_capture_time_s = 0;

    %% Phase 2: proximity operations
    run.phase2.dt_s = 1.0;
    run.phase2.S2_m = [0; -5000; 0];
    run.phase2.S4_R_abs_m = 30;
    run.phase2.initial_S2_tol_m = 50.0;
    run.phase2.tof_initial_s2_s = 1200;
    run.phase2.delta_R_cycloid_m = 400;
    run.phase2.vbar_burn_sign = -1;
    run.phase2.vbar_cross_tol_m = 1.0;
    run.phase2.max_cycloid_orbits = 4.0;
    run.phase2.rbar_hop_count = 8;
    run.phase2.tof_hop_s = 300;
    run.phase2.capture_pos_tol_m = 0.25;
    run.phase2.max_terminal_refines = 4;
    run.phase2.tof_terminal_refine_s = 180;
    run.phase2.Isp_fallback_s = 220;

    %% Phase 3: de-orbit / re-entry setup
    % "HOHMANN" or "R_BAR_200_FPA".
    run.phase3.mode = "HOHMANN";
    run.phase3.parking_altitude_km = defaults.h_reentry / 1e3;
    run.phase3.entry_interface_altitude_km = defaults.h_entry_interface / 1e3;
    run.phase3.flight_path_angle_deg = rad2deg(defaults.reentry_flight_path_angle);
    run.phase3.charge_final_reentry_fuel = false;
    run.phase3.target_radius_tol_m = 100;
    run.phase3.dt_reentry_coast_s = 2;
    run.phase3.max_reentry_coast_time_s = [];
    run.phase3.dt_rbar_wait_s = 30;
    run.phase3.max_rbar_wait_s = [];
    run.phase3.rbar_vbar_tol_m = 10;
    run.phase3.rbar_radial_tol_m = 50e3;
    % Drag-aware HOHMANN deorbit design:
    %   "AUTO"/"FILE": when orbital drag is enabled, load a Python-generated
    %                  single-retrograde-burn deorbit design from file.
    %   "OFF"        : keep the legacy dragless conic FPA injection.
    % Default DragDeorbitDesigner.py output is impulsive. Finite-burn deorbit
    % JSONs are experimental and must be generated explicitly.
    run.phase3.drag_deorbit_design.mode = "AUTO";
    run.phase3.drag_deorbit_design.file = "configs/latest_drag_deorbit_solution.json";
    run.phase3.drag_deorbit_design.delta_v_m_s = [];
    run.phase3.drag_deorbit_design.max_coast_time_s = [];

    %% Phase 4: atmospheric entry vehicle and diagnostics
    run.reentry.shape = defaults.reentry_vehicle.selected_shape; % COMPROMISE, HEATLOAD_MIN, PAYLOAD_MAX, TPS_MIN
    run.reentry.dt_s = defaults.reentry_vehicle.dt;
    run.reentry.max_time_s = defaults.reentry_vehicle.max_time;
    run.reentry.terminal_altitude_m = defaults.reentry_vehicle.terminal_altitude;
    run.reentry.lift_enabled = defaults.reentry_vehicle.lift_enabled;
    % Leave aoa_deg = [] to use the selected shape default. The current
    % entry model holds AoA and bank angle constant during propagation.
    run.reentry.aoa_deg = [];
    run.reentry.bank_angle_deg = defaults.reentry_vehicle.bank_angle_deg;
    run.reentry.los_margin_altitude_m = defaults.reentry_vehicle.los_margin_altitude;

    %% Environment
    % apply_from_phase:
    %   "PHASE1" : apply orbital drag from the beginning of the mission
    %   "PHASE3" : keep Phase 1/2 drag-free, then apply drag after berthing
    % The separated re-entry vehicle always uses the atmospheric entry model
    % in Reentry_Propagator.m.
    run.environment.atmospheric_drag.enabled = defaults.environment.atmospheric_drag.enabled;
    run.environment.atmospheric_drag.apply_from_phase = "PHASE3";
    run.environment.atmospheric_drag.model = defaults.environment.atmospheric_drag.model;
    run.environment.atmospheric_drag.use_matlab_atmosisa = defaults.environment.atmospheric_drag.use_matlab_atmosisa;
    run.environment.atmospheric_drag.co_rotate_atmosphere = defaults.environment.atmospheric_drag.co_rotate_atmosphere;
end
