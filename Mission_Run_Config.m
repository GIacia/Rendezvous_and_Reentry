function run = Mission_Run_Config()
    % Mission_Run_Config
    % Single user-facing control panel for normal mission runs.
    %
    % Typical workflow:
    %   1) edit this file
    %   2) run Main_Mission_Simulator
    %
    % Lower-level constants and vehicle shape tables remain in Mission_Config.m.

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
    run.scenario.h_insert_km = 300;
    run.scenario.h_target_km = 500;
    run.scenario.h_wait_km = 495;
    run.scenario.initial_chaser_angle_deg = 0;
    run.scenario.initial_phase_angle_deg = 90;

    %% Maneuver execution
    % "IMPULSIVE" or "FINITE_BURN".
    run.maneuver.burn_model = "FINITE_BURN";
    run.maneuver.finite_burn_thrust_N = 300;
    run.maneuver.finite_burn_isp_s = 200;
    run.maneuver.finite_burn_dt_s = 0.1;
    run.maneuver.max_single_burn_duration_s = 120;
    run.maneuver.max_single_burn_delta_v_m_s = inf;
    run.maneuver.use_thrust_noise = false;
    run.maneuver.multi_hohmann_legs = [];

    %% Phase 1: phasing / homing
    run.phase1.mode = "CUSTOM_IMPULSE"; % "CUSTOM_IMPULSE", "HOHMANN", "MULTI_HOHMANN"

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
    run.phase3.parking_altitude_km = 200;
    run.phase3.entry_interface_altitude_km = 120;
    run.phase3.flight_path_angle_deg = 4;
    run.phase3.charge_final_reentry_fuel = false;
    run.phase3.target_radius_tol_m = 100;
    run.phase3.dt_reentry_coast_s = 2;
    run.phase3.max_reentry_coast_time_s = [];
    run.phase3.dt_rbar_wait_s = 30;
    run.phase3.max_rbar_wait_s = [];
    run.phase3.rbar_vbar_tol_m = 10;
    run.phase3.rbar_radial_tol_m = 50e3;

    %% Phase 4: atmospheric entry vehicle and diagnostics
    run.reentry.shape = "COMPROMISE"; % COMPROMISE, HEATLOAD_MIN, PAYLOAD_MAX, TPS_MIN
    run.reentry.dt_s = 0.5;
    run.reentry.max_time_s = 2500;
    run.reentry.terminal_altitude_m = 20e3;
    run.reentry.lift_enabled = true;
    run.reentry.aoa_deg = [];
    run.reentry.bank_angle_deg = 0;
    run.reentry.los_margin_altitude_m = 0;

    %% Environment
    run.environment.atmospheric_drag.enabled = false;
    run.environment.atmospheric_drag.model = "ISA76";
    run.environment.atmospheric_drag.use_matlab_atmosisa = true;
    run.environment.atmospheric_drag.co_rotate_atmosphere = true;
end
