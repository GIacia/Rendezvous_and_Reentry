function sys = Mission_Config()
    % Environment constants.
    sys.mu = 3.986004418e14;   % Earth gravitational parameter [m^3/s^2]
    sys.Re = 6378137;          % Earth equatorial radius [m]
    sys.J2 = 1.08263e-3;       % J2 coefficient
    % sys.J2 = 0;              % Set to zero for two-body debugging.
    sys.g0 = 9.80665;          % Standard gravity [m/s^2]

    % Spacecraft properties.
    sys.Target_Mass = 2000;    % [kg]
    sys.Chaser_Mass_Init = 2000; % Base chaser wet mass [kg], excluding an add-on capsule by default
    sys.Inertia = diag([800, 800, 600]); % Moment of inertia [kg*m^2]

    % Orbit and mission altitudes.
    sys.h_target = 500e3;      % 500 km
    sys.h_insert = 300e3;      % 300 km
    sys.h_wait   = 495e3;      % 495 km (Waiting Point)
    sys.h_reentry = 200e3;     % 200 km
    sys.h_entry_interface = 120e3; % 120 km atmospheric entry interface
    sys.reentry_flight_path_angle = 4 * pi/180; % 120 km interface FPA magnitude [rad]
    sys.inc = pi/2;            % Polar orbit inclination [rad]

    r_insert = sys.Re + sys.h_insert;
    r_wait = sys.Re + sys.h_wait;
    a_transfer = 0.5 * (r_insert + r_wait);
    sys.phase = pi * (1 - (a_transfer / r_wait)^1.5); % two-body Hohmann phase angle [rad]

    % Propulsion model.
    sys.Isp = 200;             % Specific impulse [s]
    sys.Thrust_Impulsive = 300;% Reference finite-burn thrust [N]

    % Maneuver execution options.
    % default_burn_model: "IMPULSIVE" preserves the baseline instantaneous delta-V behavior. Use "FINITE_BURN" only for explicit finite-duration firing studies.
    sys.maneuver.default_burn_model = "IMPULSIVE";
    sys.maneuver.direction_mode = "LOCAL_TANGENTIAL_RADIAL";
    sys.maneuver.finite_burn_thrust = sys.Thrust_Impulsive; % [N]
    sys.maneuver.finite_burn_isp = sys.Isp;                 % [s]
    sys.maneuver.finite_burn_dt = 0.1;                      % [s]
    sys.maneuver.max_single_burn_duration = inf;            % [s], disabled unless a validated thermal/engine limit is supplied
    sys.maneuver.max_single_burn_delta_v = inf;             % [m/s]

    % Environment model options.
    % atmospheric_drag.enabled = false preserves the previous J2-only propagation. Set enabled = true and model = "ISA76" to add drag.
    sys.environment.atmospheric_drag.enabled = false;
    sys.environment.atmospheric_drag.model = "ISA76";
    sys.environment.atmospheric_drag.use_matlab_atmosisa = true;
    sys.environment.atmospheric_drag.co_rotate_atmosphere = true;
    sys.environment.atmospheric_drag.earth_rotation_rad_s = 7.2921159e-5;
    sys.environment.atmospheric_drag.chaser_cd = 2.2;
    sys.environment.atmospheric_drag.chaser_area_m2 = 4.0;
    sys.environment.atmospheric_drag.target_cd = 2.2;
    sys.environment.atmospheric_drag.target_area_m2 = 4.0;

    % Atmospheric re-entry vehicle options.
    % Shape dimensions come from external re-entry vehicle data.
    % L/D values are approximate values from external data.
    sys.reentry_vehicle.vehicle_mode = "SPACEPLANE"; % SPACEPLANE or CAPSULE
    sys.reentry_vehicle.selected_shape = "COMPROMISE";
    sys.reentry_vehicle.dt = 0.5;                    % [s]
    sys.reentry_vehicle.max_time = 2500;             % [s]
    sys.reentry_vehicle.terminal_altitude = 20e3;    % [m]
    sys.reentry_vehicle.lift_enabled = true;
    sys.reentry_vehicle.gravity_model = "CENTRAL_SPHERICAL"; % CENTRAL_SPHERICAL or J2
    sys.reentry_vehicle.bank_angle_deg = 0;
    sys.reentry_vehicle.los_margin_altitude = 0;     % [m] Earth-limb clearance margin
    sys.reentry_vehicle.sutton_graves_k = 1.83e-4;   % Earth entry, W/m^2 with SI inputs
    sys.reentry_vehicle.uncertainty.density_scale = 1.0;
    sys.reentry_vehicle.uncertainty.cd_scale = 1.0;
    sys.reentry_vehicle.uncertainty.ld_scale = 1.0;

    % Paper-based reduced-order SPACEPLANE model.
    % Zhang et al., "Entry trajectory optimization considering blackout
    % zone communication constraint", Advances in Space Research (2026),
    % DOI: 10.1016/j.asr.2025.11.062.
    sys.reentry_vehicle.spaceplane.aero_model = "PAPER_RLV_POLYNOMIAL";
    sys.reentry_vehicle.spaceplane.aoa_profile_mode = "SPEED_SCHEDULE";
    sys.reentry_vehicle.spaceplane.aoa_speed_grid_m_s = [0, 2000, 5000, 8000];
    sys.reentry_vehicle.spaceplane.aoa_values_deg = [15, 15, 40, 40];
    sys.reentry_vehicle.spaceplane.cl_polynomial = [-0.041065, 0.016292, 0.0002602];
    sys.reentry_vehicle.spaceplane.cd_from_cl_polynomial = [0.080505, -0.03026, 0.86495];
    sys.reentry_vehicle.spaceplane.constraints.max_dynamic_pressure_Pa = 20e3;
    sys.reentry_vehicle.spaceplane.constraints.max_g_load = 2.5;
    sys.reentry_vehicle.spaceplane.constraints.max_heat_flux_W_m2 = 160e3;
    sys.reentry_vehicle.spaceplane.constraints.max_bank_angle_deg = 60;
    sys.reentry_vehicle.spaceplane.constraints.max_bank_rate_deg_s = 10;
    sys.reentry_vehicle.spaceplane.communication.enabled = true;
    % The paper uses [0;1;0] (upper-body boresight). AFT implements the
    % requested rear-facing antenna while retaining the same RAAP geometry.
    sys.reentry_vehicle.spaceplane.communication.antenna_mount = "AFT"; % AFT or PAPER_TOP
    sys.reentry_vehicle.spaceplane.communication.paper_top_boresight_body = [0;1;0];
    sys.reentry_vehicle.spaceplane.communication.aft_boresight_body = [-1;0;0];
    sys.reentry_vehicle.spaceplane.communication.beam_half_angle_deg = 45;
    sys.reentry_vehicle.spaceplane.communication.min_range_m = 0;
    sys.reentry_vehicle.spaceplane.communication.max_range_m = inf; % set from a real link budget when available
    sys.reentry_vehicle.spaceplane.communication.tracking_scope = "CONTINUOUS"; % CONTINUOUS or BLACKOUT_ONLY
    sys.reentry_vehicle.spaceplane.communication.blackout_upper_altitude_m = 80e3;
    sys.reentry_vehicle.spaceplane.communication.blackout_lower_altitude_m = 60e3;
    sys.reentry_vehicle.spaceplane.communication.evaluate_bank_feasibility = true;

    % Paper-based reduced-order CAPSULE model.
    % Saito et al., "Guidance strategies for controlled Earth reentry of
    % small spacecraft in low Earth orbit", Acta Astronautica 229 (2025),
    % DOI: 10.1016/j.actaastro.2024.12.054.
    % The requested 60 kg replaces the paper's 150 kg capsule mass. The
    % paper reference area and Cd are retained explicitly.
    sys.reentry_vehicle.capsule.mass_kg = 60;
    % Set false when maneuver.initial_mass_kg already represents the full
    % chaser-plus-capsule stack.
    sys.reentry_vehicle.capsule.add_to_chaser_initial_mass = true;
    sys.reentry_vehicle.capsule.separation_mode = "ENTRY_INTERFACE";
    sys.reentry_vehicle.capsule.use_paper_entry_conditions = true;
    sys.reentry_vehicle.capsule.aero_model = "PAPER_CAPSULE_REDUCED";
    sys.reentry_vehicle.capsule.nominal_ld = 0.25; % paper range: 0.20-0.30
    sys.reentry_vehicle.capsule.trim_aoa_deg = 0;  % surrogate: paper does not publish the Mach-10 trim angle
    sys.reentry_vehicle.capsule.entry_interface_altitude_m = 120e3;
    sys.reentry_vehicle.capsule.reference_entry_speed_m_s = 7.97e3;
    sys.reentry_vehicle.capsule.reference_entry_fpa_deg = -1.16;
    sys.reentry_vehicle.capsule.guidance_activation_drag_g = 0.20;
    sys.reentry_vehicle.capsule.guidance_cutoff_mach = 3;
    sys.reentry_vehicle.capsule.parachute_deploy_speed_m_s = 240;
    sys.reentry_vehicle.capsule.density_uncertainty_fraction = 0.10;
    sys.reentry_vehicle.capsule.cd_uncertainty_fraction = 0.20;
    sys.reentry_vehicle.capsule.ld_bounds = [0.20, 0.30];
    sys.reentry_vehicle.capsule.constraints.max_heat_flux_W_m2 = 1.0e6;
    sys.reentry_vehicle.capsule.reference_total_heat_load_J_m2 = 200e6;
    % Published reference vehicle values retained for provenance only.
    sys.reentry_vehicle.capsule.paper_capsule_mass_kg = 150;
    sys.reentry_vehicle.capsule.paper_spacecraft_total_mass_kg = 330;
    sys.reentry_vehicle.capsule.paper_satellite_reference_area_m2 = 0.640;
    sys.reentry_vehicle.capsule.paper_satellite_cd = 2.2;

    ld_aoa_deg = [0 2 5 8 10 12 15 18 20 22 25 30 35 40 45 50 55 60 65 70 75 80];

    sys.reentry_vehicle.shapes.COMPROMISE = struct( ...
        'name', "Compromise", ...
        'length_m', 4.486, ...
        'max_width_m', 1.944, ...
        'max_height_m', 1.296, ...
        'base_diameter_m', 1.753, ...
        'nose_radius_m', 0.054, ...
        'aspect_ratio', 2.559, ...
        'reference_area_m2', pi * 1.944 * 1.296 / 4, ...
        'cd', 1.20, ...
        'default_aoa_deg', 21, ...
        'ld_aoa_deg', ld_aoa_deg, ...
        'ld_values', [0.54 0.60 0.72 0.86 0.96 1.04 1.12 1.18 1.19 1.19 1.17 1.09 0.99 0.89 0.78 0.67 0.57 0.48 0.39 0.30 0.21 0.12]);

    sys.reentry_vehicle.shapes.HEATLOAD_MIN = struct( ...
        'name', "HeatLoad Min", ...
        'length_m', 4.485, ...
        'max_width_m', 2.137, ...
        'max_height_m', 1.339, ...
        'base_diameter_m', 1.809, ...
        'nose_radius_m', 0.054, ...
        'aspect_ratio', 2.480, ...
        'reference_area_m2', pi * 2.137 * 1.339 / 4, ...
        'cd', 1.20, ...
        'default_aoa_deg', 21, ...
        'ld_aoa_deg', ld_aoa_deg, ...
        'ld_values', [0.52 0.59 0.72 0.86 0.95 1.03 1.11 1.16 1.19 1.19 1.18 1.10 1.00 0.90 0.79 0.68 0.58 0.49 0.40 0.31 0.22 0.13]);

    sys.reentry_vehicle.shapes.PAYLOAD_MAX = struct( ...
        'name', "Payload Max", ...
        'length_m', 4.484, ...
        'max_width_m', 1.953, ...
        'max_height_m', 1.339, ...
        'base_diameter_m', 1.808, ...
        'nose_radius_m', 0.070, ...
        'aspect_ratio', 2.480, ...
        'reference_area_m2', pi * 1.953 * 1.339 / 4, ...
        'cd', 1.20, ...
        'default_aoa_deg', 22, ...
        'ld_aoa_deg', ld_aoa_deg, ...
        'ld_values', [0.42 0.48 0.58 0.70 0.80 0.90 1.02 1.07 1.09 1.09 1.09 1.04 0.96 0.87 0.77 0.67 0.57 0.48 0.38 0.29 0.20 0.13]);

    sys.reentry_vehicle.shapes.TPS_MIN = struct( ...
        'name', "TPS Min", ...
        'length_m', 4.477, ...
        'max_width_m', 1.944, ...
        'max_height_m', 1.011, ...
        'base_diameter_m', 1.321, ...
        'nose_radius_m', 0.052, ...
        'aspect_ratio', 3.390, ...
        'reference_area_m2', pi * 1.944 * 1.011 / 4, ...
        'cd', 1.20, ...
        'default_aoa_deg', 17, ...
        'ld_aoa_deg', ld_aoa_deg, ...
        'ld_values', [0.40 0.62 0.95 1.25 1.47 1.61 1.65 1.64 1.61 1.58 1.50 1.34 1.17 1.01 0.88 0.74 0.62 0.51 0.41 0.32 0.23 0.14]);

    sys.reentry_vehicle.shapes.CAPSULE = struct( ...
        'name', "Paper Capsule 60 kg", ...
        'length_m', 0.657, ...
        'max_width_m', 0.840, ...
        'max_height_m', 0.840, ...
        'base_diameter_m', 0.840, ...
        'nose_radius_m', 0.420, ... % explicit surrogate: half diameter; paper does not state R_n
        'reference_area_m2', 0.554, ...
        'cd', 1.30, ...
        'default_aoa_deg', 0, ...
        'ld_aoa_deg', [0 90], ...
        'ld_values', [0.25 0.25]);

    % Uncertainty placeholders. Active deterministic runs leave these off
    % unless a maneuver helper explicitly enables thrust noise.
    sys.noise.pos = 2.0;       % Position sensor noise [m]
    sys.noise.vel = 0.05;      % Velocity sensor noise [m/s]
    sys.noise.att = 0.001;     % Quaternion attitude noise
    sys.noise.gyro = 1e-4;     % Angular-rate noise [rad/s]
    sys.noise.thrust_err = 0.02; % Thrust magnitude error, 1-sigma fraction

    % Legacy controller gains retained for future 6-DOF controller work.
    sys.GNC.Kp_pos = 0.05; sys.GNC.Kd_pos = 0.8;
    sys.GNC.Kp_att = 50;   sys.GNC.Kd_att = 100;

    % Capture Point Settings
    sys.capture.r_rel0 = [-5000; 0; 0];   % 5 km R-bar hold point in LVLH
    sys.capture.v_rel0 = [0; 0; 0];       % Relative velocity at capture
    sys.capture.TOF = 1800;               % Nominal post-Phase-1 capture coast [s]
end
