function [sys, provenance] = configure_system(options, cfg)
%CONFIGURE_SYSTEM Build the optional shared-core Zhang surrogate setup.
%   This function configures exact published AoA/aero/TDRS/BZC constants,
%   but vehicle, atmosphere and heat quantities remain labeled surrogates.

    if nargin < 2 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    if nargin < 1 || isempty(options)
        options = struct();
    end
    if ~isstruct(options)
        error('paperstudies:zhang:configureSystem:InvalidOptions', ...
              'options must be a struct.');
    end

    sys = Mission_Config();
    sys.Re = cfg.earth.radius_m;
    sys.reentry_vehicle.vehicle_mode = "SPACEPLANE";
    sys.reentry_vehicle.gravity_model = "CENTRAL_SPHERICAL";
    sys.reentry_vehicle.terminal_altitude = cfg.terminal.altitude_m;
    sys.reentry_vehicle.spaceplane.aero_model = "PAPER_RLV_POLYNOMIAL";
    sys.reentry_vehicle.spaceplane.aoa_profile_mode = "SPEED_SCHEDULE";
    sys.reentry_vehicle.spaceplane.aoa_speed_grid_m_s = [0, 2000, 5000, 8000];
    sys.reentry_vehicle.spaceplane.aoa_values_deg = [15, 15, 40, 40];
    sys.reentry_vehicle.spaceplane.cl_polynomial = cfg.aero.cl_polynomial;
    sys.reentry_vehicle.spaceplane.cd_from_cl_polynomial = ...
        cfg.aero.cd_from_cl_polynomial;
    sys.reentry_vehicle.spaceplane.paper_initial_altitude_m = ...
        cfg.table2.initial_altitude_m;
    sys.reentry_vehicle.spaceplane.paper_initial_scenarios = cfg.table2.values;
    sys.reentry_vehicle.spaceplane.paper_terminal_longitude_deg = ...
        cfg.terminal.longitude_deg;
    sys.reentry_vehicle.spaceplane.paper_terminal_latitude_deg = ...
        cfg.terminal.latitude_deg;
    sys.reentry_vehicle.spaceplane.paper_terminal_altitude_m = ...
        cfg.terminal.altitude_m;
    sys.reentry_vehicle.spaceplane.paper_terminal_speed_m_s = ...
        cfg.terminal.speed_m_s;
    sys.reentry_vehicle.spaceplane.constraints.max_dynamic_pressure_Pa = ...
        cfg.table3.max_dynamic_pressure_Pa;
    sys.reentry_vehicle.spaceplane.constraints.max_g_load = ...
        cfg.table3.max_load_factor;
    % Table 3 labels this value in W while Fig. 8 suggests a heat-flux
    % boundary with a different magnitude. Do not place the ambiguous raw
    % value into a W/m^2 comparison field.
    sys.reentry_vehicle.spaceplane.constraints.max_heat_flux_W_m2 = inf;
    sys.reentry_vehicle.spaceplane.constraints.max_bank_angle_deg = ...
        cfg.table3.max_bank_angle_deg;
    sys.reentry_vehicle.spaceplane.constraints.max_bank_rate_deg_s = ...
        cfg.table3.max_bank_rate_deg_s;
    sys.reentry_vehicle.spaceplane.heating_model = "SUTTON_GRAVES_SURROGATE";
    sys.reentry_vehicle.spaceplane.paper_heating_model = ...
        "ZHANG_EQ32_KQ_SQRT_RHO_V_3P15";
    sys.reentry_vehicle.spaceplane.heat_constraint_comparable = false;
    sys.reentry_vehicle.spaceplane.heat_constraint_required = true;

    comm = sys.reentry_vehicle.spaceplane.communication;
    comm.enabled = true;
    comm.relay_mode = cfg.communication.relay_mode;
    comm.antenna_mount = cfg.communication.antenna_mount;
    comm.paper_top_boresight_body = cfg.communication.antenna_boresight_body;
    comm.beam_half_angle_deg = cfg.table3.antenna_beam_half_angle_deg;
    comm.tracking_scope = cfg.communication.tracking_scope;
    comm.blackout_upper_altitude_m = ...
        cfg.communication.blackout_upper_altitude_m;
    comm.blackout_lower_altitude_m = ...
        cfg.communication.blackout_lower_altitude_m;
    comm.paper_tdrs_longitude_deg = cfg.table3.tdrs_longitude_deg;
    comm.paper_tdrs_latitude_deg = cfg.table3.tdrs_latitude_deg;
    comm.paper_tdrs_geocentric_radius_m = ...
        cfg.table3.tdrs_geocentric_radius_m;
    comm.earth_fixed_to_eci_angle_at_mission_epoch_deg = ...
        cfg.communication.earth_fixed_to_eci_angle_at_study_epoch_deg;
    comm.min_range_m = 0;
    comm.max_range_m = inf;
    comm.evaluate_bank_feasibility = false;
    sys.reentry_vehicle.spaceplane.communication = comm;

    sys.environment.atmospheric_drag.enabled = true;
    sys.environment.atmospheric_drag.model = cfg.surrogate.density_model;
    sys.environment.atmospheric_drag.use_matlab_atmosisa = false;
    sys.environment.atmospheric_drag.co_rotate_atmosphere = true;
    sys.environment.atmospheric_drag.earth_rotation_rad_s = ...
        cfg.surrogate.earth_rotation_rad_s;
    sys.reentry_vehicle.sutton_graves_k = cfg.surrogate.sutton_graves_k;

    provenance.classification = "SURROGATE_OPEN_LOOP_FORWARD";
    provenance.exact_settings = ["TABLE2", "TABLE3_NONTHERMAL", ...
        "AOA", "CL_CD", "PAPER_TDRS", "PAPER_TOP", "PAPER_BZC"];
    provenance.surrogate_settings = ["VEHICLE_MASS", "REFERENCE_AREA", ...
        "ISA76_DENSITY", "EARTH_ROTATION_NUMERIC", "SUTTON_GRAVES", ...
        "CONSTANT_BANK"];
    provenance.optimization_performed = false;
end
