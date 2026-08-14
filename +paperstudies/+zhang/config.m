function cfg = config()
%CONFIG Published Zhang et al. inputs and explicitly labeled assumptions.
%   CFG = PAPERSTUDIES.ZHANG.CONFIG() returns only values visible in the
%   paper as exact data. Values needed solely by the optional shared-core
%   forward propagation live under CFG.SURROGATE.

    cfg.schema_version = "1.0";
    cfg.study_id = "ZHANG_BZC_2026";
    cfg.source.title = "Entry trajectory optimization considering blackout zone communication constraint";
    cfg.source.journal = "Advances in Space Research";
    cfg.source.volume = 77;
    cfg.source.year = 2026;
    cfg.source.pages = "11407-11417";
    cfg.source.doi = "10.1016/j.asr.2025.11.062";

    cfg.levels.exact = "EXACT";
    cfg.levels.convention = "CONVENTION_REQUIRED";
    cfg.levels.surrogate = "SURROGATE";
    cfg.levels.unavailable = "UNAVAILABLE";

    % Paper dynamics: spherical Earth and dimensionless Eqs. (15)-(25).
    cfg.earth.radius_m = 6378e3;
    cfg.earth.model = "CENTRAL_SPHERICAL";
    cfg.earth.rotation_rate_numeric_status = "UNAVAILABLE";

    % Fig. 2 and Eqs. (26)-(27). Alpha is evaluated in degrees.
    cfg.aero.aoa_break_speed_m_s = [2000, 5000];
    cfg.aero.aoa_plateau_deg = [15, 40];
    cfg.aero.figure_displayed_speed_domain_m_s = [1000, 7000];
    cfg.aero.cl_polynomial = [-0.041065, 0.016292, 0.0002602];
    cfg.aero.cd_from_cl_polynomial = [0.080505, -0.03026, 0.86495];
    cfg.aero.regression_aoa_deg = [15, 40];
    cfg.aero.regression_cl = [0.26186, 1.026935];
    cfg.aero.regression_cd = [0.131891308421020, 0.961602319629914];

    % Table 2. Columns: longitude, latitude, speed, FPA, heading.
    cfg.table2.columns = ["longitude_deg", "latitude_deg", ...
                          "speed_m_s", "fpa_deg", "heading_deg"];
    cfg.table2.initial_altitude_m = 80e3;
    cfg.table2.values = [ ...
        55, 25, 6800, -0.5, 45; ...
        65, 10, 7000, -1.0, 50; ...
        75,  0, 7200,  0.0, 50];

    cfg.cases = repmat(struct(), 1, size(cfg.table2.values, 1));
    for i = 1:numel(cfg.cases)
        row = cfg.table2.values(i,:);
        cfg.cases(i).id = i;
        cfg.cases(i).altitude_m = cfg.table2.initial_altitude_m;
        cfg.cases(i).longitude_deg = row(1);
        cfg.cases(i).latitude_deg = row(2);
        cfg.cases(i).speed_m_s = row(3);
        cfg.cases(i).fpa_deg = row(4);
        cfg.cases(i).heading_deg = row(5);
        cfg.cases(i).initial_bank_deg = NaN;
    end

    % Published terminal box.
    cfg.terminal.longitude_deg = 111.3;
    cfg.terminal.latitude_deg = 42.2;
    cfg.terminal.altitude_m = 25e3;
    cfg.terminal.speed_m_s = 800;
    cfg.terminal.longitude_tolerance_deg = 0.5;
    cfg.terminal.latitude_tolerance_deg = 0.5;
    cfg.terminal.altitude_tolerance_m = 0.5e3;
    cfg.terminal.speed_tolerance_m_s = 60;
    cfg.terminal.fpa_status = "UNCONSTRAINED_IN_PAPER";
    cfg.terminal.heading_status = "UNCONSTRAINED_IN_PAPER";

    % Table 3.
    cfg.table3.max_dynamic_pressure_Pa = 20e3;
    cfg.table3.max_load_factor = 2.5;
    cfg.table3.max_heating_rate_value = 160e3;
    cfg.table3.max_heating_rate_unit = "W";
    cfg.table3.max_bank_angle_deg = 60;
    cfg.table3.max_bank_rate_deg_s = 10;
    cfg.table3.antenna_beam_half_angle_deg = 45;
    cfg.table3.blackout_lower_altitude_m = 60e3;
    cfg.table3.tdrs_longitude_deg = 77;
    cfg.table3.tdrs_latitude_deg = 0;
    cfg.table3.tdrs_geocentric_radius_m = 42164e3;

    cfg.communication.relay_mode = "PAPER_TDRS_STATIC_EARTH_FIXED";
    cfg.communication.antenna_mount = "PAPER_TOP";
    cfg.communication.antenna_boresight_body = [0; 1; 0];
    cfg.communication.tracking_scope = "PAPER_BZC";
    cfg.communication.blackout_upper_altitude_m = cfg.table2.initial_altitude_m;
    cfg.communication.blackout_lower_altitude_m = cfg.table3.blackout_lower_altitude_m;
    cfg.communication.geometry_scope = "RAAP_ONLY_NO_PUBLISHED_LINK_BUDGET";
    cfg.communication.earth_fixed_to_eci_angle_at_study_epoch_deg = 0;
    cfg.communication.eci_epoch_alignment_status = "CONVENTION_REQUIRED";

    % Eq. (32) and the incompatible Table 3/Fig. 8 limits.
    cfg.heating.model = "KQ_SQRT_RHO_V_POWER_3P15";
    cfg.heating.velocity_exponent = 3.15;
    cfg.heating.kq = NaN;
    cfg.heating.table3_limit_value = 160e3;
    cfg.heating.table3_limit_unit = "W";
    cfg.heating.figure8_boundary_approx_W_m2 = 1.6e6;
    cfg.heating.source_values_consistent = false;
    cfg.heating.comparison_status = "UNAVAILABLE";

    cfg.ocp.augmented_state = ["r", "longitude", "latitude", ...
                               "speed", "fpa", "heading", "bank"];
    cfg.ocp.control = "BANK_RATE";
    cfg.ocp.objective = "MINIMIZE_INTEGRAL_ABS_BANK_RATE";
    cfg.ocp.formulation = "TWO_PHASE_P2";
    cfg.ocp.solver_named_by_paper = "GPOPS-II";
    cfg.ocp.optimized_bank_history_status = "UNAVAILABLE";
    cfg.ocp.mesh_and_nlp_settings_status = "UNAVAILABLE";

    % These values are not attributed to the paper. They only make the
    % optional shared Reentry_Propagator adapter executable.
    cfg.surrogate.shape_name = "COMPROMISE";
    cfg.surrogate.mass_kg = 2000;
    cfg.surrogate.reference_area_m2 = pi * 1.944 * 1.296 / 4;
    cfg.surrogate.density_model = "PROJECT_ISA76";
    cfg.surrogate.earth_rotation_rad_s = 7.2921159e-5;
    cfg.surrogate.bank_angle_deg = 0;
    cfg.surrogate.heating_model = "SUTTON_GRAVES";
    cfg.surrogate.sutton_graves_k = 1.83e-4;
    cfg.surrogate.dt_s = 2;
    cfg.surrogate.max_time_s = 2000;
end
