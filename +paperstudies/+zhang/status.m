function s = status(cfg)
%STATUS Machine-readable Zhang reproduction capability and blockers.

    if nargin < 1 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end

    s.schema_version = cfg.schema_version;
    s.study_id = cfg.study_id;
    s.levels = cfg.levels;
    s.components.source_page_audit = "EXACT";
    s.components.table2_initial_conditions = "EXACT";
    s.components.table3_nonthermal_constraints = "EXACT";
    s.components.aoa_schedule = "EXACT";
    s.components.aerodynamic_polynomials = "EXACT";
    s.components.dimensionless_equations = "EXACT";
    s.components.augmented_bank_state_and_rate_control = "EXACT";
    s.components.tdrs_and_raap_geometry = "EXACT";
    s.components.paper_top_antenna = "EXACT";
    s.components.paper_bzc_logic = "EXACT";
    s.components.eci_earth_fixed_epoch_alignment = "CONVENTION_REQUIRED";
    s.components.terminal_box = "EXACT";
    s.components.shared_core_forward_dynamics = "SURROGATE";
    s.components.vehicle_mass = "UNAVAILABLE";
    s.components.reference_area = "UNAVAILABLE";
    s.components.density_model = "UNAVAILABLE";
    s.components.numeric_nondimensionalization_constants = "UNAVAILABLE";
    s.components.heating_coefficient_kq = "UNAVAILABLE";
    s.components.heating_limit = "UNAVAILABLE";
    s.components.initial_bank = "UNAVAILABLE";
    s.components.optimized_bank_history = "UNAVAILABLE";
    s.components.gpops_transcription_settings = "UNAVAILABLE";
    s.components.paper_trajectory_numeric_reproduction = "UNAVAILABLE";

    s.blockers = [ ...
        "vehicle_mass"; ...
        "reference_area"; ...
        "density_model"; ...
        "numeric_nondimensionalization_constants"; ...
        "heating_coefficient_kq"; ...
        "conflicting_heating_limit_and_units"; ...
        "initial_bank"; ...
        "terminal_and_switching_time_bounds"; ...
        "gpops_mesh_solver_tolerances_and_initial_guess"; ...
        "optimized_bank_and_bank_rate_history"; ...
        "raw_published_trajectory_data"];

    s.source_values_internally_consistent = false;
    s.heating_conflict = struct( ...
        'table3_value', cfg.heating.table3_limit_value, ...
        'table3_unit', cfg.heating.table3_limit_unit, ...
        'figure8_approx_W_m2', cfg.heating.figure8_boundary_approx_W_m2);
    s.full_numerical_reproduction_available = false;
    s.default_run_performs_optimization = false;
    s.forward_run_classification = "SURROGATE_OPEN_LOOP_FORWARD";
    s.valid_claim = "EQUATION_GEOMETRY_AND_CONSTRAINT_REPRODUCTION";
    s.invalid_claim = "PAPER_OPTIMAL_TRAJECTORY_REPRODUCTION";
end
