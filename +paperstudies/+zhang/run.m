function result = run(options)
%RUN Execute Zhang equation checks and optional open-loop forward cases.
%   RESULT = PAPERSTUDIES.ZHANG.RUN() is fast and deterministic; it does
%   not propagate or optimize. Set OPTIONS.FORWARD=true to run the selected
%   cases through Reentry_Propagator using explicitly labeled surrogate
%   mass, area, density, heating, and constant bank assumptions.

    if nargin < 1 || isempty(options)
        options = struct();
    elseif (islogical(options) || isnumeric(options)) && isscalar(options)
        options = struct('forward', logical(options));
    end
    if ~isstruct(options)
        error('paperstudies:zhang:run:InvalidOptions', ...
              'options must be a struct or scalar logical forward flag.');
    end

    cfg = paperstudies.zhang.config();
    case_ids = option_local(options, 'case_ids', 1:numel(cfg.cases));
    forward_value = option_local(options, 'forward', false);
    if ~(islogical(forward_value) || isnumeric(forward_value)) || ...
            ~isscalar(forward_value) || ~isfinite(double(forward_value))
        error('paperstudies:zhang:run:InvalidForwardFlag', ...
              'options.forward must be a finite scalar logical flag.');
    end
    forward = logical(forward_value);
    mass_kg = option_local(options, 'mass_kg', cfg.surrogate.mass_kg);
    area_m2 = option_local(options, 'reference_area_m2', ...
        cfg.surrogate.reference_area_m2);
    bank_deg = option_local(options, 'bank_angle_deg', ...
        cfg.surrogate.bank_angle_deg);
    dt_s = option_local(options, 'dt_s', cfg.surrogate.dt_s);
    max_time_s = option_local(options, 'max_time_s', ...
        cfg.surrogate.max_time_s);
    shape_name = string(option_local(options, 'shape_name', ...
        cfg.surrogate.shape_name));

    validateattributes(mass_kg, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(area_m2, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(bank_deg, {'numeric'}, {'scalar','real','finite'});
    validateattributes(dt_s, {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(max_time_s, {'numeric'}, {'scalar','real','finite','positive'});
    if ~isscalar(shape_name) || strlength(shape_name) == 0
        error('paperstudies:zhang:run:InvalidShape', ...
              'options.shape_name must be a nonempty scalar string.');
    end

    states = paperstudies.zhang.entry_state(case_ids, cfg, mass_kg);
    case_results = repmat(struct(), 1, numel(states));
    for i = 1:numel(states)
        [alpha, alpha_info] = paperstudies.zhang.aoa(states(i).speed_m_s, cfg);
        [cl, cd, aero_info] = paperstudies.zhang.aero(alpha, cfg);
        [beta, raap_info] = paperstudies.zhang.raap( ...
            states(i), alpha, bank_deg, cfg);
        case_results(i).case_id = states(i).case_id;
        case_results(i).state = states(i);
        case_results(i).initial_aoa_deg = alpha;
        case_results(i).initial_cl = cl;
        case_results(i).initial_cd = cd;
        case_results(i).initial_raap_deg = beta;
        case_results(i).initial_raap_within_cone = ...
            raap_info.within_paper_cone;
        case_results(i).equation_metadata.aoa = alpha_info;
        case_results(i).equation_metadata.aero = aero_info;
        case_results(i).equation_metadata.raap = raap_info;
    end

    [cl_check, cd_check] = paperstudies.zhang.aero( ...
        cfg.aero.regression_aoa_deg, cfg);
    aero_error = max(abs([cl_check-cfg.aero.regression_cl, ...
                          cd_check-cfg.aero.regression_cd]));
    result.schema_version = cfg.schema_version;
    result.study_id = cfg.study_id;
    result.mode = "EQUATION_AND_CONFIGURATION_AUDIT";
    result.config = cfg;
    result.reproduction_status = paperstudies.zhang.status(cfg);
    result.cases = case_results;
    result.regression.aero_max_abs_error = aero_error;
    result.regression.aero_passed = aero_error <= 1e-12;
    result.regression.table2_case_count_passed = numel(cfg.cases) == 3;
    result.regression.paper_setup_passed = ...
        cfg.communication.relay_mode == "PAPER_TDRS_STATIC_EARTH_FIXED" && ...
        cfg.communication.antenna_mount == "PAPER_TOP" && ...
        cfg.communication.tracking_scope == "PAPER_BZC";
    result.published_reference.terminal_bounds.longitude_deg = ...
        cfg.terminal.longitude_deg+[-1,1]*cfg.terminal.longitude_tolerance_deg;
    result.published_reference.terminal_bounds.latitude_deg = ...
        cfg.terminal.latitude_deg+[-1,1]*cfg.terminal.latitude_tolerance_deg;
    result.published_reference.terminal_bounds.altitude_m = ...
        cfg.terminal.altitude_m+[-1,1]*cfg.terminal.altitude_tolerance_m;
    result.published_reference.terminal_bounds.speed_m_s = ...
        cfg.terminal.speed_m_s+[-1,1]*cfg.terminal.speed_tolerance_m_s;
    result.published_reference.path_limits.dynamic_pressure_Pa = ...
        cfg.table3.max_dynamic_pressure_Pa;
    result.published_reference.path_limits.load_factor = ...
        cfg.table3.max_load_factor;
    result.published_reference.path_limits.bank_angle_deg = ...
        cfg.table3.max_bank_angle_deg;
    result.published_reference.path_limits.bank_rate_deg_s = ...
        cfg.table3.max_bank_rate_deg_s;
    result.published_reference.path_limits.bzc_raap_deg = ...
        cfg.table3.antenna_beam_half_angle_deg;
    result.published_reference.path_limits.heating_status = ...
        "UNAVAILABLE_CONFLICTING_SOURCE_DATA";
    result.comparison.performed = false;
    result.comparison.status = "UNAVAILABLE_NO_FORWARD_TRAJECTORY";
    result.forward.executed = false;
    result.forward.classification = "NOT_REQUESTED";

    if ~forward
        return;
    end

    [sys, provenance] = paperstudies.zhang.configure_system(options, cfg);
    result.mode = "SURROGATE_OPEN_LOOP_FORWARD";
    result.forward.executed = true;
    result.forward.classification = provenance.classification;
    result.forward.provenance = provenance;
    result.forward.assumptions.mass_kg = mass_kg;
    result.forward.assumptions.reference_area_m2 = area_m2;
    result.forward.assumptions.density_model = cfg.surrogate.density_model;
    result.forward.assumptions.shape_name = shape_name;
    result.forward.assumptions.earth_rotation_rad_s = ...
        cfg.surrogate.earth_rotation_rad_s;
    result.forward.assumptions.heating_model = ...
        cfg.surrogate.heating_model;
    result.forward.assumptions.bank_angle_deg = bank_deg;
    result.forward.assumptions.bank_history = "CONSTANT_NOT_OPTIMIZED";
    result.forward.cases = repmat(struct(), 1, numel(states));
    result.comparison.performed = true;
    result.comparison.status = "AVAILABLE_PER_FORWARD_CASE";

    for i = 1:numel(states)
        params.vehicle_mode = "SPACEPLANE";
        params.shape_name = shape_name;
        params.gravity_model = "CENTRAL_SPHERICAL";
        params.dt = dt_s;
        params.max_time = max_time_s;
        params.terminal_altitude = cfg.terminal.altitude_m;
        params.bank_angle_deg = bank_deg;
        params.lift_enabled = true;
        params.reference_area_m2 = area_m2;
        [X_final, ~, hist, summary] = Reentry_Propagator( ...
            sys, states(i).X_eci14, [], 0, params);
        result.forward.cases(i).case_id = states(i).case_id;
        result.forward.cases(i).X_final = X_final;
        result.forward.cases(i).hist = hist;
        result.forward.cases(i).summary = summary;
        result.forward.cases(i).comparison = ...
            paperstudies.zhang.compare_forward(hist, summary, cfg);
    end
end

function value = option_local(options, name, default_value)
    if isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = default_value;
    end
end
