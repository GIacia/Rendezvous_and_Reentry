function results = selftest(print_summary)
%SELFTEST Focused deterministic checks for the public Saito-study API.

    if nargin < 1 || isempty(print_summary)
        print_summary = true;
    end
    if ~(islogical(print_summary) || isnumeric(print_summary)) || ...
            ~isscalar(print_summary) || ~isfinite(double(print_summary)) || ...
            ~any(double(print_summary) == [0, 1])
        error('paperstudies:saito:selftest:InvalidFlag', ...
              'print_summary must be a scalar logical or numeric 0/1.');
    end

    cfg = paperstudies.saito.config();
    assert(cfg.tables.table1.total_mass_kg == 330);
    assert(isequal(size(cfg.tables.table4.values), [4, 7]));
    assert(isequal(size(cfg.tables.table5.values), [13, 7]));
    assert(cfg.tables.table6.total_number_of_cases == 27);
    assert(isequal(cfg.tables.table6.entry_range_to_go_bounds_km, [2369, 7031]));
    assert(all(isnan(cfg.tables.table6.individual_entry_range_to_go_km)));
    assert(~cfg.tables.table6.individual_entry_states_available);
    assert(isequal(size(cfg.tables.table7.values), [15, 4]));
    assert(isequal(size(cfg.tables.table10.values), [5, 7]));
    assert(isequal(size(cfg.tables.table11.values), [5, 7]));
    assert(cfg.tables.table5.values(11,2) == 1260.82);
    assert(cfg.tables.table10.values(5,3) == 104.26);
    assert(cfg.tables.table11.values(5,2) == -30.93);

    table6_grid = paperstudies.saito.table6_deorbit_grid(cfg);
    explicit_grid = paperstudies.saito.explicit_uncertainty_grid(cfg);
    rpc_grid = paperstudies.saito.rpc_uncertainty_grid(cfg);
    assert(table6_grid.count == 27);
    assert(table6_grid.ignition_epoch_utc(1) == "2026-07-05T22:34:00Z");
    assert(table6_grid.ignition_epoch_utc(end) == "2026-07-05T22:47:00Z");
    assert(all(isnan(table6_grid.values(:,3))));
    assert(~table6_grid.individual_range_to_go_available);
    assert(explicit_grid.count == 15);
    assert(rpc_grid.count == 225);
    assert(sum(rpc_grid.guidance_method == "RPC-1") == 75);
    assert(sum(rpc_grid.guidance_method == "RPC-2") == 75);
    assert(sum(rpc_grid.guidance_method == "RPC-3") == 75);

    entries = repmat(paperstudies.saito.entry_state(1, cfg), 1, 4);
    for case_id = 1:4
        entries(case_id) = paperstudies.saito.entry_state(case_id, cfg);
        error_values = struct2array(entries(case_id).roundtrip_error);
        assert(max(abs(error_values)) < 1e-8);
        assert(entries(case_id).range_diagnostic.inconsistency_detected);
        inferred_radius = entries(case_id).range_diagnostic.effective_radius_to_match_published_m;
        assert(inferred_radius > 6.462e6 && inferred_radius < 6.465e6);
    end

    assert(abs(paperstudies.saito.eq13_ld_command(0.2, 0.01, 10, 5) - 0.25) < 1e-14);
    [range, angle] = paperstudies.saito.eq14_range_to_go(0, 0, 0, pi/2, 10);
    assert(abs(angle - pi/2) < 1e-14);
    assert(abs(range - 5*pi) < 1e-13);
    [sigma15, cosine15] = paperstudies.saito.eq15_bank_command(0.125, 0.25);
    assert(abs(cosine15 - 0.5) < 1e-14);
    assert(abs(rad2deg(sigma15) - 60) < 1e-12);
    guidance_range = paperstudies.saito.eq16_guidance_range(100, 2, 3, 2, 0.1, -10, -20);
    assert(abs(guidance_range - 103) < 1e-14);
    equation_argument = 0.5;
    equation_speed = sqrt((1 - equation_argument) * ...
        cfg.conventions.gravity_m_s2 * cfg.conventions.earth_radius_m);
    R_ref_eq17 = paperstudies.saito.eq17_reference_range( ...
        0.25, equation_speed, cfg.conventions.gravity_m_s2, ...
        cfg.conventions.earth_radius_m);
    assert(abs(R_ref_eq17 - (-0.125 * log(equation_argument))) < 1e-14);
    D_ref_eq18 = paperstudies.saito.eq18_reference_drag( ...
        0.25, equation_speed, cfg.conventions.gravity_m_s2, ...
        cfg.conventions.earth_radius_m);
    assert(abs(D_ref_eq18 - 2) < 1e-14);
    assert(abs(paperstudies.saito.eq19_reference_altitude_rate( ...
        0.25, 8000, 7000) - 250) < 1e-14);
    assert(abs(paperstudies.saito.eq20_gain_F1( ...
        2, equation_speed, cfg.conventions.gravity_m_s2, ...
        cfg.conventions.earth_radius_m) - 4/(-2*log(0.5))) < 1e-13);
    assert(abs(paperstudies.saito.eq21_gain_F2(3, 4, 0.25) - 192) < 1e-14);
    assert(abs(paperstudies.saito.eq22_gain_F3(5, 4, 0.25) - 320) < 1e-14);
    [sigma23, cosine23] = paperstudies.saito.eq23_rpc_bank_command( ...
        deg2rad(60), deg2rad(10), 95, 100, 90);
    expected_cosine23 = cosd(60) + (cosd(70) - cosd(60)) / (-10) * (-5);
    assert(abs(cosine23 - expected_cosine23) < 1e-14);
    assert(abs(cos(sigma23) - expected_cosine23) < 1e-14);
    assert(abs(paperstudies.saito.eq24_fading_filter(1, 0.9, 1.1) - 1.01) < 1e-14);
    assert(abs(paperstudies.saito.eq25_terminal_bank_secant(0, 1, 10, -5) - 2/3) < 1e-14);
    [CD, CL] = paperstudies.saito.eq26_aero_coefficients(1.2, 0.3, 0);
    assert(abs(CD - 1.2) < 1e-14 && abs(CL - 0.3) < 1e-14);

    radius = cfg.conventions.earth_radius_m + 100e3;
    gravity = cfg.conventions.earth_mu_m3_s2 / radius^2;
    circular_speed = sqrt(gravity * radius);
    derivative = paperstudies.saito.eq10_spherical_dynamics( ...
        [radius;0;0;circular_speed;0;pi/2], [0;0;0], ...
        150, gravity, 0);
    assert(abs(derivative(4)) < 1e-12);
    assert(abs(derivative(5)) < 1e-12);

    default_result = paperstudies.saito.run();
    assert(~default_result.forward.executed);
    assert(default_result.table6_deorbit_grid.count == 27);
    assert(default_result.explicit_grid.count == 15);
    assert(default_result.rpc_grid.count == 225);
    assert(default_result.status.overall == "PARTIAL_REPRODUCTION");

    results.passed = true;
    results.entry_case_count = numel(entries);
    results.table6_case_count = table6_grid.count;
    results.explicit_case_count = explicit_grid.count;
    results.rpc_case_count = rpc_grid.count;
    results.maximum_roundtrip_error = max(abs([ ...
        struct2array(entries(1).roundtrip_error), ...
        struct2array(entries(2).roundtrip_error), ...
        struct2array(entries(3).roundtrip_error), ...
        struct2array(entries(4).roundtrip_error)]));
    results.range_discrepancy_m = arrayfun( ...
        @(entry) entry.range_diagnostic.difference_m, entries);

    if logical(print_summary)
        fprintf('Saito paper-study self-test: PASS\n');
        fprintf('  Table-4 entry cases : %d\n', results.entry_case_count);
        fprintf('  Table-6 epochs      : %d (individual entry states unavailable)\n', ...
            results.table6_case_count);
        fprintf('  Explicit/RPC cases  : %d / %d\n', ...
            results.explicit_case_count, results.rpc_case_count);
        fprintf('  RTG discrepancies   : %.3f to %.3f km\n', ...
            min(results.range_discrepancy_m)/1e3, ...
            max(results.range_discrepancy_m)/1e3);
    end
end
