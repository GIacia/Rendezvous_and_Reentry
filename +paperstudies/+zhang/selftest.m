function results = selftest(run_forward_smoke)
%SELFTEST Focused deterministic checks for the Zhang study package.

    if nargin < 1 || isempty(run_forward_smoke)
        run_forward_smoke = false;
    end
    if ~isscalar(run_forward_smoke)
        error('paperstudies:zhang:selftest:InvalidFlag', ...
              'run_forward_smoke must be scalar.');
    end

    cfg = paperstudies.zhang.config();
    audit = paperstudies.zhang.run();
    assert(~audit.forward.executed);
    assert(audit.regression.aero_passed);
    assert(audit.regression.paper_setup_passed);
    assert(numel(audit.cases) == 3);
    assert(audit.reproduction_status.components.optimized_bank_history == ...
        "UNAVAILABLE");

    for i = 1:3
        state = audit.cases(i).state;
        expected_radius = cfg.earth.radius_m+cfg.table2.initial_altitude_m;
        assert(abs(norm(state.r_earth_fixed_m)-expected_radius) < 1e-6);
        assert(abs(norm(state.v_atmosphere_relative_m_s)- ...
            cfg.cases(i).speed_m_s) < 1e-9);
        up = state.r_earth_fixed_m/norm(state.r_earth_fixed_m);
        recovered_fpa = asind(dot(state.v_atmosphere_relative_m_s,up)/ ...
            norm(state.v_atmosphere_relative_m_s));
        assert(abs(recovered_fpa-cfg.cases(i).fpa_deg) < 1e-10);
    end

    [alpha, ~] = paperstudies.zhang.aoa([0, 2000, 3500, 5000, 8000], cfg);
    assert(max(abs(alpha-[15, 15, 27.5, 40, 40])) < 1e-12);
    [cl, cd] = paperstudies.zhang.aero(cfg.aero.regression_aoa_deg, cfg);
    assert(max(abs(cl-cfg.aero.regression_cl)) < 1e-12);
    assert(max(abs(cd-cfg.aero.regression_cd)) < 1e-12);
    dx = paperstudies.zhang.augmented_dynamics( ...
        [1; 0; 0; 1; 0; 0; 0], 0.1, 0, 0, 0);
    assert(max(abs(dx-[0; 0; 1; 0; 0; 0; 0.1])) < 1e-12);

    mixed_hist.time = [0, 1];
    mixed_hist.rv_pos = repmat([cfg.earth.radius_m+70e3; 0; 0], 1, 2);
    mixed_hist.altitude = [70e3, 69e3];
    mixed_hist.speed_rel = [7000, 6900];
    mixed_hist.bank_angle_deg = [0, 0];
    mixed_hist.blackout_zone_active = [1, 1];
    mixed_hist.raap_deg = [50, NaN];
    mixed_summary.max_dynamic_pressure_Pa = 1;
    mixed_summary.max_g_load = 1;
    mixed_summary.max_bank_rate_deg_s = 0;
    mixed_report = paperstudies.zhang.compare_forward( ...
        mixed_hist, mixed_summary, cfg);
    assert(isequal(mixed_report.constraints.bzc_satisfied, false), ...
        'A known BZC violation must dominate a separate unknown sample.');

    [sys, provenance] = paperstudies.zhang.configure_system(struct(), cfg);
    comm = sys.reentry_vehicle.spaceplane.communication;
    assert(comm.relay_mode == "PAPER_TDRS_STATIC_EARTH_FIXED");
    assert(comm.antenna_mount == "PAPER_TOP");
    assert(comm.tracking_scope == "PAPER_BZC");
    assert(isinf(sys.reentry_vehicle.spaceplane.constraints. ...
        max_heat_flux_W_m2), ...
        'Ambiguous Table-3 W value must not be used as a W/m^2 limit.');
    assert(~provenance.optimization_performed);

    results.default_audit_passed = true;
    results.forward_smoke_executed = false;
    results.forward_smoke_passed = NaN;
    if logical(run_forward_smoke)
        forward_result = paperstudies.zhang.run(struct( ...
            'forward', true, 'case_ids', 1, 'dt_s', 1, 'max_time_s', 1));
        assert(forward_result.forward.executed);
        assert(forward_result.forward.classification == ...
            "SURROGATE_OPEN_LOOP_FORWARD");
        assert(~forward_result.forward.cases(1).comparison.paper_trajectory_reproduced);
        exact_raap = forward_result.cases(1).initial_raap_deg;
        shared_raap = forward_result.forward.cases(1).hist.raap_deg(1);
        assert(abs(exact_raap-shared_raap) < 1e-10);
        results.forward_smoke_executed = true;
        results.forward_smoke_passed = true;
        results.initial_raap_adapter_error_deg = shared_raap-exact_raap;
    end
end
