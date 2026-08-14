function report = Run_Paper_Reproduction_Suite(options)
%RUN_PAPER_REPRODUCTION_SUITE Audit Zhang and Saito paper-study models.
%   REPORT = RUN_PAPER_REPRODUCTION_SUITE() performs a fast deterministic
%   audit of published tables, equations, geometry, entry-state construction,
%   uncertainty grids, and provenance. It does not optimize or propagate.
%
%   REPORT = RUN_PAPER_REPRODUCTION_SUITE(STRUCT('forward',true)) also runs
%   the optional open-loop surrogate adapters through Reentry_Propagator.
%   Forward results never constitute full numerical paper reproduction.
%
%   Paper-specific options may be supplied under OPTIONS.ZHANG and
%   OPTIONS.SAITO. The top-level FORWARD value is used only when the nested
%   structure does not explicitly contain its own forward field. Selected
%   equation/coordinate self-tests run by default and can be disabled with
%   OPTIONS.RUN_SELFTESTS=false for data inspection only.

    if nargin < 1 || isempty(options)
        options = struct();
    end
    if ~isstruct(options) || ~isscalar(options)
        error('paperstudies:suite:InvalidOptions', ...
              'options must be a scalar struct.');
    end

    defaults = struct( ...
        'forward', false, ...
        'run_selftests', true, ...
        'print_summary', true, ...
        'zhang', struct(), ...
        'saito', struct());
    options = merge_options_local(defaults, options);
    validate_flag_local(options.forward, 'forward');
    validate_flag_local(options.run_selftests, 'run_selftests');
    validate_flag_local(options.print_summary, 'print_summary');
    if ~isstruct(options.zhang) || ~isscalar(options.zhang) || ...
            ~isstruct(options.saito) || ~isscalar(options.saito)
        error('paperstudies:suite:InvalidNestedOptions', ...
              'options.zhang and options.saito must be scalar structs.');
    end

    zhang_options = apply_forward_default_local(options.zhang, options.forward);
    saito_options = apply_forward_default_local(options.saito, options.forward);
    zhang = paperstudies.zhang.run(zhang_options);
    saito = paperstudies.saito.run(saito_options);

    if logical(options.run_selftests)
        zhang_selftest = paperstudies.zhang.selftest(false);
        saito_selftest = paperstudies.saito.selftest(false);
        zhang_selftest_passed = zhang_selftest.default_audit_passed;
        saito_selftest_passed = saito_selftest.passed;
    else
        zhang_selftest = struct('executed', false);
        saito_selftest = struct('executed', false);
        zhang_selftest_passed = NaN;
        saito_selftest_passed = NaN;
    end

    zhang_audit_passed = zhang.regression.aero_passed && ...
        zhang.regression.table2_case_count_passed && ...
        zhang.regression.paper_setup_passed && ...
        (~logical(options.run_selftests) || zhang_selftest_passed);
    saito_audit_passed = numel(saito.entry_cases) == 4 && ...
        saito.table6_deorbit_grid.count == 27 && ...
        saito.explicit_grid.count == 15 && ...
        saito.rpc_grid.count == 225 && ...
        saito.range_to_go_inconsistency_detected && ...
        (~logical(options.run_selftests) || saito_selftest_passed);

    report.schema_version = "1.0";
    report.classification = ...
        "PUBLISHED_STRUCTURE_AND_SELECTED_REGRESSION_AUDIT";
    report.zhang = zhang;
    report.saito = saito;
    report.audit.zhang_passed = zhang_audit_passed;
    report.audit.saito_passed = saito_audit_passed;
    report.audit.all_passed = zhang_audit_passed && saito_audit_passed;
    report.audit.scope = ...
        "SELECTED_ANCHORS_EQUATIONS_COORDINATES_AND_CASE_ENUMERATION";
    report.audit.full_source_transcription_proven = false;
    report.audit.selftests_requested = logical(options.run_selftests);
    report.audit.selftests.zhang = zhang_selftest;
    report.audit.selftests.saito = saito_selftest;
    report.evidence.shared_kernel_package_available = ...
        ~isempty(which('reentry_core.evaluate_state'));
    if zhang.forward.executed && saito.forward.executed
        report.evidence.same_public_propagator_used_by_surrogates = true;
    else
        report.evidence.same_public_propagator_used_by_surrogates = NaN;
    end
    report.evidence.core_equivalence_regression = ...
        "VALIDATE_REENTRY_CORE_EQUIVALENCE";
    report.evidence.forward_surrogates_executed = ...
        zhang.forward.executed && saito.forward.executed;
    report.evidence.forward_adapter_smoke_regression = ...
        "VALIDATE_PAPER_REPRODUCTION";
    report.evidence.full_paper_trajectory_reproduction_established = false;
    report.evidence.external_validation_established = false;
    report.valid_claim = ...
        "SHARED_KERNEL_AND_SELECTED_PUBLISHED_STRUCTURES_ARE_REGRESSION_TESTABLE";
    report.invalid_claim = ...
        "FULL_NUMERICAL_TRAJECTORY_OR_LANDING_DISPERSION_REPRODUCTION";
    report.research_use = ...
        "PRELIMINARY_PARAMETERS_AND_SENSITIVITY_WITH_EXPLICIT_PROVENANCE";
    report.blockers.zhang = zhang.reproduction_status.blockers;
    report.blockers.saito = string({saito.status.blockers.id}).';

    if logical(options.print_summary)
        print_summary_local(report);
    end
end

function options = merge_options_local(defaults, supplied)
    options = defaults;
    supplied_names = fieldnames(supplied);
    allowed_names = fieldnames(defaults);
    for i = 1:numel(supplied_names)
        name = supplied_names{i};
        if ~any(strcmp(name, allowed_names))
            error('paperstudies:suite:UnknownOption', ...
                  'Unknown option: %s', name);
        end
        options.(name) = supplied.(name);
    end
end

function nested = apply_forward_default_local(nested, forward)
    if ~isfield(nested, 'forward') || isempty(nested.forward)
        nested.forward = logical(forward);
    end
end

function validate_flag_local(value, name)
    if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
            ~isfinite(double(value)) || ~any(double(value) == [0, 1])
        error('paperstudies:suite:InvalidFlag', ...
              'options.%s must be a scalar logical or numeric 0/1.', name);
    end
end

function print_summary_local(report)
    fprintf('\nPaper reproduction selected-regression audit\n');
    fprintf('  Zhang selected condition checks : %s\n', ...
        pass_text_local(report.audit.zhang_passed));
    fprintf('  Saito selected condition checks : %s\n', ...
        pass_text_local(report.audit.saito_passed));
    fprintf('  Shared kernel package available : %s\n', ...
        pass_text_local(report.evidence.shared_kernel_package_available));
    fprintf('  Full numerical paper reproduction: NOT ESTABLISHED\n');
    fprintf(['  Run validation/Validate_Reentry_Core_Equivalence.m for the ' ...
             'kernel extraction regression.\n\n']);
end

function value = pass_text_local(flag)
    if flag
        value = 'PASS';
    else
        value = 'FAIL';
    end
end
