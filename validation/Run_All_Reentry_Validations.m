function results = Run_All_Reentry_Validations()
%RUN_ALL_REENTRY_VALIDATIONS Run the integrated re-entry regression suite.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    validation_dir = fileparts(mfilename('fullpath'));
    addpath(project_root);
    addpath(validation_dir);

    fprintf('\n=== Re-entry validation suite ===\n');
    results.code_analyzer = Check_Reentry_Code();
    results.core_equivalence = Validate_Reentry_Core_Equivalence();
    results.integrated_propagator = Validate_Reentry_Propagator();
    results.paper_reproduction = Validate_Paper_Reproduction();
    results.passed = true;
    fprintf('=== All re-entry validations: PASS ===\n\n');
end
