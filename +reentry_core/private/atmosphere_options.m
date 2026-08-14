function atmosphere_opts = atmosphere_options(sys)
    atmosphere_opts = struct();
    if isfield(sys, 'environment') && isfield(sys.environment, 'atmospheric_drag')
        atmosphere_opts = sys.environment.atmospheric_drag;
    end
    atmosphere_opts.enabled = true;
    atmosphere_opts.model = "ISA76";
end
