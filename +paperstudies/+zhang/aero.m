function [cl, cd, info] = aero(alpha_deg, cfg)
%AERO Zhang Eqs. (26)-(27), with alpha evaluated in degrees.

    if nargin < 2 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    if ~isnumeric(alpha_deg) || ~isreal(alpha_deg) || ...
            any(~isfinite(alpha_deg(:)))
        error('paperstudies:zhang:aero:InvalidAoA', ...
              'Angle of attack must contain finite real values in degrees.');
    end

    c = cfg.aero.cl_polynomial;
    d = cfg.aero.cd_from_cl_polynomial;
    cl = c(1) + c(2).*alpha_deg + c(3).*alpha_deg.^2;
    cd = d(1) + d(2).*cl + d(3).*cl.^2;

    info.classification = "EXACT";
    info.alpha_unit = "deg";
    info.cl_equation = "ZHANG_EQ26";
    info.cd_equation = "ZHANG_EQ27";
end
