function [alpha_deg, info] = aoa(speed_m_s, cfg)
%AOA Zhang Fig. 2 prescribed angle-of-attack schedule.

    if nargin < 2 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    if ~isnumeric(speed_m_s) || ~isreal(speed_m_s) || ...
            any(~isfinite(speed_m_s(:))) || any(speed_m_s(:) < 0)
        error('paperstudies:zhang:aoa:InvalidSpeed', ...
              'Speed must contain finite nonnegative real values.');
    end

    v1 = cfg.aero.aoa_break_speed_m_s(1);
    v2 = cfg.aero.aoa_break_speed_m_s(2);
    a1 = cfg.aero.aoa_plateau_deg(1);
    a2 = cfg.aero.aoa_plateau_deg(2);
    alpha_deg = a1 * ones(size(speed_m_s));
    transition = speed_m_s > v1 & speed_m_s < v2;
    alpha_deg(transition) = a1 + (a2-a1) .* ...
        (speed_m_s(transition)-v1) ./ (v2-v1);
    alpha_deg(speed_m_s >= v2) = a2;

    info.classification = "EXACT";
    info.figure_displayed_domain_m_s = ...
        cfg.aero.figure_displayed_speed_domain_m_s;
    info.extrapolated_outside_figure = ...
        speed_m_s < info.figure_displayed_domain_m_s(1) | ...
        speed_m_s > info.figure_displayed_domain_m_s(2);
end
