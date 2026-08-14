function [beta_deg, detail] = raap(state, alpha_deg, bank_deg, cfg)
%RAAP Implement Zhang Eqs. (1)-(14) with the paper +Y_B boresight.

    if nargin < 4 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    required = {'altitude_m','longitude_deg','latitude_deg', ...
                'fpa_deg','heading_deg'};
    for i = 1:numel(required)
        if ~isstruct(state) || ~isfield(state, required{i}) || ...
                ~isscalar(state.(required{i})) || ...
                ~isfinite(state.(required{i}))
            error('paperstudies:zhang:raap:InvalidState', ...
                  'State must contain finite scalar spherical-entry fields.');
        end
    end
    if ~isscalar(alpha_deg) || ~isfinite(alpha_deg) || ...
            ~isscalar(bank_deg) || ~isfinite(bank_deg)
        error('paperstudies:zhang:raap:InvalidAttitude', ...
              'AoA and bank must be finite scalar angles in degrees.');
    end

    theta = deg2rad(state.longitude_deg);
    phi = deg2rad(state.latitude_deg);
    gamma = deg2rad(state.fpa_deg);
    psi = deg2rad(state.heading_deg);
    alpha = deg2rad(alpha_deg);
    sigma = deg2rad(bank_deg);
    theta_s = deg2rad(cfg.table3.tdrs_longitude_deg);
    phi_s = deg2rad(cfg.table3.tdrs_latitude_deg);

    r = cfg.earth.radius_m + state.altitude_m;
    r_s = cfg.table3.tdrs_geocentric_radius_m;
    P_o = [r*sin(phi); r*cos(phi)*cos(theta); r*cos(phi)*sin(theta)];
    S_o = [r_s*sin(phi_s); r_s*cos(phi_s)*cos(theta_s); ...
           r_s*cos(phi_s)*sin(theta_s)];
    W_o = S_o-P_o;

    G_oe = [cos(-phi), sin(-phi), 0; ...
           -sin(-phi), cos(-phi), 0; 0, 0, 1] * ...
           [1, 0, 0; 0, cos(theta), sin(theta); ...
            0, -sin(theta), cos(theta)];
    G_et = [cos(gamma), sin(gamma), 0; ...
           -sin(gamma), cos(gamma), 0; 0, 0, 1] * ...
           [cos(-psi), 0, -sin(-psi); 0, 1, 0; ...
            sin(-psi), 0, cos(-psi)];
    G_tb = [cos(alpha), sin(alpha), 0; ...
           -sin(alpha), cos(alpha), 0; 0, 0, 1] * ...
           [1, 0, 0; 0, cos(sigma), sin(sigma); ...
            0, -sin(sigma), cos(sigma)];
    W_b = G_tb*G_et*G_oe*W_o;
    A_b = cfg.communication.antenna_boresight_body;
    cosine_beta = dot(A_b,W_b)/(norm(A_b)*norm(W_b));
    cosine_beta = min(1,max(-1,cosine_beta));
    beta_deg = acosd(cosine_beta);

    detail.classification = "EXACT";
    detail.W_geocentric = W_o;
    detail.W_body = W_b;
    detail.G_oe = G_oe;
    detail.G_et = G_et;
    detail.G_tb = G_tb;
    detail.boresight_body = A_b;
    detail.within_paper_cone = ...
        beta_deg <= cfg.table3.antenna_beam_half_angle_deg;
end
