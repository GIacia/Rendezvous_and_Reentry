function forces_wind = eq12_force_transform( ...
    forces_body, alpha_rad, beta_rad, sigma_rad)
%EQ12_FORCE_TRANSFORM Saito Eq. (12), body to wind-axis force vector.
%   FORCES_BODY is [D; F; L] in the notation used by the paper.
    if numel(forces_body) ~= 3
        error('paperstudies:saito:eq12:InvalidForce', ...
              'forces_body must be [D; F; L].');
    end
    ca = cos(alpha_rad); sa = sin(alpha_rad);
    cb = cos(beta_rad);  sb = sin(beta_rad);
    cs = cos(sigma_rad); ss = sin(sigma_rad);
    transform = [ ...
        ca*cb,              ca*sb,              -sa; ...
        ss*sa*cb-cs*sb,     ss*sa*sb+cs*cb,      ss*ca; ...
        cs*sa*cb+ss*sb,     cs*sa*sb-ss*cb,      cs*ca];
    forces_wind = transform * forces_body(:);
end
