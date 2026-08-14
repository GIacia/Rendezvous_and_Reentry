function D_ref = eq18_reference_drag(ld_ratio, speed, gravity, earth_radius)
%EQ18_REFERENCE_DRAG Exact public algebra of Saito Eq. (18).
    if any(~isfinite(ld_ratio(:))) || any(ld_ratio(:) == 0)
        error('paperstudies:saito:eq18:InvalidLD', ...
              'ld_ratio must contain finite nonzero values.');
    end
    D_ref = (1 ./ ld_ratio) .* ...
        (1 - speed.^2 ./ (gravity .* earth_radius));
end
