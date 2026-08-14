function F3 = eq22_gain_F3(K3, R_ref, ld_ratio)
%EQ22_GAIN_F3 Exact public algebra of Saito Eq. (22).
%   K3 is intentionally caller-supplied; the paper does not publish it.
    if any(~isfinite(ld_ratio(:))) || any(ld_ratio(:) == 0)
        error('paperstudies:saito:eq22:InvalidLD', ...
              'ld_ratio must contain finite nonzero values.');
    end
    F3 = K3 .* R_ref.^2 ./ ld_ratio;
end
