function F2 = eq21_gain_F2(K2, R_ref, ld_ratio)
%EQ21_GAIN_F2 Exact public algebra of Saito Eq. (21).
%   K2 is intentionally caller-supplied; the paper does not publish it.
    if any(~isfinite(ld_ratio(:))) || any(ld_ratio(:) == 0)
        error('paperstudies:saito:eq21:InvalidLD', ...
              'ld_ratio must contain finite nonzero values.');
    end
    F2 = K2 .* R_ref.^2 ./ ld_ratio;
end
