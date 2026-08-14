function F1 = eq20_gain_F1(K1, speed, gravity, earth_radius)
%EQ20_GAIN_F1 Exact public algebra of Saito Eq. (20).
%   K1 is intentionally caller-supplied; the paper does not publish it.
    argument = 1 - speed.^2 ./ (gravity .* earth_radius);
    denominator = -K1 .* log(argument);
    if any(~isfinite(argument(:))) || any(argument(:) <= 0) || ...
            any(~isfinite(denominator(:))) || any(denominator(:) == 0)
        error('paperstudies:saito:eq20:InvalidInputs', ...
              'Eq. (20) requires a positive log argument and nonzero denominator.');
    end
    F1 = 4 ./ denominator;
end
