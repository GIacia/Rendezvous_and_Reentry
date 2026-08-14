function R_ref_printed = eq17_reference_range( ...
    ld_ratio, speed, gravity, earth_radius)
%EQ17_REFERENCE_RANGE Exact printed algebra of Saito Eq. (17).
%   The result as printed is dimensionless although Table 5 labels Rref in
%   km. No unprinted Earth-radius scale is applied here; callers must make
%   any redimensionalization explicit and classify it CONVENTION_REQUIRED.

    argument = 1 - speed.^2 ./ (gravity .* earth_radius);
    if any(~isfinite(argument(:))) || any(argument(:) <= 0)
        error('paperstudies:saito:eq17:InvalidLogArgument', ...
              'Eq. (17) requires 1 - V^2/(g*r_E) to be finite and positive.');
    end
    R_ref_printed = -0.5 .* ld_ratio .* log(argument);
end
