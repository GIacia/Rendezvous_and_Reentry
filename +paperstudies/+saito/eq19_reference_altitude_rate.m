function hdot_ref = eq19_reference_altitude_rate(ld_ratio, V0, speed)
%EQ19_REFERENCE_ALTITUDE_RATE Exact public algebra of Saito Eq. (19).
%   V0 is exposed because its numerical value/convention is not supplied
%   sufficiently for reproduction of Table 5.
    hdot_ref = ld_ratio .* (V0 - speed);
end
