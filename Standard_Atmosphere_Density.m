function [rho, T, P, a_sound] = Standard_Atmosphere_Density(altitude_m, atmosphere_opts)
    % Standard atmosphere helper.
    %
    % Returns density [kg/m^3], temperature [K], pressure [Pa], and speed of
    % sound [m/s]. MATLAB atmosisa is used when available below 86 km; above
    % that, density falls back to a log-interpolated 1976 standard atmosphere
    % table suitable for entry-interface studies.

    if nargin < 2 || isempty(atmosphere_opts)
        atmosphere_opts = struct();
    end

    altitude_m = max(0, altitude_m);

    use_atmosisa = get_atmos_bool_local(atmosphere_opts, 'use_matlab_atmosisa', true);
    if use_atmosisa && altitude_m <= 84852 && exist('atmosisa', 'file') == 2
        try
            [T, a_sound, P, rho] = atmosisa(altitude_m);
            if all(isfinite([rho, T, P, a_sound])) && rho >= 0
                return;
            end
        catch
            % Fall through to the built-in model.
        end
    end

    if altitude_m <= 84852
        [rho, T, P] = isa_lower_atmosphere_local(altitude_m);
    else
        rho = isa76_upper_density_table_local(altitude_m);
        T = 186.946;
        P = rho * 287.05287 * T;
    end

    a_sound = sqrt(1.4 * 287.05287 * T);
end

function [rho, T, P] = isa_lower_atmosphere_local(h)
    hb = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];
    Tb = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946];
    Pb = [101325, 22632.06, 5474.889, 868.019, 110.906, 66.939, 3.9564, 0.3734];
    Lb = [-0.0065, 0, 0.0010, 0.0028, 0, -0.0028, -0.0020];

    idx = find(h >= hb, 1, 'last');
    idx = min(idx, numel(Lb));

    g0 = 9.80665;
    R = 287.05287;
    dh = h - hb(idx);
    L = Lb(idx);
    T0 = Tb(idx);
    P0 = Pb(idx);

    if abs(L) > eps
        T = T0 + L * dh;
        P = P0 * (T0 / T)^(g0 / (R * L));
    else
        T = T0;
        P = P0 * exp(-g0 * dh / (R * T));
    end

    rho = P / (R * T);
end

function rho = isa76_upper_density_table_local(altitude_m)
    alt_m = [ ...
        80000, 90000, 100000, 110000, 120000, 130000, 140000, ...
        150000, 160000, 170000, 180000, 190000, 200000, 250000, ...
        300000, 350000, 400000, 450000, 500000, 600000, 700000, ...
        800000, 900000, 1000000 ];

    rho_table = [ ...
        1.905e-5, 3.396e-6, 5.604e-7, 9.708e-8, 2.222e-8, ...
        8.152e-9, 3.831e-9, 2.076e-9, 1.234e-9, 7.824e-10, ...
        5.194e-10, 3.581e-10, 2.541e-10, 6.073e-11, 1.916e-11, ...
        7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 1.137e-13, ...
        3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15 ];

    if altitude_m <= alt_m(1)
        rho = rho_table(1);
    elseif altitude_m > alt_m(end)
        rho = 0;
    else
        rho = exp(interp1(alt_m, log(rho_table), altitude_m, 'linear'));
    end
end

function value = get_atmos_bool_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end

    if islogical(value)
        return;
    end
    if isnumeric(value)
        value = value ~= 0;
        return;
    end

    txt = lower(strtrim(string(value)));
    value = txt == "true" || txt == "1" || txt == "yes" || txt == "on";
end
