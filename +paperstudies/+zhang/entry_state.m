function states = entry_state(case_ids, cfg, mass_kg)
%ENTRY_STATE Construct the three Table 2 states on a spherical Earth.
%   The paper state and atmosphere-relative velocity are exact. X_ECI14 is
%   an adapter for the shared project propagator and therefore uses the
%   explicitly labeled surrogate Earth rotation rate and vehicle mass.

    if nargin < 2 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    if nargin < 1 || isempty(case_ids)
        case_ids = 1:numel(cfg.cases);
    end
    if nargin < 3 || isempty(mass_kg)
        mass_kg = cfg.surrogate.mass_kg;
    end
    if ~isnumeric(case_ids) || any(case_ids ~= fix(case_ids)) || ...
            any(case_ids < 1) || any(case_ids > numel(cfg.cases))
        error('paperstudies:zhang:entryState:InvalidCase', ...
              'Case identifiers must be integers from 1 through 3.');
    end
    if ~isscalar(mass_kg) || ~isfinite(mass_kg) || mass_kg <= 0
        error('paperstudies:zhang:entryState:InvalidMass', ...
              'The surrogate adapter mass must be a finite positive scalar.');
    end

    template = struct('case_id', [], 'classification', "", ...
        'altitude_m', [], 'longitude_deg', [], 'latitude_deg', [], ...
        'speed_m_s', [], 'fpa_deg', [], 'heading_deg', [], ...
        'r_earth_fixed_m', [], 'v_atmosphere_relative_m_s', [], ...
        'v_inertial_m_s', [], 'X_eci14', [], ...
        'mass_kg', [], 'mass_classification', "");
    states = repmat(template, 1, numel(case_ids));
    omega_vec = [0; 0; cfg.surrogate.earth_rotation_rad_s];

    for k = 1:numel(case_ids)
        c = cfg.cases(case_ids(k));
        lon = deg2rad(c.longitude_deg);
        lat = deg2rad(c.latitude_deg);
        gamma = deg2rad(c.fpa_deg);
        heading = deg2rad(c.heading_deg);

        up = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
        north = [-sin(lat)*cos(lon); -sin(lat)*sin(lon); cos(lat)];
        east = [-sin(lon); cos(lon); 0];
        horizontal = cos(heading)*north + sin(heading)*east;
        r_fixed = (cfg.earth.radius_m+c.altitude_m)*up;
        v_relative = c.speed_m_s * ...
            (sin(gamma)*up + cos(gamma)*horizontal);
        v_inertial = v_relative + cross(omega_vec, r_fixed);

        X = zeros(14,1);
        X(1:3) = r_fixed;
        X(4:6) = v_inertial;
        X(7:10) = [0; 0; 0; 1];
        X(14) = mass_kg;

        states(k).case_id = c.id;
        states(k).classification = "EXACT_SPHERICAL_STATE_WITH_SURROGATE_ECI_ADAPTER";
        states(k).altitude_m = c.altitude_m;
        states(k).longitude_deg = c.longitude_deg;
        states(k).latitude_deg = c.latitude_deg;
        states(k).speed_m_s = c.speed_m_s;
        states(k).fpa_deg = c.fpa_deg;
        states(k).heading_deg = c.heading_deg;
        states(k).r_earth_fixed_m = r_fixed;
        states(k).v_atmosphere_relative_m_s = v_relative;
        states(k).v_inertial_m_s = v_inertial;
        states(k).X_eci14 = X;
        states(k).mass_kg = mass_kg;
        states(k).mass_classification = "SURROGATE";
    end
end
