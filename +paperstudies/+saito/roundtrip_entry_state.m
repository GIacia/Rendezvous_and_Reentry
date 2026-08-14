function state = roundtrip_entry_state(X14, cfg)
%ROUNDTRIP_ENTRY_STATE Recover spherical entry coordinates from ECI state.

    if nargin < 2 || isempty(cfg)
        cfg = paperstudies.saito.config();
    end
    if numel(X14) < 6 || any(~isfinite(X14(1:6)))
        error('paperstudies:saito:roundtrip_entry_state:InvalidState', ...
              'X14 must contain at least six finite Cartesian states.');
    end

    position = X14(1:3);
    velocity_eci = X14(4:6);
    radius = norm(position);
    if radius <= 0
        error('paperstudies:saito:roundtrip_entry_state:ZeroRadius', ...
              'Position radius must be positive.');
    end
    r_hat = position / radius;
    latitude_rad = asin(min(1, max(-1, r_hat(3))));
    longitude_rad = atan2(r_hat(2), r_hat(1));
    north_hat = [-sin(latitude_rad) * cos(longitude_rad); ...
                 -sin(latitude_rad) * sin(longitude_rad); ...
                  cos(latitude_rad)];
    east_hat = [-sin(longitude_rad); cos(longitude_rad); 0];

    omega_earth = [0; 0; cfg.conventions.earth_rotation_rad_s];
    velocity_relative = velocity_eci - cross(omega_earth, position);
    speed_relative = norm(velocity_relative);
    if speed_relative <= 0
        error('paperstudies:saito:roundtrip_entry_state:ZeroSpeed', ...
              'Atmosphere-relative speed must be positive.');
    end
    radial_speed = dot(velocity_relative, r_hat);
    horizontal_velocity = velocity_relative - radial_speed * r_hat;
    fpa_rad = asin(min(1, max(-1, radial_speed / speed_relative)));
    azimuth_rad = atan2(dot(horizontal_velocity, east_hat), ...
                        dot(horizontal_velocity, north_hat));

    state.altitude_m = radius - cfg.conventions.earth_radius_m;
    state.longitude_deg = rad2deg(longitude_rad);
    state.latitude_deg = rad2deg(latitude_rad);
    state.speed_relative_m_s = speed_relative;
    state.fpa_deg = rad2deg(fpa_rad);
    state.azimuth_deg = rad2deg(azimuth_rad);
end
