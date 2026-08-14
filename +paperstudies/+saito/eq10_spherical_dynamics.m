function derivative = eq10_spherical_dynamics( ...
    state, forces_wind_N, mass_kg, gravity_m_s2, earth_rotation_rad_s)
%EQ10_SPHERICAL_DYNAMICS Saito Eq. (10) in a rotating spherical frame.
%   STATE = [r; phi; theta; v; gamma; psi] uses SI units and radians.
%   FORCES_WIND_N = [D_W; F_W; L_W]. Gravity is supplied explicitly so
%   callers can choose constant g or mu/r^2 without hiding a convention.

    if numel(state) ~= 6 || numel(forces_wind_N) ~= 3
        error('paperstudies:saito:eq10:InvalidDimensions', ...
              'state must have six elements and forces_wind_N three.');
    end
    if ~isscalar(mass_kg) || ~isfinite(mass_kg) || mass_kg <= 0
        error('paperstudies:saito:eq10:InvalidMass', ...
              'mass_kg must be a finite positive scalar.');
    end

    r = state(1);
    phi = state(2);
    v = state(4);
    gamma = state(5);
    psi = state(6);
    D_w = forces_wind_N(1);
    F_w = forces_wind_N(2);
    L_w = forces_wind_N(3);
    omega = earth_rotation_rad_s;

    if r <= 0 || v <= 0 || abs(cos(phi)) <= 1e-12 || abs(cos(gamma)) <= 1e-12
        error('paperstudies:saito:eq10:SingularState', ...
              'Eq. (10) requires positive r/v and nonsingular cos(phi/gamma).');
    end

    r_dot = v * sin(gamma);
    phi_dot = v * cos(gamma) * cos(psi) / r;
    theta_dot = v * cos(gamma) * sin(psi) / (r * cos(phi));

    v_dot = -D_w / mass_kg - gravity_m_s2 * sin(gamma) - ...
        omega^2 * r * cos(phi) * ...
        (cos(gamma) * cos(psi) * sin(phi) - sin(gamma) * cos(phi));

    psi_dot = (F_w + mass_kg * v^2 / r * cos(gamma)^2 * ...
        sin(psi) * tan(phi) + mass_kg * omega^2 * r * ...
        sin(psi) * sin(phi) * cos(phi) - 2 * mass_kg * omega * v * ...
        (sin(gamma) * cos(psi) * cos(phi) - cos(gamma) * sin(phi))) / ...
        (mass_kg * v * cos(gamma));

    gamma_dot = (L_w - mass_kg * gravity_m_s2 * cos(gamma) + ...
        mass_kg * v^2 / r * cos(gamma) + ...
        2 * mass_kg * omega * v * sin(psi) * cos(phi) + ...
        mass_kg * omega^2 * r * cos(phi) * ...
        (sin(gamma) * cos(psi) * sin(phi) + cos(gamma) * cos(phi))) / ...
        (mass_kg * v);

    derivative = [r_dot; phi_dot; theta_dot; v_dot; gamma_dot; psi_dot];
end
