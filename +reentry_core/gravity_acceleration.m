function a = gravity_acceleration(r, sys, shape)
    r_norm = norm(r);
    a_g = -sys.mu / r_norm^3 * r;
    gravity_model = upper(string(get_field(shape, 'gravity_model', "CENTRAL_SPHERICAL")));
    if gravity_model == "CENTRAL_SPHERICAL" || gravity_model == "CENTRAL"
        a = a_g;
        return;
    end

    z2 = (r(3)/r_norm)^2;
    factor = 1.5 * sys.J2 * (sys.mu/r_norm^2) * (sys.Re/r_norm)^2;
    a_j2 = factor * [ (r(1)/r_norm)*(5*z2 - 1); ...
                      (r(2)/r_norm)*(5*z2 - 1); ...
                      (r(3)/r_norm)*(5*z2 - 3) ];
    a = a_g + a_j2;
end
