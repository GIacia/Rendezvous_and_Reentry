function fpa_deg = flight_path_angle_deg(r, v)
    r_hat = r / norm(r);
    v_radial = dot(v, r_hat);
    v_horizontal = norm(v - v_radial * r_hat);
    fpa_deg = rad2deg(atan2(v_radial, v_horizontal));
end
