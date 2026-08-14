function raap_deg = raap_deg(r, v_rel, r_relay, aoa_deg, bank_deg, boresight_body)
    raap_deg = NaN;
    if norm(v_rel) <= eps || norm(r_relay - r) <= eps
        return;
    end

    x_t = v_rel / norm(v_rel);
    r_hat = r / norm(r);
    y_t = r_hat - dot(r_hat, x_t) * x_t;
    if norm(y_t) <= eps
        return;
    end
    y_t = y_t / norm(y_t);
    z_t = cross(x_t, y_t);
    z_t = z_t / norm(z_t);

    w_t = [x_t, y_t, z_t]' * (r_relay - r);
    alpha = deg2rad(aoa_deg);
    sigma = deg2rad(bank_deg);
    g_alpha = [cos(alpha), sin(alpha), 0; ...
              -sin(alpha), cos(alpha), 0; ...
               0, 0, 1];
    g_bank = [1, 0, 0; ...
              0, cos(sigma), sin(sigma); ...
              0, -sin(sigma), cos(sigma)];
    w_body = g_alpha * g_bank * w_t;

    boresight = boresight_body(:);
    cosine = dot(boresight, w_body) / (norm(boresight) * norm(w_body));
    raap_deg = rad2deg(acos(min(1, max(-1, cosine))));
end
