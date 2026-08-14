function dir = lift_direction(r, v_rel, bank_angle_deg)
    dir = zeros(3,1);
    h = cross(r, v_rel);
    if norm(h) <= eps || norm(v_rel) <= eps
        return;
    end

    vhat = v_rel / norm(v_rel);
    hhat = h / norm(h);
    dir = cross(vhat, hhat);
    if norm(dir) <= eps
        dir = zeros(3,1);
        return;
    end
    dir = dir / norm(dir);

    bank = deg2rad(bank_angle_deg);
    dir = dir*cos(bank) + cross(vhat, dir)*sin(bank) + ...
        vhat*dot(vhat, dir)*(1-cos(bank));
    dir = dir / norm(dir);
end
