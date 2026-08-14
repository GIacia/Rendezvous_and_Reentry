function [clear, clearance, elevation_deg] = line_of_sight_geometry(r_relay, r_rv, sys, margin_altitude)
    los = r_relay - r_rv;
    los_norm = norm(los);
    if los_norm <= eps
        clearance = norm(r_rv) - sys.Re;
        clear = clearance > margin_altitude;
        elevation_deg = 90;
        return;
    end

    u = r_relay - r_rv;
    tau = -dot(r_rv, u) / dot(u, u);
    tau = min(1, max(0, tau));
    closest = r_rv + tau * u;
    clearance = norm(closest) - sys.Re;
    clear = clearance > margin_altitude;

    up = r_rv / norm(r_rv);
    elevation_deg = rad2deg(asin(dot(los / los_norm, up)));
end
