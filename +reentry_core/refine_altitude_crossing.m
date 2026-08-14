function [X_cross, t_cross] = refine_altitude_crossing(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, target_altitude, dt_window)
%REFINE_ALTITUDE_CROSSING Refine a bracketed downward altitude crossing.
    validateattributes(dt_window, {'numeric'}, ...
        {'real','scalar','finite','positive'}, mfilename, 'dt_window');
    validateattributes(target_altitude, {'numeric'}, ...
        {'real','scalar','finite'}, mfilename, 'target_altitude');
    if ~isnumeric(X0) || numel(X0) ~= 7 || any(~isfinite(X0(:)))
        error('reentry_core:refineAltitude:InvalidState', ...
              'X0 must be a finite seven-element state.');
    end
    X0 = X0(:);

    lo = 0;
    hi = dt_window;
    altitude_lo = norm(X0(1:3)) - sys.Re;
    X_cross = reentry_core.rk4_step( ...
        X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, hi);
    altitude_hi = norm(X_cross(1:3)) - sys.Re;
    tolerance = 1e-9;
    if abs(altitude_lo-target_altitude) <= tolerance
        X_cross = X0;
        t_cross = 0;
        return;
    end
    if altitude_lo < target_altitude || altitude_hi > target_altitude
        error('reentry_core:refineAltitude:UnbracketedCrossing', ...
              ['The interval must bracket a downward crossing: initial ' ...
               'altitude >= target and final altitude <= target.']);
    end

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        X_mid = reentry_core.rk4_step( ...
            X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, mid);
        altitude_mid = norm(X_mid(1:3)) - sys.Re;

        if altitude_mid > target_altitude
            lo = mid;
        else
            hi = mid;
            X_cross = X_mid;
        end

        if hi - lo <= 1e-7
            break;
        end
    end

    t_cross = hi;
end
