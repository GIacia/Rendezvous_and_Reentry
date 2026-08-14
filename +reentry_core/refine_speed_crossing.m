function [X_cross, t_cross] = refine_speed_crossing(X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, target_speed, dt_window)
%REFINE_SPEED_CROSSING Refine a bracketed downward relative-speed crossing.
    validateattributes(dt_window, {'numeric'}, ...
        {'real','scalar','finite','positive'}, mfilename, 'dt_window');
    validateattributes(target_speed, {'numeric'}, ...
        {'real','scalar','finite','nonnegative'}, mfilename, 'target_speed');
    if ~isnumeric(X0) || numel(X0) ~= 7 || any(~isfinite(X0(:)))
        error('reentry_core:refineSpeed:InvalidState', ...
              'X0 must be a finite seven-element state.');
    end
    X0 = X0(:);

    lo = 0;
    hi = dt_window;
    aux_lo = reentry_core.evaluate_state( ...
        X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    X_cross = reentry_core.rk4_step( ...
        X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, hi);
    aux_hi = reentry_core.evaluate_state( ...
        X_cross, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    tolerance = 1e-9;
    if abs(aux_lo.speed_rel-target_speed) <= tolerance
        X_cross = X0;
        t_cross = 0;
        return;
    end
    if aux_lo.speed_rel < target_speed || aux_hi.speed_rel > target_speed
        error('reentry_core:refineSpeed:UnbracketedCrossing', ...
              ['The interval must bracket a downward crossing: initial ' ...
               'relative speed >= target and final speed <= target.']);
    end

    for iter = 1:50
        mid = 0.5 * (lo + hi);
        X_mid = reentry_core.rk4_step( ...
            X0, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, mid);
        aux_mid = reentry_core.evaluate_state( ...
            X_mid, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);

        if aux_mid.speed_rel > target_speed
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
