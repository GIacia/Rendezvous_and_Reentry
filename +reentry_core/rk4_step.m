function X_next = rk4_step(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k, dt)
    validateattributes(dt, {'numeric'}, ...
        {'real','scalar','finite','nonnegative'}, mfilename, 'dt');
    X = X(:);
    f = @(X_now) reentry_core.dynamics( ...
        X_now, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    k1 = f(X);
    k2 = f(X + 0.5 * dt * k1);
    k3 = f(X + 0.5 * dt * k2);
    k4 = f(X + dt * k3);
    X_next = X + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end
