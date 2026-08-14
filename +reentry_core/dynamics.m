function dX = dynamics(X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k)
    X = X(:);
    aux = reentry_core.evaluate_state( ...
        X, sys, shape, aoa_deg, bank_angle_deg, lift_enabled, heat_k);
    dX = [X(4:6); aux.a_gravity + aux.a_aero; 0];
end
