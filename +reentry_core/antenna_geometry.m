function antenna = antenna_geometry(X_rv, X_relay, aux, shape, los_clear)
    antenna = struct('raap_deg', NaN, 'range_m', NaN, 'constraint_active', false, ...
                     'blackout_active', false, ...
                     'in_beam', false, 'within_cone', false, ...
                     'within_range', false, 'communication_available', false, ...
                     'best_raap_deg', NaN, ...
                     'best_bank_deg', NaN);
    if ~get_bool(shape, 'antenna_enabled', false)
        return;
    end

    if isempty(X_relay) || numel(X_relay) < 3 || aux.speed_rel <= eps
        return;
    end

    boresight = get_field(shape, 'antenna_boresight_body', [-1;0;0]);
    antenna.range_m = norm(X_relay(1:3) - X_rv(1:3));
    antenna.raap_deg = reentry_core.raap_deg( ...
        X_rv(1:3), aux.v_rel, X_relay(1:3), ...
        aux.aoa_deg, aux.bank_angle_deg, boresight);

    beam_half_angle = get_field(shape, 'antenna_beam_half_angle_deg', 45);
    angle_tolerance = get_field(shape, 'antenna_angle_tolerance_deg', 1e-10);
    min_range = get_field(shape, 'antenna_min_range_m', 0);
    max_range = get_field(shape, 'antenna_max_range_m', inf);
    antenna.within_cone = ...
        antenna.raap_deg <= beam_half_angle + angle_tolerance;
    antenna.in_beam = antenna.within_cone;
    antenna.within_range = ...
        antenna.range_m >= min_range && antenna.range_m <= max_range;
    antenna.communication_available = ...
        antenna.within_cone && antenna.within_range && ...
        isfinite(los_clear) && los_clear ~= 0;

    bank_scan_enabled = get_bool(shape, 'evaluate_bank_feasibility', false);
    if bank_scan_enabled
        constraints = get_field(shape, 'constraints', struct());
        max_bank = get_field(constraints, 'max_bank_angle_deg', 60);
        bank_grid = linspace(-max_bank, max_bank, ...
            max(3, ceil(2*max_bank) + 1));
        raap_grid = nan(size(bank_grid));
        for ii = 1:numel(bank_grid)
            raap_grid(ii) = reentry_core.raap_deg( ...
                X_rv(1:3), aux.v_rel, X_relay(1:3), ...
                aux.aoa_deg, bank_grid(ii), boresight);
        end
        [antenna.best_raap_deg, idx] = min(raap_grid);
        antenna.best_bank_deg = bank_grid(idx);
    else
        antenna.best_raap_deg = NaN;
        antenna.best_bank_deg = NaN;
    end
    if bank_scan_enabled && (~isfinite(los_clear) || los_clear == 0)
        antenna.best_raap_deg = inf;
    end
end
