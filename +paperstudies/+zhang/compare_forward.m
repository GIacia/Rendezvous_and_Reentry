function report = compare_forward(hist, summary, cfg)
%COMPARE_FORWARD Compare an open-loop surrogate with published bounds.

    if nargin < 3 || isempty(cfg)
        cfg = paperstudies.zhang.config();
    end
    if isempty(hist.time)
        error('paperstudies:zhang:compareForward:EmptyHistory', ...
              'A nonempty propagation history is required.');
    end

    elapsed = hist.time(end);
    angle = cfg.surrogate.earth_rotation_rad_s*elapsed;
    R_eci_to_fixed = [cos(angle), sin(angle), 0; ...
                     -sin(angle), cos(angle), 0; 0, 0, 1];
    r_fixed = R_eci_to_fixed*hist.rv_pos(:,end);
    radius = norm(r_fixed);
    longitude_deg = mod(atan2d(r_fixed(2),r_fixed(1))+180,360)-180;
    latitude_deg = asind(r_fixed(3)/radius);

    report.classification = "SURROGATE_OPEN_LOOP_FORWARD";
    report.paper_trajectory_reproduced = false;
    report.final.longitude_deg = longitude_deg;
    report.final.latitude_deg = latitude_deg;
    report.final.altitude_m = hist.altitude(end);
    report.final.speed_m_s = hist.speed_rel(end);
    longitude_error_deg = mod(longitude_deg - ...
        cfg.terminal.longitude_deg + 180, 360) - 180;
    report.terminal.longitude_satisfied = abs(longitude_error_deg) <= ...
        cfg.terminal.longitude_tolerance_deg;
    report.terminal.latitude_satisfied = abs(latitude_deg- ...
        cfg.terminal.latitude_deg) <= cfg.terminal.latitude_tolerance_deg;
    report.terminal.altitude_satisfied = abs(hist.altitude(end)- ...
        cfg.terminal.altitude_m) <= cfg.terminal.altitude_tolerance_m;
    report.terminal.speed_satisfied = abs(hist.speed_rel(end)- ...
        cfg.terminal.speed_m_s) <= cfg.terminal.speed_tolerance_m_s;
    terminal_values = struct2cell(report.terminal);
    report.terminal.all_satisfied = all(cellfun(@(x) logical(x), terminal_values));

    report.constraints.dynamic_pressure_satisfied = ...
        summary.max_dynamic_pressure_Pa <= cfg.table3.max_dynamic_pressure_Pa;
    report.constraints.load_factor_satisfied = ...
        summary.max_g_load <= cfg.table3.max_load_factor;
    report.constraints.bank_angle_satisfied = ...
        all(abs(hist.bank_angle_deg) <= cfg.table3.max_bank_angle_deg+1e-10);
    report.constraints.bank_rate_satisfied = ...
        summary.max_bank_rate_deg_s <= cfg.table3.max_bank_rate_deg_s+1e-10;

    active = hist.blackout_zone_active ~= 0;
    valid_raap = active & isfinite(hist.raap_deg);
    report.constraints.bzc_interval_present = any(active);
    report.constraints.bzc_fully_evaluated = any(active) && all(valid_raap(active));
    bzc_any_evaluated_violation = any(valid_raap & ...
        hist.raap_deg > cfg.table3.antenna_beam_half_angle_deg+1e-10);
    report.constraints.bzc_any_evaluated_violation = ...
        bzc_any_evaluated_violation;
    if any(valid_raap)
        report.constraints.max_bzc_raap_deg = max(hist.raap_deg(valid_raap));
    else
        report.constraints.max_bzc_raap_deg = NaN;
    end
    if bzc_any_evaluated_violation
        report.constraints.bzc_satisfied = false;
    elseif report.constraints.bzc_fully_evaluated
        report.constraints.bzc_satisfied = true;
    else
        report.constraints.bzc_satisfied = NaN;
    end
    report.constraints.heating_satisfied = NaN;
    report.constraints.heating_status = "UNAVAILABLE_CONFLICTING_SOURCE_DATA";

    evaluated = [report.constraints.dynamic_pressure_satisfied, ...
                 report.constraints.load_factor_satisfied, ...
                 report.constraints.bank_angle_satisfied, ...
                 report.constraints.bank_rate_satisfied];
    if isequal(report.constraints.bzc_satisfied, false)
        evaluated(end+1) = false;
    elseif isequal(report.constraints.bzc_satisfied, true)
        evaluated(end+1) = true;
    end
    report.constraints.evaluated_constraints_satisfied = all(evaluated);
    if ~report.constraints.evaluated_constraints_satisfied
        report.constraints.all_published_constraints_satisfied = false;
    else
        report.constraints.all_published_constraints_satisfied = NaN;
    end
    report.optimizer.status = "UNAVAILABLE";
    report.optimizer.bank_history_source = "CONSTANT_SURROGATE_NOT_OPTIMIZED";
end
