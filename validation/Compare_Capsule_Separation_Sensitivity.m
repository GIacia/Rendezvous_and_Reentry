function comparison = Compare_Capsule_Separation_Sensitivity(attached_mass_kg, duration_s)
    % Compare a separated 60 kg capsule with an attached high-mass case.
    %
    % Both cases intentionally use the same capsule reference area, Cd, and
    % L/D. The attached case is therefore a ballistic-coefficient sensitivity
    % case, not a validated aerodynamic model of the chaser-capsule stack.

    project_root = fileparts(fileparts(mfilename('fullpath')));
    addpath(project_root);
    sys = Mission_Config();

    if nargin < 1 || isempty(attached_mass_kg)
        attached_mass_kg = sys.Chaser_Mass_Init;
        if sys.reentry_vehicle.capsule.add_to_chaser_initial_mass
            attached_mass_kg = attached_mass_kg + sys.reentry_vehicle.capsule.mass_kg;
        end
    end
    if nargin < 2 || isempty(duration_s)
        duration_s = 600;
    end

    X0 = paper_capsule_entry_state_local(sys, attached_mass_kg);
    common = struct();
    common.vehicle_mode = "CAPSULE";
    common.dt = 1;
    common.max_time = duration_s;
    common.terminal_altitude = 20e3;
    common.lift_enabled = true;
    common.bank_angle_deg = 0;

    separated = common;
    separated.separation_mode = "ENTRY_INTERFACE";
    separated.capsule_mass_kg = sys.reentry_vehicle.capsule.mass_kg;
    [X_sep, ~, hist_sep, summary_sep] = Reentry_Propagator(sys, X0, [], 0, separated);

    attached = common;
    attached.separation_mode = "ATTACHED";
    [X_att, ~, hist_att, summary_att] = Reentry_Propagator(sys, X0, [], 0, attached);

    stack_surrogate = common;
    stack_surrogate.separation_mode = "ATTACHED";
    stack_surrogate.reference_area_m2 = sys.environment.atmospheric_drag.chaser_area_m2;
    stack_surrogate.cd = sys.environment.atmospheric_drag.chaser_cd;
    stack_surrogate.nominal_ld = 0;
    stack_surrogate.lift_enabled = false;
    [X_stack, ~, hist_stack, summary_stack] = ...
        Reentry_Propagator(sys, X0, [], 0, stack_surrogate);

    comparison = struct();
    comparison.duration_requested_s = duration_s;
    comparison.separated_mass_kg = X_sep(14);
    comparison.attached_mass_kg = X_att(14);
    comparison.separated_beta_kg_m2 = summary_sep.ballistic_coefficient_kg_m2;
    comparison.attached_beta_kg_m2 = summary_att.ballistic_coefficient_kg_m2;
    comparison.common_comparison_time_s = min(hist_sep.time(end), hist_att.time(end));
    sep_pos_common = interp1(hist_sep.time', hist_sep.rv_pos', comparison.common_comparison_time_s, 'linear')';
    att_pos_common = interp1(hist_att.time', hist_att.rv_pos', comparison.common_comparison_time_s, 'linear')';
    sep_vel_common = interp1(hist_sep.time', hist_sep.rv_vel', comparison.common_comparison_time_s, 'linear')';
    att_vel_common = interp1(hist_att.time', hist_att.rv_vel', comparison.common_comparison_time_s, 'linear')';
    comparison.final_position_difference_m = norm(sep_pos_common - att_pos_common);
    comparison.final_velocity_difference_m_s = norm(sep_vel_common - att_vel_common);
    comparison.final_speed_difference_m_s = abs(norm(sep_vel_common) - norm(att_vel_common));
    comparison.separated_altitude_at_common_time_m = norm(sep_pos_common) - sys.Re;
    comparison.attached_altitude_at_common_time_m = norm(att_pos_common) - sys.Re;
    comparison.separated_final_altitude_m = summary_sep.final_altitude_m;
    comparison.attached_final_altitude_m = summary_att.final_altitude_m;
    comparison.separated_peak_dynamic_pressure_Pa = summary_sep.max_dynamic_pressure_Pa;
    comparison.attached_peak_dynamic_pressure_Pa = summary_att.max_dynamic_pressure_Pa;
    comparison.separated_peak_heat_flux_W_m2 = summary_sep.max_heat_flux_W_m2;
    comparison.attached_peak_heat_flux_W_m2 = summary_att.max_heat_flux_W_m2;
    comparison.separated_termination_reason = summary_sep.termination_reason;
    comparison.attached_termination_reason = summary_att.termination_reason;
    comparison.stack_surrogate_mass_kg = X_stack(14);
    comparison.stack_surrogate_area_m2 = stack_surrogate.reference_area_m2;
    comparison.stack_surrogate_cd = stack_surrogate.cd;
    comparison.stack_surrogate_beta_kg_m2 = summary_stack.ballistic_coefficient_kg_m2;
    comparison.stack_surrogate_common_time_s = min(hist_sep.time(end), hist_stack.time(end));
    sep_pos_stack_time = interp1(hist_sep.time', hist_sep.rv_pos', ...
                                 comparison.stack_surrogate_common_time_s, 'linear')';
    stack_pos_common = interp1(hist_stack.time', hist_stack.rv_pos', ...
                               comparison.stack_surrogate_common_time_s, 'linear')';
    sep_vel_stack_time = interp1(hist_sep.time', hist_sep.rv_vel', ...
                                 comparison.stack_surrogate_common_time_s, 'linear')';
    stack_vel_common = interp1(hist_stack.time', hist_stack.rv_vel', ...
                               comparison.stack_surrogate_common_time_s, 'linear')';
    comparison.stack_surrogate_position_difference_m = norm(sep_pos_stack_time - stack_pos_common);
    comparison.stack_surrogate_velocity_difference_m_s = norm(sep_vel_stack_time - stack_vel_common);
    comparison.separated_altitude_at_stack_common_time_m = norm(sep_pos_stack_time) - sys.Re;
    comparison.stack_surrogate_altitude_at_common_time_m = norm(stack_pos_common) - sys.Re;
    comparison.stack_surrogate_termination_reason = summary_stack.termination_reason;
    comparison.separated_history = hist_sep;
    comparison.attached_history = hist_att;
    comparison.stack_surrogate_history = hist_stack;

    fprintf('Capsule separation sensitivity (same assumed capsule aerodynamics)\n');
    fprintf('  mass / beta: %.1f kg / %.2f kg/m^2 versus %.1f kg / %.2f kg/m^2\n', ...
            comparison.separated_mass_kg, comparison.separated_beta_kg_m2, ...
            comparison.attached_mass_kg, comparison.attached_beta_kg_m2);
    fprintf('  common-time |dr|: %.3f km, |dv|: %.3f m/s at %.1f s\n', ...
            comparison.final_position_difference_m/1e3, ...
            comparison.final_velocity_difference_m_s, comparison.common_comparison_time_s);
    fprintf('  common-time altitude: separated %.3f km, attached %.3f km\n', ...
            comparison.separated_altitude_at_common_time_m/1e3, ...
            comparison.attached_altitude_at_common_time_m/1e3);
    fprintf('  individual endpoints: %s at %.1f s versus %s at %.1f s\n', ...
            char(comparison.separated_termination_reason), hist_sep.time(end), ...
            char(comparison.attached_termination_reason), hist_att.time(end));
    fprintf('  peak q: separated %.3f kPa, attached %.3f kPa\n', ...
            comparison.separated_peak_dynamic_pressure_Pa/1e3, ...
            comparison.attached_peak_dynamic_pressure_Pa/1e3);
    fprintf('  current-chaser CdA surrogate: %.1f kg, Cd %.2f, A %.2f m^2, beta %.2f kg/m^2\n', ...
            comparison.stack_surrogate_mass_kg, comparison.stack_surrogate_cd, ...
            comparison.stack_surrogate_area_m2, comparison.stack_surrogate_beta_kg_m2);
    fprintf('  separated versus CdA surrogate at %.1f s: |dr| %.3f km, |dv| %.3f m/s, altitude %.3f/%.3f km\n', ...
            comparison.stack_surrogate_common_time_s, ...
            comparison.stack_surrogate_position_difference_m/1e3, ...
            comparison.stack_surrogate_velocity_difference_m_s, ...
            comparison.separated_altitude_at_stack_common_time_m/1e3, ...
            comparison.stack_surrogate_altitude_at_common_time_m/1e3);
end

function X = paper_capsule_entry_state_local(sys, mass_kg)
    radius = sys.Re + sys.reentry_vehicle.capsule.entry_interface_altitude_m;
    speed = sys.reentry_vehicle.capsule.reference_entry_speed_m_s;
    fpa = deg2rad(sys.reentry_vehicle.capsule.reference_entry_fpa_deg);

    X = zeros(14,1);
    X(1:3) = [radius;0;0];
    X(4:6) = [speed*sin(fpa);speed*cos(fpa);0];
    X(7:10) = [0;0;0;1];
    X(14) = mass_kg;
end
