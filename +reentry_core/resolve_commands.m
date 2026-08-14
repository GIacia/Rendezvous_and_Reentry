function [aoa_deg, bank_deg] = resolve_commands(shape, aoa_fallback_deg, bank_fallback_deg, speed_rel, mach)
    profile_mode = upper(string(get_field(shape, 'aoa_profile_mode', "CONSTANT")));
    bank_deg = bank_fallback_deg;

    if profile_mode == "SPEED_SCHEDULE"
        speed_grid = get_field(shape, 'aoa_speed_grid_m_s', [0 8000]);
        aoa_grid = get_field(shape, 'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        speed_query = min(max(speed_rel, speed_grid(1)), speed_grid(end));
        aoa_deg = interp1(speed_grid, aoa_grid, speed_query, 'linear');
    elseif profile_mode == "MACH_SCHEDULE"
        mach_grid = get_field(shape, 'aoa_mach_grid', [0 30]);
        aoa_grid = get_field(shape, 'aoa_values_deg', [aoa_fallback_deg aoa_fallback_deg]);
        mach_query = min(max(mach, mach_grid(1)), mach_grid(end));
        aoa_deg = interp1(mach_grid, aoa_grid, mach_query, 'linear');
    else
        aoa_deg = aoa_fallback_deg;
    end
end
