function ld = lookup_ld(shape, aoa_deg, ~)
    aoa_grid = get_field(shape, 'ld_aoa_deg', [0 20 80]);
    ld_grid = get_field(shape, 'ld_values', [0 1 0]);
    aoa = min(max(aoa_deg, min(aoa_grid)), max(aoa_grid));
    ld = interp1(aoa_grid, ld_grid, aoa, 'pchip');
    ld = max(0, ld);
end
