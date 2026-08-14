function [cd, cl, ld] = aerodynamic_coefficients(shape, aoa_deg, mach)
    aero_model = upper(string(get_field(shape, 'aero_model', "LEGACY_SHAPE_LD")));

    if aero_model == "PAPER_RLV_POLYNOMIAL"
        cl_coeff = get_field(shape, 'cl_polynomial', [-0.041065, 0.016292, 0.0002602]);
        cd_coeff = get_field(shape, 'cd_from_cl_polynomial', [0.080505, -0.03026, 0.86495]);
        cl = cl_coeff(1) + cl_coeff(2)*aoa_deg + cl_coeff(3)*aoa_deg^2;
        cd = cd_coeff(1) + cd_coeff(2)*cl + cd_coeff(3)*cl^2;
        cd = max(cd, 1e-6);
        ld = cl / cd;
    elseif aero_model == "PAPER_CAPSULE_REDUCED"
        cd = get_field(shape, 'cd', 1.3);
        ld = get_field(shape, 'nominal_ld', 0.25);
    else
        cd = get_field(shape, 'cd', 1.2);
        ld = lookup_ld(shape, aoa_deg, mach);
    end

    cd = cd * get_field(shape, 'cd_scale', 1.0);
    ld = ld * get_field(shape, 'ld_scale', 1.0);
    cl = cd * ld;
end
