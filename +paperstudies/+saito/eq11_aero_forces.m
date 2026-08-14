function [lift_N, drag_N] = eq11_aero_forces( ...
    density_kg_m3, speed_m_s, reference_area_m2, CL, CD)
%EQ11_AERO_FORCES Saito Eq. (11).
    dynamic_pressure_Pa = 0.5 .* density_kg_m3 .* speed_m_s.^2;
    lift_N = dynamic_pressure_Pa .* CL .* reference_area_m2;
    drag_N = dynamic_pressure_Pa .* CD .* reference_area_m2;
end
