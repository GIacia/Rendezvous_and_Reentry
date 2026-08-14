function [CD, CL] = eq26_aero_coefficients(CA, CN, alpha_rad)
%EQ26_AERO_COEFFICIENTS Public CA/CN-to-CD/CL algebra from Eq. (26).
%   This helper does not supply the proprietary HSRC CA/CN database.
    CD = CA .* cos(alpha_rad) + CN .* sin(alpha_rad);
    CL = -CA .* sin(alpha_rad) + CN .* cos(alpha_rad);
end
