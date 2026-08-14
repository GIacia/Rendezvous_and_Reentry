function [R_current, central_angle_rad] = eq14_range_to_go( ...
    latitude1_rad, longitude1_rad, latitude2_rad, longitude2_rad, radius)
%EQ14_RANGE_TO_GO Saito Eq. (14) and optional distance conversion.
%   Without RADIUS, R_CURRENT is the central angle printed by Eq. (14).
%   With RADIUS, R_CURRENT is RADIUS times that angle. The latter
%   convention is needed for the paper's km-valued guidance tables but is
%   not written in Eq. (14).

    argument = sin(latitude1_rad) .* sin(latitude2_rad) + ...
        cos(latitude1_rad) .* cos(latitude2_rad) .* ...
        cos(longitude2_rad - longitude1_rad);
    argument = min(1, max(-1, argument));
    central_angle_rad = acos(argument);
    if nargin < 5 || isempty(radius)
        R_current = central_angle_rad;
    else
        R_current = radius .* central_angle_rad;
    end
end
