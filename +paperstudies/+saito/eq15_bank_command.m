function [sigma_command_rad, cos_sigma_command] = ...
    eq15_bank_command(ld_command, ld_vehicle)
%EQ15_BANK_COMMAND Principal bank-angle magnitude from Saito Eq. (15).
%   The paper's separate cross-range logic assigns bank sign; that schedule
%   is unpublished and is intentionally not inferred here.

    if any(~isfinite(ld_vehicle(:))) || any(ld_vehicle(:) == 0)
        error('paperstudies:saito:eq15:InvalidVehicleLD', ...
              'ld_vehicle must contain finite nonzero values.');
    end
    cos_sigma_command = ld_command ./ ld_vehicle;
    tolerance = 100 * eps(max(1, max(abs(cos_sigma_command(:)))));
    if any(abs(cos_sigma_command(:)) > 1 + tolerance)
        error('paperstudies:saito:eq15:InfeasibleCommand', ...
              '|ld_command/ld_vehicle| exceeds one.');
    end
    cos_sigma_command = min(1, max(-1, cos_sigma_command));
    sigma_command_rad = acos(cos_sigma_command);
end
