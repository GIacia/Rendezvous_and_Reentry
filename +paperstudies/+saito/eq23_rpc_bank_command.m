function [sigma_command_rad, cos_sigma_command] = ...
    eq23_rpc_bank_command(sigma_current_rad, delta_sigma_rad, ...
                          R_current, R_predict1, R_predict2)
%EQ23_RPC_BANK_COMMAND Saito Eq. (23), principal unsigned bank result.

    denominator = R_predict2 - R_predict1;
    if any(~isfinite(denominator(:))) || any(abs(denominator(:)) <= eps)
        error('paperstudies:saito:eq23:SingularSensitivity', ...
              'R_predict2 - R_predict1 must be finite and nonzero.');
    end
    cos_current = cos(sigma_current_rad);
    cos_perturbed = cos(sigma_current_rad + delta_sigma_rad);
    cos_sigma_command = cos_current + ...
        (cos_perturbed - cos_current) ./ denominator .* ...
        (R_current - R_predict1);

    tolerance = 100 * eps(max(1, max(abs(cos_sigma_command(:)))));
    if any(abs(cos_sigma_command(:)) > 1 + tolerance)
        error('paperstudies:saito:eq23:InfeasibleCommand', ...
              'Eq. (23) produced a cosine outside [-1,1].');
    end
    cos_sigma_command = min(1, max(-1, cos_sigma_command));
    sigma_command_rad = acos(cos_sigma_command);
end
