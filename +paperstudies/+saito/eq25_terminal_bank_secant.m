function sigma_next_rad = eq25_terminal_bank_secant( ...
    sigma_i_rad, sigma_ip1_rad, R_predict_i, R_predict_ip1)
%EQ25_TERMINAL_BANK_SECANT One Saito Eq. (25) secant iteration.

    denominator = R_predict_ip1 - R_predict_i;
    if any(~isfinite(denominator(:))) || any(abs(denominator(:)) <= eps)
        error('paperstudies:saito:eq25:SingularSecant', ...
              'Predicted errors must have a finite nonzero difference.');
    end
    sigma_next_rad = sigma_ip1_rad - R_predict_ip1 .* ...
        (sigma_ip1_rad - sigma_i_rad) ./ denominator;
end
