function estimate_next = eq24_fading_filter(estimate_current, beta, measured_ratio)
%EQ24_FADING_FILTER Saito Eq. (24) for lift or drag scale estimates.
    if ~isscalar(beta) || ~isfinite(beta) || beta < 0 || beta > 1
        error('paperstudies:saito:eq24:InvalidBeta', ...
              'beta must be a finite scalar in [0,1].');
    end
    estimate_next = beta .* estimate_current + (1 - beta) .* measured_ratio;
end
