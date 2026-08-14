function dx_augmented = augmented_dynamics(x_augmented, bank_rate, lift, drag, omega)
%AUGMENTED_DYNAMICS Zhang Eqs. (37)-(38) without an optimizer.
%   The first six derivatives use Eqs. (15)-(23), while the seventh is
%   bank-rate control u. This helper does not provide or infer an optimized
%   bank history.

    if ~isnumeric(x_augmented) || ~isreal(x_augmented) || ...
            numel(x_augmented) ~= 7 || any(~isfinite(x_augmented(:)))
        error('paperstudies:zhang:augmentedDynamics:InvalidState', ...
              'The augmented state must be a finite seven-element vector.');
    end
    if ~isscalar(bank_rate) || ~isfinite(bank_rate) || ~isreal(bank_rate)
        error('paperstudies:zhang:augmentedDynamics:InvalidControl', ...
              'Bank-rate control must be a finite real scalar.');
    end
    x_augmented = x_augmented(:);
    dx_augmented = [paperstudies.zhang.dimensionless_dynamics( ...
        x_augmented(1:6), x_augmented(7), lift, drag, omega); bank_rate];
end
