function dx = dimensionless_dynamics(x, bank_rad, lift, drag, omega)
%DIMENSIONLESS_DYNAMICS Zhang Eqs. (15)-(23).
%   X = [r; longitude; latitude; speed; fpa; heading]. LIFT, DRAG and
%   OMEGA must use the paper's dimensionless normalization. The paper does
%   not publish enough constants to construct those numerical quantities;
%   this helper preserves the exact equation structure for unit tests and
%   future OCP work.

    if ~isnumeric(x) || ~isreal(x) || numel(x) ~= 6 || any(~isfinite(x(:)))
        error('paperstudies:zhang:dynamics:InvalidState', ...
              'x must be a finite six-element real vector.');
    end
    if ~isscalar(bank_rad) || ~isscalar(lift) || ~isscalar(drag) || ...
            ~isscalar(omega)
        error('paperstudies:zhang:dynamics:InvalidInput', ...
              'Bank, lift, drag, and omega must be scalar values.');
    end
    values = [bank_rad, lift, drag, omega];
    if any(~isfinite(values)) || ~isreal(values)
        error('paperstudies:zhang:dynamics:InvalidInput', ...
              'Bank, lift, drag, and omega must be finite real scalars.');
    end

    x = x(:);
    r = x(1);
    phi = x(3);
    v = x(4);
    gamma = x(5);
    psi = x(6);
    if r <= 0 || v <= 0 || abs(cos(phi)) < 1e-12 || abs(cos(gamma)) < 1e-12
        error('paperstudies:zhang:dynamics:SingularState', ...
              'State is singular for the published spherical equations.');
    end

    p_v = omega^2*r*cos(phi) * ...
        (sin(gamma)*cos(phi) - cos(gamma)*sin(phi)*cos(psi));
    p_gamma = 2*omega*cos(phi)*sin(psi) + ...
        omega^2*r*cos(phi) * ...
        (cos(gamma)*cos(phi) + sin(gamma)*sin(phi)*cos(psi)) / v;
    p_psi = -2*omega*(tan(gamma)*cos(phi)*cos(psi) - sin(phi)) + ...
        omega^2*r*sin(phi)*cos(phi)*sin(psi) / (v*cos(gamma));

    dx = zeros(6,1);
    dx(1) = v*sin(gamma);
    dx(2) = v*cos(gamma)*sin(psi)/(r*cos(phi));
    dx(3) = v*cos(gamma)*cos(psi)/r;
    dx(4) = -drag - sin(gamma)/r^2 + p_v;
    dx(5) = lift*cos(bank_rad)/v - ...
        (1/r-v^2)*cos(gamma)/(v*r) + p_gamma;
    dx(6) = lift*sin(bank_rad)/(v*cos(gamma)) + ...
        v*cos(gamma)*sin(psi)*tan(phi)/r + p_psi;
end
