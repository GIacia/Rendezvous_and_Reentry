function a_drag = Atmospheric_Drag_Acceleration(r, v, mass, sys, vehicle_role)
    % Atmospheric drag acceleration in ECI coordinates.
    %
    % Inputs use the simulator's native MATLAB units:
    %   r [m], v [m/s], mass [kg], output [m/s^2].
    %
    % The ISA76 model uses MATLAB atmosisa when available below 86 km and a
    % log-interpolated 1976 standard atmosphere density table above it.

    if nargin < 5 || isempty(vehicle_role)
        vehicle_role = "chaser";
    end

    a_drag = zeros(3,1);

    drag = get_drag_config_local(sys);
    if isempty(drag) || ~get_drag_bool_local(drag, 'enabled', false)
        return;
    end

    model = upper(string(get_drag_field_local(drag, 'model', "ISA76")));
    if model == "OFF" || model == "NONE" || model == "DISABLED"
        return;
    end

    if mass <= 0
        return;
    end

    role = lower(string(vehicle_role));
    if role == "target"
        cd = get_drag_field_local(drag, 'target_cd', get_drag_field_local(drag, 'cd', 2.2));
        area = get_drag_field_local(drag, 'target_area_m2', get_drag_field_local(drag, 'area_m2', 4.0));
    else
        cd = get_drag_field_local(drag, 'chaser_cd', get_drag_field_local(drag, 'cd', 2.2));
        area = get_drag_field_local(drag, 'chaser_area_m2', get_drag_field_local(drag, 'area_m2', 4.0));
    end

    if cd <= 0 || area <= 0
        return;
    end

    altitude_m = norm(r) - sys.Re;
    rho = Standard_Atmosphere_Density(altitude_m, drag);
    if rho <= 0 || ~isfinite(rho)
        return;
    end

    if get_drag_bool_local(drag, 'co_rotate_atmosphere', true)
        omega = get_drag_field_local(drag, 'earth_rotation_rad_s', 7.2921159e-5);
        if isscalar(omega)
            omega_vec = [0; 0; omega];
        else
            omega_vec = omega(:);
        end
        v_atm = cross(omega_vec, r(:));
    else
        v_atm = zeros(3,1);
    end

    v_rel = v(:) - v_atm;
    v_rel_norm = norm(v_rel);
    if v_rel_norm == 0
        return;
    end

    a_drag = -0.5 * rho * cd * area / mass * v_rel_norm * v_rel;
end

function drag = get_drag_config_local(sys)
    drag = [];
    if isstruct(sys) && isfield(sys, 'environment') && ...
            isfield(sys.environment, 'atmospheric_drag')
        drag = sys.environment.atmospheric_drag;
    end
end

function value = get_drag_field_local(s, name, default_value)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    else
        value = default_value;
    end
end

function value = get_drag_bool_local(s, name, default_value)
    value = get_drag_field_local(s, name, default_value);
    if islogical(value)
        return;
    end
    if isnumeric(value)
        value = value ~= 0;
        return;
    end

    txt = lower(strtrim(string(value)));
    value = txt == "true" || txt == "1" || txt == "yes" || txt == "on";
end
