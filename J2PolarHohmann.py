import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar, minimize_scalar


# ============================================================
# Constants
# ============================================================

@dataclass
class EarthJ2:
    mu: float = 398600.4418          # km^3/s^2
    radius: float = 6378.137         # km
    J2: float = 1.08262668e-3        
    # J2: float = 0.0


G0_M_S2 = 9.80665


# ============================================================
# Basic utilities
# ============================================================

def wrap_to_2pi(angle):
    return angle % (2.0 * np.pi)


def signed_angle_error(angle, reference):
    """
    Returns angle - reference wrapped to [-pi, pi].
    """
    return np.arctan2(np.sin(angle - reference), np.cos(angle - reference))


# Polar orbital plane:
# ECI x-z plane.
# r = +x at argument of latitude u = 0.
# motion is +x -> +z, so the spacecraft passes over the north pole.
P_HAT = np.array([1.0, 0.0, 0.0])
Q_HAT = np.array([0.0, 0.0, 1.0])


def argument_of_latitude_in_plane(r_eci):
    """
    Argument of latitude measured inside the fixed polar x-z orbital plane.
    """
    x_comp = np.dot(r_eci, P_HAT)
    z_comp = np.dot(r_eci, Q_HAT)
    return np.arctan2(z_comp, x_comp)


def circular_polar_state(radius_km, u_rad, earth=EarthJ2()):
    """
    Circular polar orbit state in ECI frame.

    Orbit plane: x-z plane.
    Direction: +x -> +z.
    """
    r = radius_km * (np.cos(u_rad) * P_HAT + np.sin(u_rad) * Q_HAT)
    
    v_mag_np = np.sqrt(earth.mu / radius_km)
    v_mag = np.sqrt(earth.mu / radius_km * (1 - earth.J2 * (earth.radius / radius_km)**2 * (3*np.sin(u_rad)**2 - 1)))
    v = v_mag * (-np.sin(u_rad) * P_HAT + np.cos(u_rad) * Q_HAT)

    return np.concatenate([r, v])

# ============================================================
# LVLH / ECI coordinate transformation
# ============================================================

def lvlh_basis_from_target_state(target_state_eci):
    """
    Builds target-centered LVLH basis vectors expressed in ECI frame.

    LVLH convention used here:
        x_LVLH : radial outward direction from Earth to target
        y_LVLH : along-track / prograde direction
        z_LVLH : orbit normal direction, r x v

    Parameters
    ----------
    target_state_eci : array-like, shape (6,)
        [r_x, r_y, r_z, v_x, v_y, v_z] in ECI.
        Units: km, km/s.

    Returns
    -------
    C_L2I : ndarray, shape (3, 3)
        Rotation matrix from LVLH coordinates to ECI coordinates.

        rel_eci = C_L2I @ rel_lvlh

    C_I2L : ndarray, shape (3, 3)
        Rotation matrix from ECI coordinates to LVLH coordinates.

        rel_lvlh = C_I2L @ rel_eci
    """
    r_t = np.array(target_state_eci[0:3], dtype=float)
    v_t = np.array(target_state_eci[3:6], dtype=float)

    r_norm = np.linalg.norm(r_t)
    h_vec = np.cross(r_t, v_t)
    h_norm = np.linalg.norm(h_vec)

    if r_norm == 0.0:
        raise ValueError("Target position vector has zero norm.")

    if h_norm == 0.0:
        raise ValueError("Target angular momentum vector has zero norm.")

    x_lvlh_in_eci = r_t / r_norm
    z_lvlh_in_eci = h_vec / h_norm
    y_lvlh_in_eci = np.cross(z_lvlh_in_eci, x_lvlh_in_eci)
    y_lvlh_in_eci = y_lvlh_in_eci / np.linalg.norm(y_lvlh_in_eci)

    C_L2I = np.column_stack([
        x_lvlh_in_eci,
        y_lvlh_in_eci,
        z_lvlh_in_eci
    ])

    C_I2L = C_L2I.T

    return C_L2I, C_I2L


def lvlh_position_to_eci_point(target_state_eci, rel_pos_lvlh):
    """
    Converts a target-centered LVLH relative position into an absolute ECI point.

    Parameters
    ----------
    target_state_eci : array-like, shape (6,)
        Target state in ECI.
        Units: km, km/s.

    rel_pos_lvlh : array-like, shape (3,)
        Desired relative position from target in LVLH frame.
        Units: km.

        Convention:
            [x, y, z] = [radial outward, along-track, orbit normal]

    Returns
    -------
    point_eci : ndarray, shape (3,)
        Absolute ECI position of the desired point.
        Units: km.
    """
    target_state_eci = np.array(target_state_eci, dtype=float)
    rel_pos_lvlh = np.array(rel_pos_lvlh, dtype=float)

    r_t = target_state_eci[0:3]

    C_L2I, _ = lvlh_basis_from_target_state(target_state_eci)

    rel_pos_eci = C_L2I @ rel_pos_lvlh
    point_eci = r_t + rel_pos_eci

    return point_eci


def eci_point_to_lvlh_position(target_state_eci, point_eci):
    """
    Converts an absolute ECI point into target-centered LVLH relative position.

    Parameters
    ----------
    target_state_eci : array-like, shape (6,)
        Target state in ECI.
        Units: km, km/s.

    point_eci : array-like, shape (3,)
        Absolute ECI position of some point.
        Units: km.

    Returns
    -------
    rel_pos_lvlh : ndarray, shape (3,)
        Relative position from target to point in LVLH frame.
        Units: km.

        Convention:
            [x, y, z] = [radial outward, along-track, orbit normal]
    """
    target_state_eci = np.array(target_state_eci, dtype=float)
    point_eci = np.array(point_eci, dtype=float)

    r_t = target_state_eci[0:3]

    _, C_I2L = lvlh_basis_from_target_state(target_state_eci)

    rel_pos_eci = point_eci - r_t
    rel_pos_lvlh = C_I2L @ rel_pos_eci

    return rel_pos_lvlh


def relative_position_eci_to_lvlh(target_state_eci, chaser_position_eci):
    """
    Convenience function:
    Converts chaser ECI position into target-centered LVLH relative position.
    """
    return eci_point_to_lvlh_position(
        target_state_eci=target_state_eci,
        point_eci=chaser_position_eci
    )


def relative_position_lvlh_to_eci(target_state_eci, rel_pos_lvlh):
    """
    Convenience function:
    Converts target-centered LVLH relative position into ECI relative vector.

    Unlike lvlh_position_to_eci_point(), this returns only the relative vector,
    not the absolute ECI point.
    """
    C_L2I, _ = lvlh_basis_from_target_state(target_state_eci)

    rel_pos_lvlh = np.array(rel_pos_lvlh, dtype=float)
    rel_pos_eci = C_L2I @ rel_pos_lvlh

    return rel_pos_eci

# ============================================================
# Dynamics with J2
# ============================================================

def acceleration_j2_eci(r_eci, earth=EarthJ2()):
    """
    Two-body gravity + J2 perturbation in ECI frame.
    Units:
        r: km
        acceleration: km/s^2
    """
    x, y, z = r_eci
    r2 = np.dot(r_eci, r_eci)
    r = np.sqrt(r2)

    a_two_body = -earth.mu * r_eci / r**3

    z_over_r = z / r
    factor = 1.5 * earth.J2 * earth.mu * earth.radius**2 / r**5

    a_j2 = factor * np.array([
        x * (5.0 * z_over_r**2 - 1.0),
        y * (5.0 * z_over_r**2 - 1.0),
        z * (5.0 * z_over_r**2 - 3.0)
    ])

    return a_two_body + a_j2


def target_chaser_ode_j2(t, y, earth=EarthJ2()):
    """
    Combined target-chaser propagation.

    y =
    [
        target r(3), target v(3),
        chaser r(3), chaser v(3)
    ]
    """
    dydt = np.zeros_like(y)

    r_t = y[0:3]
    v_t = y[3:6]
    r_c = y[6:9]
    v_c = y[9:12]

    dydt[0:3] = v_t
    dydt[3:6] = acceleration_j2_eci(r_t, earth)

    dydt[6:9] = v_c
    dydt[9:12] = acceleration_j2_eci(r_c, earth)

    return dydt


# ============================================================
# Non-J2 Hohmann reference values
# ============================================================

def non_j2_hohmann_defaults(
    h_chaser_km=300.0,
    h_target_km=500.0,
    earth=EarthJ2()
):
    """
    Non-J2 circular-to-circular Hohmann reference.

    Chaser starts lower, target is higher.

    phase_angle definition:
        phase_angle = target argument - chaser argument,
        wrapped to [0, 2pi).

    For ideal Hohmann:
        target must be ahead by
        pi - n_target * transfer_time.
    """
    r1 = earth.radius + h_chaser_km
    r2 = earth.radius + h_target_km

    a_transfer = 0.5 * (r1 + r2)
    transfer_time = np.pi * np.sqrt(a_transfer**3 / earth.mu)

    v_circ_1 = np.sqrt(earth.mu / r1)
    v_transfer_perigee = np.sqrt(earth.mu * (2.0 / r1 - 1.0 / a_transfer))

    delta_v_1 = v_transfer_perigee - v_circ_1

    n_target = np.sqrt(earth.mu / r2**3)

    phase_angle = wrap_to_2pi(np.pi - n_target * transfer_time)

    return {
        "r_chaser_initial_km": r1,
        "r_target_initial_km": r2,
        "a_transfer_km": a_transfer,
        "transfer_time_s": transfer_time,
        "delta_v_1_km_s": delta_v_1,
        "phase_angle_rad": phase_angle,
        "phase_angle_deg": np.degrees(phase_angle)
    }


# ============================================================
# Phase angle waiting logic
# ============================================================

def current_phase_angle(y):
    """
    phase angle = target argument - chaser argument,
    wrapped to [0, 2pi).
    """
    r_t = y[0:3]
    r_c = y[6:9]

    u_t = argument_of_latitude_in_plane(r_t)
    u_c = argument_of_latitude_in_plane(r_c)

    return wrap_to_2pi(u_t - u_c)


def find_first_phase_crossing(
    sol,
    t_start,
    t_end,
    desired_phase_angle,
    sample_step_s=60.0
):
    """
    Finds the first time when the target-chaser phase angle reaches desired_phase_angle.

    Uses dense output from solve_ivp.
    """
    ts = np.arange(t_start, t_end + sample_step_s, sample_step_s)

    if ts[-1] > t_end:
        ts[-1] = t_end

    ys = sol.sol(ts)

    u_t_list = []
    u_c_list = []

    for i in range(ys.shape[1]):
        u_t_list.append(argument_of_latitude_in_plane(ys[0:3, i]))
        u_c_list.append(argument_of_latitude_in_plane(ys[6:9, i]))

    u_t_unwrapped = np.unwrap(np.array(u_t_list))
    u_c_unwrapped = np.unwrap(np.array(u_c_list))

    phase_unwrapped = u_t_unwrapped - u_c_unwrapped

    phase_start = phase_unwrapped[0]
    phase_end = phase_unwrapped[-1]

    lo = min(phase_start, phase_end)
    hi = max(phase_start, phase_end)

    k_min = int(np.floor((lo - desired_phase_angle) / (2.0 * np.pi))) - 1
    k_max = int(np.ceil((hi - desired_phase_angle) / (2.0 * np.pi))) + 1

    candidates = []

    for k in range(k_min, k_max + 1):
        branch_phase = desired_phase_angle + 2.0 * np.pi * k

        if not (lo <= branch_phase <= hi):
            continue

        values = phase_unwrapped - branch_phase

        for i in range(len(ts) - 1):
            if values[i] == 0.0 and ts[i] > t_start:
                candidates.append((ts[i], branch_phase, i, i))
                break

            if values[i] * values[i + 1] <= 0.0:
                candidates.append((ts[i], branch_phase, i, i + 1))
                break

    if len(candidates) == 0:
        return None

    candidates.sort(key=lambda item: item[0])

    _, branch_phase, i0, i1 = candidates[0]

    if i0 == i1:
        return ts[i0]

    t_left = ts[i0]
    t_right = ts[i1]

    def root_function(t):
        y = sol.sol(t)

        u_t = argument_of_latitude_in_plane(y[0:3])
        u_c = argument_of_latitude_in_plane(y[6:9])

        phase_now = u_t - u_c

        return signed_angle_error(phase_now, desired_phase_angle)

    try:
        f_left = root_function(t_left)
        f_right = root_function(t_right)

        if f_left * f_right > 0.0:
            # Fallback: linear interpolation on unwrapped sampled phase.
            v_left = phase_unwrapped[i0] - branch_phase
            v_right = phase_unwrapped[i1] - branch_phase
            alpha = abs(v_left) / (abs(v_left) + abs(v_right))
            return t_left + alpha * (t_right - t_left)

        result = root_scalar(
            root_function,
            bracket=[t_left, t_right],
            xtol=1e-9,
            rtol=1e-11
        )

        return result.root

    except Exception:
        # Robust fallback.
        v_left = phase_unwrapped[i0] - branch_phase
        v_right = phase_unwrapped[i1] - branch_phase
        alpha = abs(v_left) / (abs(v_left) + abs(v_right))
        return t_left + alpha * (t_right - t_left)


# ============================================================
# Burn model
# ============================================================

def normalize_burn_model(burn_model):
    model = str(burn_model).strip().lower()

    aliases = {
        "instant": "impulsive",
        "instantaneous": "impulsive",
        "custom_impulse": "impulsive",
        "finite": "finite_burn",
        "finite_impulse": "finite_burn",
        "continuous": "finite_burn",
    }

    model = aliases.get(model, model)

    if model not in {"impulsive", "finite_burn"}:
        raise ValueError("burn_model must be 'impulsive' or 'finite_burn'.")

    return model


def local_tangential_radial_burn_direction(r_c, v_c, gamma_rad):
    r_hat = r_c / np.linalg.norm(r_c)

    h_vec = np.cross(r_c, v_c)
    h_hat = h_vec / np.linalg.norm(h_vec)

    t_hat = np.cross(h_hat, r_hat)
    t_hat = t_hat / np.linalg.norm(t_hat)

    burn_dir = (
        np.cos(gamma_rad) * t_hat
        + np.sin(gamma_rad) * r_hat
    )

    return burn_dir / np.linalg.norm(burn_dir)


def apply_impulsive_burn_to_chaser(y, delta_v_km_s, gamma_rad):
    """
    Applies one impulsive burn to the chaser.

    gamma definition:
        gamma = angle between delta-V and local tangential direction.

        gamma = 0:
            pure prograde tangential burn.

        gamma > 0:
            burn is tilted toward outward radial direction.

        gamma < 0:
            burn is tilted toward inward radial direction.
    """
    y_new = y.copy()

    r_c = y_new[6:9]
    v_c = y_new[9:12]

    delta_v_vec = delta_v_km_s * local_tangential_radial_burn_direction(
        r_c=r_c,
        v_c=v_c,
        gamma_rad=gamma_rad
    )

    y_new[9:12] += delta_v_vec

    return y_new, delta_v_vec


def finite_burn_duration_s(delta_v_km_s, mass_kg, thrust_N, isp_s):
    if delta_v_km_s < 0.0:
        raise ValueError("delta_v must be non-negative.")

    if delta_v_km_s == 0.0:
        return 0.0

    if mass_kg <= 0.0:
        raise ValueError("initial_mass_kg must be positive.")

    if thrust_N <= 0.0:
        raise ValueError("thrust_N must be positive for finite burns.")

    if isp_s <= 0.0:
        raise ValueError("isp_s must be positive for finite burns.")

    exhaust_velocity_m_s = isp_s * G0_M_S2
    delta_v_m_s = delta_v_km_s * 1000.0
    final_mass_kg = mass_kg * np.exp(-delta_v_m_s / exhaust_velocity_m_s)

    return (mass_kg - final_mass_kg) * exhaust_velocity_m_s / thrust_N


def target_chaser_ode_j2_finite_burn(t, y, burn_dir_eci, thrust_N, isp_s, earth=EarthJ2()):
    """
    Combined target/chaser propagation with a finite high-thrust burn.

    y =
    [
        target r(3), target v(3),
        chaser r(3), chaser v(3),
        chaser mass kg
    ]
    """
    dydt = np.zeros_like(y)

    r_t = y[0:3]
    v_t = y[3:6]
    r_c = y[6:9]
    v_c = y[9:12]
    mass_kg = y[12]

    burn_dir = np.array(burn_dir_eci, dtype=float)
    burn_dir = burn_dir / np.linalg.norm(burn_dir)

    thrust_accel_km_s2 = (thrust_N / mass_kg) / 1000.0

    dydt[0:3] = v_t
    dydt[3:6] = acceleration_j2_eci(r_t, earth)

    dydt[6:9] = v_c
    dydt[9:12] = acceleration_j2_eci(r_c, earth) + thrust_accel_km_s2 * burn_dir
    dydt[12] = -thrust_N / (isp_s * G0_M_S2)

    return dydt


def apply_finite_burn_to_chaser(
    y,
    delta_v_km_s,
    gamma_rad,
    thrust_N=1.0,
    isp_s=200.0,
    initial_mass_kg=2000.0,
    rtol=1e-10,
    atol=1e-12,
    max_step_burn_s=5.0,
    earth=EarthJ2()
):
    burn_duration_s = finite_burn_duration_s(
        delta_v_km_s=delta_v_km_s,
        mass_kg=initial_mass_kg,
        thrust_N=thrust_N,
        isp_s=isp_s
    )

    initial_direction = local_tangential_radial_burn_direction(
        r_c=y[6:9],
        v_c=y[9:12],
        gamma_rad=gamma_rad
    )

    if burn_duration_s == 0.0:
        return {
            "y_plus": y.copy(),
            "burn_duration_s": 0.0,
            "final_mass_kg": initial_mass_kg,
            "propellant_used_kg": 0.0,
            "delivered_delta_v_km_s": 0.0,
            "equivalent_delta_v_vec_eci_km_s": np.zeros(3),
            "initial_burn_direction_eci": initial_direction,
            "solution": None
        }

    y_aug_0 = np.concatenate([y, np.array([initial_mass_kg], dtype=float)])
    burn_dir_eci = initial_direction

    sol_burn = solve_ivp(
        fun=lambda t, y_aug: target_chaser_ode_j2_finite_burn(
            t=t,
            y=y_aug,
            burn_dir_eci=burn_dir_eci,
            thrust_N=thrust_N,
            isp_s=isp_s,
            earth=earth
        ),
        t_span=(0.0, burn_duration_s),
        y0=y_aug_0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
        max_step=max_step_burn_s
    )

    if not sol_burn.success:
        raise RuntimeError(f"Finite burn propagation failed: {sol_burn.message}")

    y_aug_plus = sol_burn.y[:, -1]
    final_mass_kg = y_aug_plus[12]
    delivered_delta_v_km_s = (
        isp_s * G0_M_S2 * np.log(initial_mass_kg / final_mass_kg) / 1000.0
    )

    return {
        "y_plus": y_aug_plus[0:12],
        "burn_duration_s": burn_duration_s,
        "final_mass_kg": final_mass_kg,
        "propellant_used_kg": initial_mass_kg - final_mass_kg,
        "delivered_delta_v_km_s": delivered_delta_v_km_s,
        "equivalent_delta_v_vec_eci_km_s": delta_v_km_s * initial_direction,
        "initial_burn_direction_eci": initial_direction,
        "solution": sol_burn
    }


# ============================================================
# Closest approach search
# ============================================================

def find_closest_approach(
    sol,
    t_start,
    t_end,
    sample_count=1000,
    desired_rel_lvlh=None
):
    """
    Finds closest approach to a desired target-centered LVLH point.

    desired_rel_lvlh:
        [radial outward, along-track, orbit normal] in km

        [0, 0, 0]      : target body
        [-5, 0, 0]     : 5 km below target, radial inward
        [0, -5, 0]     : 5 km behind target, along-track backward
    """
    if desired_rel_lvlh is None:
        desired_rel_lvlh = np.array([0.0, 0.0, 0.0])
    else:
        desired_rel_lvlh = np.array(desired_rel_lvlh, dtype=float)

    def error_at_time(t):
        y = sol.sol(t)

        target_state = y[0:6]
        chaser_state = y[6:12]

        r_t = target_state[0:3]
        v_t = target_state[3:6]
        r_c = chaser_state[0:3]
        v_c = chaser_state[3:6]

        C_L2I, _ = lvlh_basis_from_target_state(target_state)

        desired_rel_eci = C_L2I @ desired_rel_lvlh
        desired_point_eci = r_t + desired_rel_eci

        # Velocity of a point fixed in target LVLH frame
        h_vec = np.cross(r_t, v_t)
        omega_eci = h_vec / np.linalg.norm(r_t)**2
        desired_point_velocity_eci = v_t + np.cross(omega_eci, desired_rel_eci)

        position_error_eci = r_c - desired_point_eci
        position_error_lvlh = eci_point_to_lvlh_position(
            target_state_eci=target_state,
            point_eci=r_c
        ) - desired_rel_lvlh

        velocity_error_eci = v_c - desired_point_velocity_eci

        distance = np.linalg.norm(position_error_lvlh)

        return {
            "y": y,
            "target_state": target_state,
            "chaser_state": chaser_state,
            "desired_point_eci": desired_point_eci,
            "desired_point_velocity_eci": desired_point_velocity_eci,
            "position_error_eci": position_error_eci,
            "position_error_lvlh": position_error_lvlh,
            "velocity_error_eci": velocity_error_eci,
            "distance": distance
        }

    ts = np.linspace(t_start, t_end, sample_count)

    distances = np.array([
        error_at_time(t)["distance"]
        for t in ts
    ])

    idx = int(np.argmin(distances))

    if idx == 0 or idx == len(ts) - 1:
        t_min = ts[idx]
        data_min = error_at_time(t_min)

    else:
        t_left = ts[idx - 1]
        t_right = ts[idx + 1]

        def distance_function(t):
            return error_at_time(t)["distance"]

        result = minimize_scalar(
            distance_function,
            bounds=(t_left, t_right),
            method="bounded",
            options={"xatol": 1e-8}
        )

        t_min = result.x
        data_min = error_at_time(t_min)

    return {
        "t_closest_after_burn_s": t_min,
        "min_distance_km": data_min["distance"],

        "desired_rel_lvlh_km": desired_rel_lvlh,
        "desired_point_eci_km": data_min["desired_point_eci"],
        "desired_point_velocity_eci_km_s": data_min["desired_point_velocity_eci"],

        "position_error_eci_km": data_min["position_error_eci"],
        "position_error_lvlh_km": data_min["position_error_lvlh"],
        "velocity_error_eci_km_s": data_min["velocity_error_eci"],

        # Backward-compatible names
        "relative_position_eci_km": data_min["position_error_eci"],
        "relative_velocity_eci_km_s": data_min["velocity_error_eci"],
        "relative_speed_km_s": np.linalg.norm(data_min["velocity_error_eci"]),

        "target_state_eci_at_closest": data_min["target_state"],
        "chaser_state_eci_at_closest": data_min["chaser_state"]
    }

# ============================================================
# Main propagator
# ============================================================

def run_j2_polar_hohmann_rendezvous(
    phase_angle=None,
    delta_v=None,
    gamma=0.0,
    burn_model="impulsive",
    thrust_N=1.0,
    isp_s=200.0,
    initial_mass_kg=2000.0,
    angle_unit="rad",
    h_target_km=500.0,
    h_chaser_km=300.0,
    initial_phase_angle=0.0,
    initial_chaser_angle=0.0,
    max_wait_time_s=None,
    post_burn_duration_s=None,
    rtol=1e-10,
    atol=1e-12,
    max_step_wait_s=60.0,
    max_step_transfer_s=20.0,
    max_step_burn_s=5.0,
    wait_sample_step_s=60.0,
    closest_sample_count=1200,
    desired_rel_lvlh=None,
    verbose=True,
    earth=EarthJ2()
):
    """
    J2 polar Hohmann-like rendezvous propagator.

    Parameters
    ----------
    phase_angle:
        Desired phase angle at burn.
        Definition:
            phase_angle = target argument - chaser argument,
            wrapped to [0, 2pi).

        If None, uses non-J2 Hohmann phase angle.

    delta_v:
        Chaser maneuver delta-V budget.
        Unit: km/s.
        If None, uses non-J2 Hohmann first impulse.

    gamma:
        Angle between delta-V direction and tangential direction.
        gamma = 0 means pure prograde tangential burn.

    burn_model:
        "impulsive" for instantaneous delta-V, or "finite_burn" to spread
        the same delta-V direction over a short high-thrust burn using
        thrust_N, isp_s, and initial_mass_kg.

    angle_unit:
        "rad" or "deg".

    initial_phase_angle:
        Initial target-chaser phase angle at simulation start.
        Default 0 means target and chaser start on the same ECI radial line,
        but at different altitudes.

    Returns
    -------
    result: dict
    """
    defaults = non_j2_hohmann_defaults(
        h_chaser_km=h_chaser_km,
        h_target_km=h_target_km,
        earth=earth
    )

    if angle_unit.lower() == "deg":
        if phase_angle is not None:
            phase_angle = np.radians(phase_angle)

        gamma = np.radians(gamma)
        initial_phase_angle = np.radians(initial_phase_angle)
        initial_chaser_angle = np.radians(initial_chaser_angle)

    elif angle_unit.lower() != "rad":
        raise ValueError("angle_unit must be either 'rad' or 'deg'.")

    if phase_angle is None:
        phase_angle = defaults["phase_angle_rad"]

    if delta_v is None:
        delta_v = defaults["delta_v_1_km_s"]

    burn_model = normalize_burn_model(burn_model)

    if desired_rel_lvlh is None:
        desired_rel_lvlh = np.array([0, 0, 0])
    else:
        desired_rel_lvlh = np.array(desired_rel_lvlh)

    phase_angle = wrap_to_2pi(phase_angle)

    r_target = earth.radius + h_target_km
    r_chaser = earth.radius + h_chaser_km

    # Initial states
    target_state_0 = circular_polar_state(
        radius_km=r_target,
        u_rad=initial_phase_angle + initial_chaser_angle,
        earth=earth
    )

    chaser_state_0 = circular_polar_state(
        radius_km=r_chaser,
        u_rad=initial_chaser_angle,
        earth=earth
    )

    y0 = np.concatenate([target_state_0, chaser_state_0])

    phase_now = current_phase_angle(y0)

    # If already at desired phase angle, burn at t = 0.
    if abs(signed_angle_error(phase_now, phase_angle)) < 1e-10:
        t_burn = 0.0
        y_burn_minus = y0.copy()

    else:
        n_target = np.sqrt(earth.mu / r_target**3)
        n_chaser = np.sqrt(earth.mu / r_chaser**3)
        synodic_period = 2.0 * np.pi / abs(n_chaser - n_target)

        if max_wait_time_s is None:
            max_wait_time_s = 2.0 * synodic_period

        sol_wait = solve_ivp(
            fun=lambda t, y: target_chaser_ode_j2(t, y, earth),
            t_span=(0.0, max_wait_time_s),
            y0=y0,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
            max_step=max_step_wait_s
        )

        if not sol_wait.success:
            raise RuntimeError(f"Wait propagation failed: {sol_wait.message}")

        t_burn = find_first_phase_crossing(
            sol=sol_wait,
            t_start=0.0,
            t_end=max_wait_time_s,
            desired_phase_angle=phase_angle,
            sample_step_s=wait_sample_step_s
        )

        if t_burn is None:
            raise RuntimeError(
                "Could not find requested phase angle within max_wait_time_s. "
                "Try increasing max_wait_time_s."
            )

        y_burn_minus = sol_wait.sol(t_burn)

    # Apply maneuver
    if burn_model == "impulsive":
        y_burn_plus, delta_v_vec_eci = apply_impulsive_burn_to_chaser(
            y=y_burn_minus,
            delta_v_km_s=delta_v,
            gamma_rad=gamma
        )
        burn_duration_s = 0.0
        delivered_delta_v_km_s = delta_v
        final_mass_kg = initial_mass_kg * np.exp(
            -(delta_v * 1000.0) / (isp_s * G0_M_S2)
        )
        propellant_used_kg = initial_mass_kg - final_mass_kg
        initial_burn_direction_eci = (
            delta_v_vec_eci / np.linalg.norm(delta_v_vec_eci)
            if np.linalg.norm(delta_v_vec_eci) > 0.0
            else np.zeros(3)
        )
        burn_solution = None

    elif burn_model == "finite_burn":
        burn_result = apply_finite_burn_to_chaser(
            y=y_burn_minus,
            delta_v_km_s=delta_v,
            gamma_rad=gamma,
            thrust_N=thrust_N,
            isp_s=isp_s,
            initial_mass_kg=initial_mass_kg,
            rtol=rtol,
            atol=atol,
            max_step_burn_s=max_step_burn_s,
            earth=earth
        )
        y_burn_plus = burn_result["y_plus"]
        burn_duration_s = burn_result["burn_duration_s"]
        delivered_delta_v_km_s = burn_result["delivered_delta_v_km_s"]
        final_mass_kg = burn_result["final_mass_kg"]
        propellant_used_kg = burn_result["propellant_used_kg"]
        delta_v_vec_eci = burn_result["equivalent_delta_v_vec_eci_km_s"]
        initial_burn_direction_eci = burn_result["initial_burn_direction_eci"]
        burn_solution = burn_result["solution"]

    else:
        raise ValueError(f"Unsupported burn_model: {burn_model}")

    if post_burn_duration_s is None:
        post_burn_duration_s = 1.5 * defaults["transfer_time_s"]

    sol_transfer = solve_ivp(
        fun=lambda t, y: target_chaser_ode_j2(t, y, earth),
        t_span=(0.0, post_burn_duration_s),
        y0=y_burn_plus,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        dense_output=True,
        max_step=max_step_transfer_s
    )

    if not sol_transfer.success:
        raise RuntimeError(f"Post-burn propagation failed: {sol_transfer.message}")

    closest = find_closest_approach(
        sol=sol_transfer,
        t_start=0.0,
        t_end=post_burn_duration_s,
        sample_count=closest_sample_count,
        desired_rel_lvlh=desired_rel_lvlh
    )

    result = {
        "input_used": {
            "burn_model": burn_model,
            "phase_angle_rad": phase_angle,
            "phase_angle_deg": np.degrees(phase_angle),
            "delta_v_km_s": delta_v,
            "delta_v_m_s": delta_v * 1000.0,
            "gamma_rad": gamma,
            "gamma_deg": np.degrees(gamma),
            "thrust_N": thrust_N,
            "isp_s": isp_s,
            "initial_mass_kg": initial_mass_kg,
            "initial_phase_angle_rad": initial_phase_angle,
            "initial_phase_angle_deg": np.degrees(initial_phase_angle)
        },
        "non_j2_hohmann_reference": defaults,
        "burn": {
            "burn_model": burn_model,
            "t_burn_since_sim_start_s": t_burn,
            "t_burn_since_sim_start_min": t_burn / 60.0,
            "t_burn_since_sim_start_hr": t_burn / 3600.0,
            "burn_duration_s": burn_duration_s,
            "burn_end_since_sim_start_s": t_burn + burn_duration_s,
            "delivered_delta_v_km_s": delivered_delta_v_km_s,
            "delivered_delta_v_m_s": delivered_delta_v_km_s * 1000.0,
            "initial_mass_kg": initial_mass_kg,
            "final_mass_kg": final_mass_kg,
            "propellant_used_kg": propellant_used_kg,
            "thrust_N": thrust_N,
            "isp_s": isp_s,
            "target_state_eci_before_burn": y_burn_minus[0:6],
            "chaser_state_eci_before_burn": y_burn_minus[6:12],
            "target_state_eci_after_burn": y_burn_plus[0:6],
            "chaser_state_eci_after_burn": y_burn_plus[6:12],
            "initial_burn_direction_eci": initial_burn_direction_eci,
            "delta_v_vec_eci_km_s": delta_v_vec_eci,
            "finite_burn_solution": burn_solution
        },
        "closest_approach": {
            **closest,
            "t_closest_since_sim_start_s": t_burn + burn_duration_s + closest["t_closest_after_burn_s"],
            "t_closest_since_sim_start_min": (t_burn + burn_duration_s + closest["t_closest_after_burn_s"]) / 60.0,
            "t_closest_since_sim_start_hr": (t_burn + burn_duration_s + closest["t_closest_after_burn_s"]) / 3600.0
        }
    }

    if verbose:
        print("========== J2 Polar Hohmann-like Rendezvous Propagation ==========")
        print(f"Target altitude              : {h_target_km:.3f} km")
        print(f"Chaser altitude              : {h_chaser_km:.3f} km")
        print()
        print("[Non-J2 Hohmann reference]")
        print(f"Reference phase angle         : {defaults['phase_angle_rad']:.10f} rad")
        print(f"                              : {defaults['phase_angle_deg']:.6f} deg")
        print(f"Reference delta-V 1           : {defaults['delta_v_1_km_s']:.10f} km/s")
        print(f"                              : {defaults['delta_v_1_km_s'] * 1000.0:.6f} m/s")
        print(f"Reference transfer time       : {defaults['transfer_time_s']:.6f} s")
        print(f"                              : {defaults['transfer_time_s'] / 60.0:.6f} min")
        print()
        print("[Input actually used]")
        print(f"Burn model                    : {burn_model}")
        print(f"Phase angle                   : {phase_angle:.10f} rad")
        print(f"                              : {np.degrees(phase_angle):.6f} deg")
        print(f"Delta-V                       : {delta_v:.10f} km/s")
        print(f"                              : {delta_v * 1000.0:.6f} m/s")
        print(f"Gamma                         : {gamma:.10f} rad")
        print(f"                              : {np.degrees(gamma):.6f} deg")
        print()
        print("[Burn]")
        print(f"Burn time since sim start     : {t_burn:.6f} s")
        print(f"                              : {t_burn / 60.0:.6f} min")
        print(f"                              : {t_burn / 3600.0:.6f} hr")
        print(f"Burn duration                 : {burn_duration_s:.6f} s")
        print(f"Delivered delta-V budget      : {delivered_delta_v_km_s:.10f} km/s")
        print(f"                              : {delivered_delta_v_km_s * 1000.0:.6f} m/s")
        print(f"Propellant used               : {propellant_used_kg:.9f} kg")
        print(f"Equivalent delta-V vector ECI : {delta_v_vec_eci} km/s")
        print()
        print("[Closest approach]")
        print(f"Closest after burn            : {closest['t_closest_after_burn_s']:.6f} s")
        print(f"Closest since sim start       : {t_burn + burn_duration_s + closest['t_closest_after_burn_s']:.6f} s")
        print(f"                              : {(t_burn + burn_duration_s + closest['t_closest_after_burn_s']) / 60.0:.6f} min")
        print(f"                              : {(t_burn + burn_duration_s + closest['t_closest_after_burn_s']) / 3600.0:.6f} hr")
        print(f"Minimum distance              : {closest['min_distance_km']:.9f} km")
        print(f"Relative speed at closest     : {closest['relative_speed_km_s']:.10f} km/s")
        print(f"                              : {closest['relative_speed_km_s'] * 1000.0:.6f} m/s")
        print("==================================================================")

    return result


# ============================================================
# Example run
# ============================================================

if __name__ == "__main__":
    result = run_j2_polar_hohmann_rendezvous(
        phase_angle=None,     # default: non-J2 Hohmann phase angle
        delta_v=None,         # default: non-J2 Hohmann first impulse
        gamma=0.0,            # pure tangential burn
        burn_model="impulsive",
        thrust_N=1.0,
        isp_s=200.0,
        initial_mass_kg=2000.0,
        angle_unit="rad",
        initial_phase_angle=np.pi/2,
        initial_chaser_angle=0.0, # chaser inserting angle, clockwise direction from the north pole
        desired_rel_lvlh=[0,-5,0],
        verbose=True
    )
