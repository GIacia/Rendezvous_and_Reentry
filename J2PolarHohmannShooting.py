import numpy as np
import matplotlib.pyplot as plt
DESIRED_REL_LVLH = np.array([0.0, -5.0, 0.0])

from J2PolarHohmann import (
    run_j2_polar_hohmann_rendezvous,
    non_j2_hohmann_defaults
)


# ============================================================
# Parameter conversion helpers
# ============================================================

def get_default_parameter_vector():
    """
    Human-friendly optimization variables:

    p = [
        phase_angle_deg,
        delta_v_m_s,
        gamma_deg
    ]
    """
    defaults = non_j2_hohmann_defaults()

    phase_angle_deg = defaults["phase_angle_deg"]
    delta_v_m_s = defaults["delta_v_1_km_s"] * 1000.0
    gamma_deg = 0.0

    return np.array([phase_angle_deg, delta_v_m_s, gamma_deg], dtype=float)


def wrap_phase_deg(phase_deg):
    return phase_deg % 360.0


def sanitize_parameters(p):
    """
    Keeps parameters inside reasonable physical/numerical ranges.

    p[0] : phase angle [deg]
    p[1] : delta-V [m/s]
    p[2] : gamma [deg]
    """
    p = np.array(p, dtype=float).copy()

    p[0] = wrap_phase_deg(p[0])

    # Avoid negative burn magnitude.
    p[1] = max(p[1], 0.0)

    # For this first shooting attempt, keep gamma moderate.
    # You can widen this later if needed.
    p[2] = np.clip(p[2], -45.0, 45.0)

    return p


# ============================================================
# One simulation evaluation
# ============================================================

def evaluate_rendezvous_error(p, verbose=False):
    """
    Runs the J2 propagator once.

    Returns
    -------
    residual : np.ndarray
        Target-LVLH position error [radial, along-track] in km.

    distance : float
        Minimum distance to desired target-centered LVLH point in km.

    result : dict
        Full propagation result.
    """
    p = sanitize_parameters(p)

    phase_angle_deg = p[0]
    delta_v_km_s = p[1] / 1000.0
    gamma_deg = p[2]

    result = run_j2_polar_hohmann_rendezvous(
        phase_angle=phase_angle_deg,
        delta_v=delta_v_km_s,
        gamma=gamma_deg,
        angle_unit="deg",
        initial_phase_angle=90.0,
        initial_chaser_angle=0.0,
        desired_rel_lvlh=DESIRED_REL_LVLH,
        verbose=False
    )

    error_lvlh = result["closest_approach"]["position_error_lvlh_km"]

    # In the same orbital plane, the meaningful errors are:
    # x_LVLH = radial error
    # y_LVLH = along-track error
    residual = np.array([error_lvlh[0], error_lvlh[1]], dtype=float)

    distance = result["closest_approach"]["min_distance_km"]

    if verbose:
        print("p =", p)
        print("desired_rel_lvlh km =", DESIRED_REL_LVLH)
        print("error_lvlh km =", error_lvlh)
        print("residual [radial, along-track] km =", residual)
        print("distance km =", distance)

    return residual, distance, result

# ============================================================
# Finite-difference Jacobian
# ============================================================

def finite_difference_jacobian(
    p,
    finite_steps=np.array([0.01, 0.01, 0.01])
):
    """
    Numerically computes Jacobian:

        J = d residual / d parameters

    Parameters are:
        phase angle [deg]
        delta-V [m/s]
        gamma [deg]

    finite_steps:
        [deg, m/s, deg]
    """
    p = sanitize_parameters(p)

    residual_0, distance_0, result_0 = evaluate_rendezvous_error(p)

    m = len(residual_0)
    n = len(p)

    J = np.zeros((m, n))

    for j in range(n):
        h = finite_steps[j]

        p_plus = p.copy()
        p_minus = p.copy()

        p_plus[j] += h
        p_minus[j] -= h

        p_plus = sanitize_parameters(p_plus)
        p_minus = sanitize_parameters(p_minus)

        residual_plus, _, _ = evaluate_rendezvous_error(p_plus)
        residual_minus, _, _ = evaluate_rendezvous_error(p_minus)

        J[:, j] = (residual_plus - residual_minus) / (2.0 * h)

    return J, residual_0, distance_0, result_0


# ============================================================
# Damped Gauss-Newton optimizer
# ============================================================

def optimize_j2_hohmann_rendezvous(
    p0=None,
    max_iter=30,
    tolerance_km=1e-3,
    damping_initial=1e-2,
    finite_steps=np.array([0.01, 0.01, 0.01]),
    max_update=np.array([3.0, 20.0, 3.0]),
    verbose=True
):
    """
    Shooting method optimizer for:

        phase angle,
        delta-V,
        gamma

    using damped Gauss-Newton.

    Parameters
    ----------
    p0:
        Initial parameter vector:
            [phase_angle_deg, delta_v_m_s, gamma_deg]

        If None, uses non-J2 Hohmann default.

    tolerance_km:
        Stop when minimum distance becomes smaller than this.

    max_update:
        Maximum update per iteration:
            [deg, m/s, deg]

    Returns
    -------
    final_result : dict
        Includes optimized parameters, history, and final propagation result.
    """
    if p0 is None:
        p = get_default_parameter_vector()
    else:
        p = np.array(p0, dtype=float)

    p = sanitize_parameters(p)

    damping = damping_initial

    history = {
        "step": [],
        "phase_angle_deg": [],
        "delta_v_m_s": [],
        "gamma_deg": [],
        "distance_km": [],
        "residual_x_km": [],
        "residual_z_km": [],
        "closest_time_since_start_s": [],
        "closest_time_after_burn_s": []
    }

    best_p = p.copy()
    best_distance = np.inf
    best_result = None

    for k in range(max_iter):
        J, residual, distance, result = finite_difference_jacobian(
            p,
            finite_steps=finite_steps
        )

        if distance < best_distance:
            best_distance = distance
            best_p = p.copy()
            best_result = result

        history["step"].append(k)
        history["phase_angle_deg"].append(p[0])
        history["delta_v_m_s"].append(p[1])
        history["gamma_deg"].append(p[2])
        history["distance_km"].append(distance)
        history["residual_x_km"].append(residual[0])
        history["residual_z_km"].append(residual[1])
        history["closest_time_since_start_s"].append(
            result["closest_approach"]["t_closest_since_sim_start_s"]
        )
        history["closest_time_after_burn_s"].append(
            result["closest_approach"]["t_closest_after_burn_s"]
        )

        if verbose:
            print(f"\nIteration {k}")
            print("----------------------------------------")
            print(f"phase angle : {p[0]:.10f} deg")
            print(f"delta-V     : {p[1]:.10f} m/s")
            print(f"gamma       : {p[2]:.10f} deg")
            print(f"distance    : {distance:.12f} km")
            print(f"residual    : [{residual[0]:.12f}, {residual[1]:.12f}] km")
            print(f"damping     : {damping:.3e}")

        if distance < tolerance_km:
            if verbose:
                print("\nConverged.")
            break

        # Damped Gauss-Newton step:
        #
        #   (J^T J + lambda I) dp = -J^T r
        #
        A = J.T @ J + damping * np.eye(3)
        b = -J.T @ residual

        try:
            dp = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            dp = np.linalg.lstsq(A, b, rcond=None)[0]

        # Prevent one iteration from jumping too far.
        dp = np.clip(dp, -max_update, max_update)

        # Line search.
        # Accept only if distance decreases.
        accepted = False
        alpha = 1.0

        for _ in range(10):
            p_trial = sanitize_parameters(p + alpha * dp)

            residual_trial, distance_trial, result_trial = evaluate_rendezvous_error(
                p_trial
            )

            if distance_trial < distance:
                p = p_trial
                damping = max(damping * 0.5, 1e-8)
                accepted = True

                if verbose:
                    print(f"accepted alpha : {alpha:.5f}")
                    print(f"new distance   : {distance_trial:.12f} km")

                break

            alpha *= 0.5

        if not accepted:
            damping *= 10.0

            if verbose:
                print("step rejected; increasing damping.")

    final_residual, final_distance, final_result = evaluate_rendezvous_error(best_p)

    final_output = {
        "optimized_parameters": {
            "phase_angle_deg": best_p[0],
            "phase_angle_rad": np.radians(best_p[0]),
            "delta_v_m_s": best_p[1],
            "delta_v_km_s": best_p[1] / 1000.0,
            "gamma_deg": best_p[2],
            "gamma_rad": np.radians(best_p[2])
        },
        "best_distance_km": final_distance,
        "best_residual_xz_km": final_residual,
        "history": history,
        "final_propagation_result": final_result
    }

    return final_output


# ============================================================
# Plotting
# ============================================================

def plot_distance_history(history, log_scale=True):
    steps = history["step"]
    distances = history["distance_km"]

    plt.figure(figsize=(8, 5))
    plt.plot(steps, distances, marker="o")
    plt.xlabel("Step")
    plt.ylabel("Distance error [km]")
    plt.title("J2 Polar Hohmann Shooting: Distance Error History")
    plt.grid(True)

    if log_scale:
        plt.yscale("log")

    plt.tight_layout()
    plt.show()


def plot_parameter_history(history):
    steps = history["step"]

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["phase_angle_deg"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Phase angle [deg]")
    plt.title("Phase Angle History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["delta_v_m_s"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Delta-V [m/s]")
    plt.title("Delta-V History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["gamma_deg"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Gamma [deg]")
    plt.title("Gamma History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

from scipy.optimize import minimize


# ============================================================
# Extract parameter vector from previous optimizer output
# ============================================================

def extract_parameter_vector_from_output(output):
    """
    Converts output from optimize_j2_hohmann_rendezvous()
    into parameter vector:

        p = [phase_angle_deg, delta_v_m_s, gamma_deg]
    """
    params = output["optimized_parameters"]

    return np.array([
        params["phase_angle_deg"],
        params["delta_v_m_s"],
        params["gamma_deg"]
    ], dtype=float)


# ============================================================
# Constrained delta-V minimization
# ============================================================

def minimize_delta_v_on_zero_distance_manifold(
    p_start=None,
    objective_mode="first_burn",
    feasibility_tolerance_km=1e-3,
    residual_scale_km=1.0,
    delta_v_scale_m_s=100.0,
    stage1_max_iter=25,
    stage2_max_iter=80,
    stage1_verbose=True,
    stage2_verbose=True,
    bounds=None,
    ftol=1e-9
):
    """
    Minimize delta-V while keeping rendezvous position error zero.

    Variables
    ---------
    p = [
        phase_angle_deg,
        delta_v_m_s,
        gamma_deg
    ]

    Optimization problem
    --------------------
    first_burn mode:

        minimize     delta_v_1

        subject to   rel_pos_x = 0
                     rel_pos_z = 0

    two_impulse_total mode:

        minimize     delta_v_1 + terminal_relative_speed

        subject to   rel_pos_x = 0
                     rel_pos_z = 0

    Notes
    -----
    The constraint uses the relative position at closest approach.
    Since the orbit is a polar orbit in the fixed x-z ECI plane,
    the meaningful residual is [rel_x, rel_z].

    Returns
    -------
    output : dict
    """
    if objective_mode not in ["first_burn", "two_impulse_total"]:
        raise ValueError(
            "objective_mode must be either 'first_burn' or 'two_impulse_total'."
        )

    if bounds is None:
        bounds = [
            (0.0, 360.0),     # phase angle [deg]
            (0.0, 500.0),     # delta-V [m/s]
            (-60.0, 60.0)     # gamma [deg]
        ]

    # ------------------------------------------------------------
    # Stage 0: initial guess
    # ------------------------------------------------------------
    if p_start is None:
        p_start = get_default_parameter_vector()
    else:
        p_start = np.array(p_start, dtype=float)

    p_start = sanitize_parameters(p_start)

    residual_start, distance_start, result_start = evaluate_rendezvous_error(
        p_start
    )

    # ------------------------------------------------------------
    # Stage 1: make position feasible first
    # ------------------------------------------------------------
    if distance_start > feasibility_tolerance_km:
        if stage1_verbose:
            print("\n========== Stage 1: Feasibility shooting ==========")
            print("Initial distance is not small enough.")
            print(f"Initial distance: {distance_start:.12f} km")
            print("Running position-only shooting first.")
            print("===================================================\n")

        feasibility_output = optimize_j2_hohmann_rendezvous(
            p0=p_start,
            max_iter=stage1_max_iter,
            tolerance_km=feasibility_tolerance_km,
            verbose=stage1_verbose
        )

        p_feasible = extract_parameter_vector_from_output(feasibility_output)

    else:
        feasibility_output = None
        p_feasible = p_start.copy()

    p_feasible = sanitize_parameters(p_feasible)

    residual_feasible, distance_feasible, result_feasible = evaluate_rendezvous_error(
        p_feasible
    )

    if stage2_verbose:
        print("\n========== Stage 2: Delta-V minimization ==========")
        print("Starting constrained optimization from:")
        print(f"phase angle : {p_feasible[0]:.10f} deg")
        print(f"delta-V     : {p_feasible[1]:.10f} m/s")
        print(f"gamma       : {p_feasible[2]:.10f} deg")
        print(f"distance    : {distance_feasible:.12f} km")
        print("===================================================\n")

    # ------------------------------------------------------------
    # Evaluation cache
    # ------------------------------------------------------------
    cache = {}

    def cache_key(p):
        p_clean = sanitize_parameters(p)
        return tuple(np.round(p_clean, decimals=10))

    def eval_cached(p):
        key = cache_key(p)

        if key not in cache:
            p_clean = sanitize_parameters(p)

            try:
                residual, distance, result = evaluate_rendezvous_error(
                    p_clean
                )

                cache[key] = {
                    "p": p_clean,
                    "residual": residual,
                    "distance": distance,
                    "result": result,
                    "failed": False
                }

            except Exception as e:
                cache[key] = {
                    "p": p_clean,
                    "residual": np.array([1e6, 1e6], dtype=float),
                    "distance": 1e6,
                    "result": None,
                    "failed": True,
                    "error": e
                }

        return cache[key]

    # ------------------------------------------------------------
    # Objective function
    # ------------------------------------------------------------
    def objective(p):
        data = eval_cached(p)
        p_clean = data["p"]

        delta_v_1_m_s = p_clean[1]

        if objective_mode == "first_burn":
            return delta_v_1_m_s / delta_v_scale_m_s

        if objective_mode == "two_impulse_total":
            result = data["result"]

            if result is None:
                return 1e9

            terminal_relative_speed_m_s = (
                result["closest_approach"]["relative_speed_km_s"] * 1000.0
            )

            total_delta_v_m_s = delta_v_1_m_s + terminal_relative_speed_m_s

            return total_delta_v_m_s / delta_v_scale_m_s

    # ------------------------------------------------------------
    # Equality / Inequality constraint
    # ------------------------------------------------------------
    constraint_tolerance_km = 0.0001

    def equality_constraints(p):
        data = eval_cached(p)
        residual = data["residual"]

        return residual / residual_scale_km
    
    def inequality_constraints(p):
        """
        SLSQP inequality constraint:
            fun(p) >= 0

        Here:
            constraint_tolerance_km - distance >= 0

        means:
            distance <= constraint_tolerance_km
        """
        data = eval_cached(p)
        distance = data["distance"]

        return constraint_tolerance_km - distance

    constraints = [
        {
            "type": "eq",
            "fun": equality_constraints
        }
    ]

    # ------------------------------------------------------------
    # History recording
    # ------------------------------------------------------------
    history = {
        "step": [],
        "phase_angle_deg": [],
        "delta_v_m_s": [],
        "gamma_deg": [],
        "distance_km": [],
        "residual_x_km": [],
        "residual_z_km": [],
        "relative_speed_m_s": [],
        "objective_value": []
    }

    step_counter = {"k": 0}

    def record_history(p):
        data = eval_cached(p)
        p_clean = data["p"]
        residual = data["residual"]
        distance = data["distance"]
        result = data["result"]

        if result is None:
            relative_speed_m_s = np.nan
        else:
            relative_speed_m_s = (
                result["closest_approach"]["relative_speed_km_s"] * 1000.0
            )

        history["step"].append(step_counter["k"])
        history["phase_angle_deg"].append(p_clean[0])
        history["delta_v_m_s"].append(p_clean[1])
        history["gamma_deg"].append(p_clean[2])
        history["distance_km"].append(distance)
        history["residual_x_km"].append(residual[0])
        history["residual_z_km"].append(residual[1])
        history["relative_speed_m_s"].append(relative_speed_m_s)
        history["objective_value"].append(objective(p_clean))

        step_counter["k"] += 1

    def callback(p):
        record_history(p)

        if stage2_verbose:
            data = eval_cached(p)
            p_clean = data["p"]

            print(f"SQP step {step_counter['k'] - 1}")
            print("----------------------------------------")
            print(f"phase angle : {p_clean[0]:.10f} deg")
            print(f"delta-V     : {p_clean[1]:.10f} m/s")
            print(f"gamma       : {p_clean[2]:.10f} deg")
            print(f"distance    : {data['distance']:.12f} km")
            print(f"residual    : {data['residual']} km")
            print(f"objective   : {objective(p_clean):.12f}")
            print()

    # Record initial feasible point.
    record_history(p_feasible)

    # ------------------------------------------------------------
    # Stage 2: SLSQP constrained minimization
    # ------------------------------------------------------------
    solution = minimize(
        fun=objective,
        x0=p_feasible,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        callback=callback,
        options={
            "maxiter": stage2_max_iter,
            "ftol": ftol,
            "disp": stage2_verbose
        }
    )

    p_final = sanitize_parameters(solution.x)

    final_data = eval_cached(p_final)
    final_residual = final_data["residual"]
    final_distance = final_data["distance"]
    final_result = final_data["result"]

    if final_result is None:
        terminal_relative_speed_m_s = np.nan
        total_delta_v_m_s = np.nan
    else:
        terminal_relative_speed_m_s = (
            final_result["closest_approach"]["relative_speed_km_s"] * 1000.0
        )
        total_delta_v_m_s = p_final[1] + terminal_relative_speed_m_s

    output = {
        "success": solution.success,
        "message": solution.message,
        "objective_mode": objective_mode,
        "optimized_parameters": {
            "phase_angle_deg": p_final[0],
            "phase_angle_rad": np.radians(p_final[0]),
            "delta_v_1_m_s": p_final[1],
            "delta_v_1_km_s": p_final[1] / 1000.0,
            "gamma_deg": p_final[2],
            "gamma_rad": np.radians(p_final[2])
        },
        "final_distance_km": final_distance,
        "final_residual_xz_km": final_residual,
        "terminal_relative_speed_m_s": terminal_relative_speed_m_s,
        "total_two_impulse_delta_v_m_s": total_delta_v_m_s,
        "history": history,
        "feasibility_output": feasibility_output,
        "final_propagation_result": final_result,
        "raw_scipy_solution": solution
    }

    if stage2_verbose:
        print("\n========== Constrained Delta-V Result ==========")
        print(f"success                      : {solution.success}")
        print(f"message                      : {solution.message}")
        print(f"objective mode               : {objective_mode}")
        print()
        print(f"phase angle                  : {p_final[0]:.10f} deg")
        print(f"delta-V 1                    : {p_final[1]:.10f} m/s")
        print(f"gamma                        : {p_final[2]:.10f} deg")
        print()
        print(f"final distance               : {final_distance:.12f} km")
        print(f"final residual [x, z]        : {final_residual} km")
        print(f"terminal relative speed      : {terminal_relative_speed_m_s:.10f} m/s")
        print(f"two-impulse total delta-V    : {total_delta_v_m_s:.10f} m/s")
        print("================================================")

    return output


# ============================================================
# Plotting for constrained optimization
# ============================================================
def plot_constrained_delta_v_history(history, log_distance=True):
    steps = history["step"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(steps, history["distance_km"], marker="o")
    ax.set_xlabel("Step")
    ax.set_ylabel("Distance error [km]")
    ax.set_title("Distance Error")
    ax.grid(True)
    if log_distance:
        ax.set_yscale("log")

    ax = axes[0, 1]
    ax.plot(steps, history["delta_v_m_s"], marker="o")
    ax.set_xlabel("Step")
    ax.set_ylabel("Delta-V 1 [m/s]")
    ax.set_title("Delta-V History")
    ax.grid(True)

    ax = axes[1, 0]
    ax.plot(steps, history["gamma_deg"], marker="o")
    ax.set_xlabel("Step")
    ax.set_ylabel("Gamma [deg]")
    ax.set_title("Gamma History")
    ax.grid(True)

    ax = axes[1, 1]
    ax.plot(steps, history["relative_speed_m_s"], marker="o")
    ax.set_xlabel("Step")
    ax.set_ylabel("Relative speed [m/s]")
    ax.set_title("Terminal Relative Speed")
    ax.grid(True)

    fig.suptitle("Constrained J2 Polar Hohmann Optimization History", fontsize=14)
    plt.tight_layout()
    plt.show()
'''
def plot_constrained_delta_v_history(history, log_distance=True):
    steps = history["step"]

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["distance_km"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Distance error [km]")
    plt.title("Constrained Optimization: Distance Error History")
    plt.grid(True)

    if log_distance:
        plt.yscale("log")

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["delta_v_m_s"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Delta-V 1 [m/s]")
    plt.title("Constrained Optimization: Delta-V History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["gamma_deg"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Gamma [deg]")
    plt.title("Constrained Optimization: Gamma History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(8, 5))
    plt.plot(steps, history["relative_speed_m_s"], marker="o")
    plt.xlabel("Step")
    plt.ylabel("Relative speed at closest approach [m/s]")
    plt.title("Terminal Relative Speed History")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
'''


# ============================================================
# Example run
# ============================================================

if __name__ == "__main__":
    output = minimize_delta_v_on_zero_distance_manifold(
        p_start=None,
        objective_mode="two_impulse_total",
        stage1_verbose=True,
        stage2_verbose=True
    )

    plot_constrained_delta_v_history(
        output["history"],
        log_distance=True
    )