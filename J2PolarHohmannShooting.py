import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import matplotlib.pyplot as plt
DESIRED_REL_LVLH = np.array([0.0, -5.0, 0.0])

DEFAULT_MANEUVER_CONFIG = {
    "burn_model": "impulsive",
    "thrust_N": 300.0,
    "isp_s": 200.0,
    "initial_mass_kg": 2000.0,
    "max_step_burn_s": 0.5
}

DEFAULT_SCENARIO_CONFIG = {
    "h_chaser_km": 300.0,
    "h_target_km": 500.0,
    "initial_phase_angle": 90.0,
    "initial_chaser_angle": 0.0,
    "desired_rel_lvlh": DESIRED_REL_LVLH,
    "max_wait_time_s": None,
    "post_burn_duration_s": None,
    "closest_sample_count": 1200,
}

DEFAULT_ENVIRONMENT_CONFIG = {
    "atmospheric_drag": {
        "enabled": False,
        "model": "ISA76",
        "use_matlab_atmosisa": True,
        "co_rotate_atmosphere": True,
        "earth_rotation_rad_s": 7.2921159e-5,
        "chaser_cd": 2.2,
        "chaser_area_m2": 4.0,
        "target_cd": 2.2,
        "target_area_m2": 4.0,
        "target_mass_kg": 2000.0
    }
}

from J2PolarHohmann import (
    run_j2_polar_hohmann_rendezvous,
    non_j2_hohmann_defaults
)


def make_maneuver_config(
    burn_model="impulsive",
    thrust_N=300.0,
    isp_s=200.0,
    initial_mass_kg=2000.0,
    max_step_burn_s=0.5
):
    return {
        "burn_model": burn_model,
        "thrust_N": thrust_N,
        "isp_s": isp_s,
        "initial_mass_kg": initial_mass_kg,
        "max_step_burn_s": max_step_burn_s
    }


def resolve_maneuver_config(maneuver_config=None):
    config = DEFAULT_MANEUVER_CONFIG.copy()
    if maneuver_config is not None:
        config.update(maneuver_config)
    return config


def make_scenario_config(
    h_chaser_km=300.0,
    h_target_km=500.0,
    initial_phase_angle=90.0,
    initial_chaser_angle=0.0,
    desired_rel_lvlh=DESIRED_REL_LVLH,
    max_wait_time_s=None,
    post_burn_duration_s=None,
    closest_sample_count=1200
):
    return {
        "h_chaser_km": h_chaser_km,
        "h_target_km": h_target_km,
        "initial_phase_angle": initial_phase_angle,
        "initial_chaser_angle": initial_chaser_angle,
        "desired_rel_lvlh": np.array(desired_rel_lvlh, dtype=float),
        "max_wait_time_s": max_wait_time_s,
        "post_burn_duration_s": post_burn_duration_s,
        "closest_sample_count": closest_sample_count
    }


def resolve_scenario_config(scenario_config=None):
    config = DEFAULT_SCENARIO_CONFIG.copy()
    if scenario_config is not None:
        config.update(scenario_config)
    config["desired_rel_lvlh"] = np.array(config["desired_rel_lvlh"], dtype=float)
    return config


def make_environment_config(
    atmospheric_drag="off",
    chaser_cd=2.2,
    chaser_area_m2=4.0,
    target_cd=2.2,
    target_area_m2=4.0,
    target_mass_kg=2000.0,
    co_rotate_atmosphere=True,
    use_matlab_atmosisa=True
):
    model = str(atmospheric_drag).strip().lower()
    enabled = model not in {"off", "none", "disabled", "0"}
    if not enabled:
        model = "ISA76"

    return {
        "atmospheric_drag": {
            "enabled": enabled,
            "model": model.upper(),
            "use_matlab_atmosisa": bool(use_matlab_atmosisa),
            "co_rotate_atmosphere": bool(co_rotate_atmosphere),
            "earth_rotation_rad_s": 7.2921159e-5,
            "chaser_cd": chaser_cd,
            "chaser_area_m2": chaser_area_m2,
            "target_cd": target_cd,
            "target_area_m2": target_area_m2,
            "target_mass_kg": target_mass_kg
        }
    }


def resolve_environment_config(environment_config=None):
    config = {
        "atmospheric_drag": DEFAULT_ENVIRONMENT_CONFIG["atmospheric_drag"].copy()
    }
    if environment_config is None:
        return config

    if "atmospheric_drag" in environment_config:
        source_drag = environment_config["atmospheric_drag"]
        explicit_enabled = "enabled" in source_drag
        config["atmospheric_drag"].update(source_drag)
    else:
        explicit_enabled = "enabled" in environment_config
        config["atmospheric_drag"].update(environment_config)

    model = str(config["atmospheric_drag"].get("model", "ISA76")).upper()
    if not explicit_enabled and model not in {"OFF", "NONE", "DISABLED", "0"}:
        config["atmospheric_drag"]["enabled"] = True

    enabled = config["atmospheric_drag"].get("enabled", False)
    if isinstance(enabled, str):
        enabled = enabled.strip().lower() in {"true", "1", "yes", "on"}
    config["atmospheric_drag"]["enabled"] = bool(enabled)
    config["atmospheric_drag"]["model"] = model
    return config


def to_jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def canonical_json_text(value):
    return json.dumps(to_jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def short_hash(value, length=12):
    return hashlib.sha256(canonical_json_text(value).encode("utf-8")).hexdigest()[:length]


def config_settings_signature(config):
    return {
        "schema_version": config.get("schema_version"),
        "source": config.get("source"),
        "scenario": config.get("scenario", {}),
        "phase1_mode": (config.get("phase1", {}) or {}).get("mode"),
        "burn_model": (config.get("phase1", {}) or {}).get("burn_model"),
        "desired_rel_lvlh_m": (config.get("phase1", {}) or {}).get("desired_rel_lvlh_m"),
        "maneuver": config.get("maneuver", {}),
        "environment": config.get("environment", {}),
    }


def config_result_signature(config):
    return {
        "settings": config_settings_signature(config),
        "phase1_solution": config.get("phase1", {}),
        "optimizer": config.get("optimizer", {}),
    }


def safe_slug(value):
    text = str(value).strip().lower()
    chars = []
    for ch in text:
        if ch.isalnum():
            chars.append(ch)
        elif ch in {"_", "-", "."}:
            chars.append(ch)
        else:
            chars.append("-")
    slug = "".join(chars).strip("-")
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug or "case"


def archive_labels(config):
    scenario = config.get("scenario", {}) or {}
    phase1 = config.get("phase1", {}) or {}
    drag = ((config.get("environment", {}) or {}).get("atmospheric_drag", {}) or {})

    burn = safe_slug(phase1.get("burn_model", "unknown"))
    drag_enabled = bool(drag.get("enabled", False))
    drag_model = safe_slug(drag.get("model", "off")) if drag_enabled else "drag-off"
    h_chaser = scenario.get("h_chaser_km", "x")
    h_target = scenario.get("h_target_km", "x")
    phase0 = scenario.get("initial_phase_angle_deg", "x")

    scenario_label = f"h{safe_slug(h_chaser)}-{safe_slug(h_target)}_phase{safe_slug(phase0)}"
    return burn, drag_model, scenario_label


def attach_archive_metadata(config, created_at=None):
    config = to_jsonable(config)
    created_at = created_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    settings_hash = short_hash(config_settings_signature(config), 16)
    result_hash = short_hash(config_result_signature(config), 16)
    burn, drag_model, scenario_label = archive_labels(config)
    timestamp = created_at.replace("-", "").replace(":", "").replace("T", "_").replace("Z", "Z")
    case_id = f"{burn}_{drag_model}_{scenario_label}_{settings_hash[:8]}_{result_hash[:8]}_{timestamp}"

    config["archive"] = {
        "schema_version": 1,
        "created_at": created_at,
        "case_id": case_id,
        "settings_hash": settings_hash,
        "result_hash": result_hash,
        "settings_signature": config_settings_signature(config)
    }
    return config


def relative_posix_path(path, base_dir):
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def build_index_entry(config, archive_path, config_dir, latest_path=None):
    archive = config["archive"]
    phase1 = config.get("phase1", {}) or {}
    scenario = config.get("scenario", {}) or {}
    drag = ((config.get("environment", {}) or {}).get("atmospheric_drag", {}) or {})
    optimizer = config.get("optimizer", {}) or {}

    entry = {
        "case_id": archive["case_id"],
        "settings_hash": archive["settings_hash"],
        "result_hash": archive["result_hash"],
        "created_at": archive["created_at"],
        "path": relative_posix_path(archive_path, config_dir),
        "phase1": {
            "mode": phase1.get("mode"),
            "burn_model": phase1.get("burn_model"),
            "phase_angle_deg": phase1.get("phase_angle_deg"),
            "delta_v_m_s": phase1.get("delta_v_m_s"),
            "gamma_deg": phase1.get("gamma_deg"),
        },
        "scenario": scenario,
        "environment": {
            "atmospheric_drag_enabled": bool(drag.get("enabled", False)),
            "atmospheric_drag_model": drag.get("model", "ISA76"),
        },
        "optimizer": {
            "success": optimizer.get("success"),
            "objective_mode": optimizer.get("objective_mode"),
            "final_distance_km": optimizer.get("final_distance_km"),
            "total_two_impulse_delta_v_m_s": optimizer.get("total_two_impulse_delta_v_m_s"),
            "burn_duration_s": optimizer.get("burn_duration_s"),
        }
    }
    if latest_path is not None:
        entry["latest_alias_path"] = relative_posix_path(latest_path, config_dir)
    return entry


def update_solution_index(index_path, config, archive_path, latest_path=None):
    config_dir = index_path.parent
    if index_path.exists():
        try:
            index = json.loads(index_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            index = {}
    else:
        index = {}

    cases = index.get("cases", [])
    if not isinstance(cases, list):
        cases = []

    entry = build_index_entry(config, archive_path, config_dir, latest_path=latest_path)
    cases = [case for case in cases if case.get("case_id") != entry["case_id"]]
    cases.append(entry)
    cases.sort(key=lambda case: str(case.get("created_at", "")))

    index = {
        "schema_version": 1,
        "source": "J2PolarHohmannShooting.py",
        "updated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "latest_case_id": entry["case_id"],
        "cases": cases
    }
    index_path.write_text(json.dumps(index, indent=2, sort_keys=True), encoding="utf-8")


def matlab_burn_model_name(burn_model):
    model = str(burn_model).strip().lower()
    if model == "continuous":
        raise ValueError(
            "burn_model='continuous' is not supported. "
            "Use 'finite_burn' for a short finite-duration execution of an impulsive delta-V."
        )
    if model in {"finite", "finite_burn", "finite_impulse"}:
        return "FINITE_BURN"
    if model in {"impulsive", "instant", "instantaneous", "custom_impulse"}:
        return "IMPULSIVE"
    raise ValueError("burn_model must be 'impulsive' or 'finite_burn'.")


def get_optimized_delta_v_m_s(params):
    if "delta_v_1_m_s" in params:
        return params["delta_v_1_m_s"]
    return params["delta_v_m_s"]


def build_matlab_mission_config(output, phase1_mode="CUSTOM_IMPULSE"):
    scenario = resolve_scenario_config(output.get("scenario_config"))
    maneuver = resolve_maneuver_config(output.get("maneuver_config"))
    environment = resolve_environment_config(output.get("environment_config"))
    params = output["optimized_parameters"]
    final_result = output.get("final_propagation_result") or {}
    burn = final_result.get("burn", {})

    return {
        "schema_version": 1,
        "source": "J2PolarHohmannShooting.py",
        "scenario": {
            "h_chaser_km": scenario["h_chaser_km"],
            "h_target_km": scenario["h_target_km"],
            "initial_phase_angle_deg": scenario["initial_phase_angle"],
            "initial_chaser_angle_deg": scenario["initial_chaser_angle"]
        },
        "phase1": {
            "mode": phase1_mode,
            "burn_model": matlab_burn_model_name(maneuver["burn_model"]),
            "phase_angle_deg": params["phase_angle_deg"],
            "delta_v_m_s": get_optimized_delta_v_m_s(params),
            "gamma_deg": params["gamma_deg"],
            "desired_rel_lvlh_m": (scenario["desired_rel_lvlh"] * 1000.0).tolist()
        },
        "maneuver": {
            "finite_burn_thrust_N": maneuver["thrust_N"],
            "finite_burn_isp_s": maneuver["isp_s"],
            "finite_burn_dt_s": maneuver["max_step_burn_s"],
            "initial_mass_kg": maneuver["initial_mass_kg"]
        },
        "environment": environment,
        "optimizer": {
            "success": bool(output.get("success", True)),
            "message": str(output.get("message", "")),
            "objective_mode": output.get("objective_mode", ""),
            "final_distance_km": output.get("final_distance_km", output.get("best_distance_km")),
            "terminal_relative_speed_m_s": output.get("terminal_relative_speed_m_s"),
            "total_two_impulse_delta_v_m_s": output.get("total_two_impulse_delta_v_m_s"),
            "burn_duration_s": burn.get("burn_duration_s"),
            "delivered_delta_v_m_s": burn.get("delivered_delta_v_m_s")
        }
    }


def archive_matlab_mission_config(config, latest_path=None, archive_dir=None, index_path=None, created_at=None):
    latest_path = Path(latest_path) if latest_path is not None else Path("configs/latest_python_solution.json")
    config_dir = latest_path.parent
    archive_dir = Path(archive_dir) if archive_dir is not None else config_dir / "python_runs"
    index_path = Path(index_path) if index_path is not None else config_dir / "python_solution_index.json"

    archive_dir.mkdir(parents=True, exist_ok=True)
    index_path.parent.mkdir(parents=True, exist_ok=True)

    config = attach_archive_metadata(config, created_at=created_at)
    archive_path = archive_dir / f"{config['archive']['case_id']}.json"
    archive_path.write_text(
        json.dumps(to_jsonable(config), indent=2, sort_keys=True),
        encoding="utf-8"
    )
    update_solution_index(index_path, config, archive_path, latest_path=latest_path)
    return config, archive_path, index_path


def save_matlab_mission_config(output, path="configs/latest_python_solution.json", phase1_mode="CUSTOM_IMPULSE"):
    config = build_matlab_mission_config(output, phase1_mode=phase1_mode)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    config, archive_path, index_path = archive_matlab_mission_config(config, latest_path=path)
    path.write_text(
        json.dumps(to_jsonable(config), indent=2, sort_keys=True),
        encoding="utf-8"
    )
    print(f"Archived MATLAB mission config: {archive_path}")
    print(f"Updated MATLAB config index: {index_path}")
    return path


# ============================================================
# Parameter conversion helpers
# ============================================================

def get_default_parameter_vector(scenario_config=None):
    """
    Human-friendly optimization variables:

    p = [
        phase_angle_deg,
        delta_v_m_s,
        gamma_deg
    ]
    """
    scenario_config = resolve_scenario_config(scenario_config)
    defaults = non_j2_hohmann_defaults(
        h_chaser_km=scenario_config["h_chaser_km"],
        h_target_km=scenario_config["h_target_km"]
    )

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

def evaluate_rendezvous_error(
    p,
    verbose=False,
    maneuver_config=None,
    scenario_config=None,
    environment_config=None
):
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
    maneuver_config = resolve_maneuver_config(maneuver_config)
    scenario_config = resolve_scenario_config(scenario_config)
    environment_config = resolve_environment_config(environment_config)

    phase_angle_deg = p[0]
    delta_v_km_s = p[1] / 1000.0
    gamma_deg = p[2]

    result = run_j2_polar_hohmann_rendezvous(
        phase_angle=phase_angle_deg,
        delta_v=delta_v_km_s,
        gamma=gamma_deg,
        angle_unit="deg",
        h_target_km=scenario_config["h_target_km"],
        h_chaser_km=scenario_config["h_chaser_km"],
        initial_phase_angle=scenario_config["initial_phase_angle"],
        initial_chaser_angle=scenario_config["initial_chaser_angle"],
        max_wait_time_s=scenario_config["max_wait_time_s"],
        post_burn_duration_s=scenario_config["post_burn_duration_s"],
        closest_sample_count=scenario_config["closest_sample_count"],
        desired_rel_lvlh=scenario_config["desired_rel_lvlh"],
        environment_config=environment_config,
        verbose=False,
        **maneuver_config
    )

    error_lvlh = result["closest_approach"]["position_error_lvlh_km"]

    # In the same orbital plane, the meaningful errors are:
    # x_LVLH = radial error
    # y_LVLH = along-track error
    residual = np.array([error_lvlh[0], error_lvlh[1]], dtype=float)

    distance = result["closest_approach"]["min_distance_km"]

    if verbose:
        print("p =", p)
        print("scenario =", scenario_config)
        print("error_lvlh km =", error_lvlh)
        print("residual [radial, along-track] km =", residual)
        print("distance km =", distance)

    return residual, distance, result

# ============================================================
# Finite-difference Jacobian
# ============================================================

def finite_difference_jacobian(
    p,
    finite_steps=np.array([0.01, 0.01, 0.01]),
    maneuver_config=None,
    scenario_config=None,
    environment_config=None
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

    residual_0, distance_0, result_0 = evaluate_rendezvous_error(
        p,
        maneuver_config=maneuver_config,
        scenario_config=scenario_config,
        environment_config=environment_config
    )

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

        residual_plus, _, _ = evaluate_rendezvous_error(
            p_plus,
            maneuver_config=maneuver_config,
            scenario_config=scenario_config,
            environment_config=environment_config
        )
        residual_minus, _, _ = evaluate_rendezvous_error(
            p_minus,
            maneuver_config=maneuver_config,
            scenario_config=scenario_config,
            environment_config=environment_config
        )

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
    maneuver_config=None,
    scenario_config=None,
    environment_config=None,
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
    maneuver_config = resolve_maneuver_config(maneuver_config)
    scenario_config = resolve_scenario_config(scenario_config)
    environment_config = resolve_environment_config(environment_config)

    if p0 is None:
        p = get_default_parameter_vector(scenario_config=scenario_config)
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
            finite_steps=finite_steps,
            maneuver_config=maneuver_config,
            scenario_config=scenario_config,
            environment_config=environment_config
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
                p_trial,
                maneuver_config=maneuver_config,
                scenario_config=scenario_config,
                environment_config=environment_config
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

    final_residual, final_distance, final_result = evaluate_rendezvous_error(
        best_p,
        maneuver_config=maneuver_config,
        scenario_config=scenario_config,
        environment_config=environment_config
    )

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
        "maneuver_config": maneuver_config,
        "scenario_config": scenario_config,
        "environment_config": environment_config,
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
    maneuver_config=None,
    scenario_config=None,
    environment_config=None,
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

    maneuver_config = resolve_maneuver_config(maneuver_config)
    scenario_config = resolve_scenario_config(scenario_config)
    environment_config = resolve_environment_config(environment_config)

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
        p_start = get_default_parameter_vector(scenario_config=scenario_config)
    else:
        p_start = np.array(p_start, dtype=float)

    p_start = sanitize_parameters(p_start)

    residual_start, distance_start, result_start = evaluate_rendezvous_error(
        p_start,
        maneuver_config=maneuver_config,
        scenario_config=scenario_config,
        environment_config=environment_config
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
            maneuver_config=maneuver_config,
            scenario_config=scenario_config,
            environment_config=environment_config,
            verbose=stage1_verbose
        )

        p_feasible = extract_parameter_vector_from_output(feasibility_output)

    else:
        feasibility_output = None
        p_feasible = p_start.copy()

    p_feasible = sanitize_parameters(p_feasible)

    residual_feasible, distance_feasible, result_feasible = evaluate_rendezvous_error(
        p_feasible,
        maneuver_config=maneuver_config,
        scenario_config=scenario_config,
        environment_config=environment_config
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
                    p_clean,
                    maneuver_config=maneuver_config,
                    scenario_config=scenario_config,
                    environment_config=environment_config
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
        "maneuver_config": maneuver_config,
        "scenario_config": scenario_config,
        "environment_config": environment_config,
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
    parser = argparse.ArgumentParser(
        description="Optimize J2 polar rendezvous parameters and export a MATLAB mission JSON."
    )
    parser.add_argument(
        "--matlab-config-out",
        default="configs/latest_python_solution.json",
        help="Path for the MATLAB-readable JSON config."
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip plotting after optimization."
    )
    parser.add_argument(
        "--burn-model",
        default="impulsive",
        choices=["impulsive", "finite_burn"],
        help="Maneuver execution model used during Python optimization."
    )
    parser.add_argument(
        "--atmospheric-drag",
        default="off",
        choices=["off", "isa76"],
        help="Atmospheric drag model used during Python optimization and exported to MATLAB."
    )
    parser.add_argument(
        "--chaser-cd",
        type=float,
        default=2.2,
        help="Chaser drag coefficient."
    )
    parser.add_argument(
        "--chaser-area-m2",
        type=float,
        default=4.0,
        help="Chaser reference area for drag."
    )
    parser.add_argument(
        "--target-cd",
        type=float,
        default=2.2,
        help="Target drag coefficient."
    )
    parser.add_argument(
        "--target-area-m2",
        type=float,
        default=4.0,
        help="Target reference area for drag."
    )
    args = parser.parse_args()

    maneuver_config = make_maneuver_config(
        burn_model=args.burn_model,
        thrust_N=300.0,
        isp_s=200.0,
        initial_mass_kg=2000.0,
        max_step_burn_s=0.5
    )

    scenario_config = make_scenario_config(
        h_chaser_km=300.0,
        h_target_km=500.0,
        initial_phase_angle=90.0,
        initial_chaser_angle=0.0,
        desired_rel_lvlh=DESIRED_REL_LVLH
    )

    environment_config = make_environment_config(
        atmospheric_drag=args.atmospheric_drag,
        chaser_cd=args.chaser_cd,
        chaser_area_m2=args.chaser_area_m2,
        target_cd=args.target_cd,
        target_area_m2=args.target_area_m2,
        target_mass_kg=2000.0
    )

    output = minimize_delta_v_on_zero_distance_manifold(
        p_start=None,
        objective_mode="two_impulse_total",
        maneuver_config=maneuver_config,
        scenario_config=scenario_config,
        environment_config=environment_config,
        stage1_verbose=True,
        stage2_verbose=True
    )

    config_path = save_matlab_mission_config(
        output,
        path=args.matlab_config_out,
        phase1_mode="CUSTOM_IMPULSE"
    )
    print(f"\nSaved MATLAB mission config: {config_path}")

    if not args.no_plot:
        plot_constrained_delta_v_history(
            output["history"],
            log_distance=True
        )
