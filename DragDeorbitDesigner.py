import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, root_scalar

from J2PolarHohmann import (
    EarthJ2,
    acceleration_j2_eci,
    atmospheric_drag_acceleration_eci,
    circular_polar_state,
    resolve_environment_config,
)


G0_M_S2 = 9.80665


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


def safe_slug(value):
    text = str(value).strip().lower()
    chars = []
    for ch in text:
        if ch.isalnum() or ch in {"_", "-", "."}:
            chars.append(ch)
        else:
            chars.append("-")
    slug = "".join(chars).strip("-")
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug or "case"


def build_environment_config(args):
    environment = resolve_environment_config({
        "atmospheric_drag": {
            "enabled": True,
            "model": "ISA76",
            "co_rotate_atmosphere": bool(args.co_rotate_atmosphere),
            "earth_rotation_rad_s": args.earth_rotation_rad_s,
            "chaser_cd": args.chaser_cd,
            "chaser_area_m2": args.chaser_area_m2,
            "target_cd": args.target_cd,
            "target_area_m2": args.target_area_m2,
            "target_mass_kg": args.target_mass_kg,
        }
    })
    environment["atmospheric_drag"]["model"] = str(environment["atmospheric_drag"]["model"]).upper()
    return environment


def normalize_burn_model(value):
    model = str(value).strip().lower()
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
        raise ValueError("burn model must be 'impulsive', 'finite_burn', or 'continuous'.")
    return model


def normalize_burn_steering(value):
    steering = str(value).strip().lower().replace("-", "_")
    aliases = {
        "retrograde": "velocity_retrograde",
        "velocity": "velocity_retrograde",
        "velocity_following": "velocity_retrograde",
        "fixed": "fixed_initial_retrograde",
        "fixed_start": "fixed_initial_retrograde",
        "fixed_retrograde": "fixed_initial_retrograde",
    }
    steering = aliases.get(steering, steering)
    if steering not in {"velocity_retrograde", "fixed_initial_retrograde"}:
        raise ValueError("burn steering must be 'velocity_retrograde' or 'fixed_initial_retrograde'.")
    return steering


def finite_burn_mass_and_duration(delta_v_m_s, initial_mass_kg, thrust_N, isp_s):
    if delta_v_m_s < 0:
        raise ValueError("delta-V must be non-negative.")
    if initial_mass_kg <= 0:
        raise ValueError("initial mass must be positive.")
    if thrust_N <= 0:
        raise ValueError("finite burn thrust must be positive.")
    if isp_s <= 0:
        raise ValueError("finite burn Isp must be positive.")

    if delta_v_m_s == 0:
        return initial_mass_kg, 0.0, 0.0

    exhaust_velocity_m_s = isp_s * G0_M_S2
    final_mass_kg = initial_mass_kg * np.exp(-delta_v_m_s / exhaust_velocity_m_s)
    propellant_used_kg = initial_mass_kg - final_mass_kg
    duration_s = propellant_used_kg * exhaust_velocity_m_s / thrust_N
    return final_mass_kg, propellant_used_kg, duration_s


def flight_path_angle_deg(r_km, v_km_s):
    r_hat = r_km / np.linalg.norm(r_km)
    v_radial = float(np.dot(v_km_s, r_hat))
    v_horizontal = float(np.linalg.norm(v_km_s - v_radial * r_hat))
    return float(np.degrees(np.arctan2(v_radial, v_horizontal)))


def single_body_ode(t, y, earth, environment_config, mass_kg):
    r = y[:3]
    v = y[3:]
    acceleration = (
        acceleration_j2_eci(r, earth)
        + atmospheric_drag_acceleration_eci(
            r,
            v,
            mass_kg,
            environment_config,
            vehicle_role="chaser",
            earth=earth,
        )
    )
    return np.concatenate([v, acceleration])


def finite_burn_ode(t, y, earth, environment_config, thrust_N, isp_s, steering, fixed_direction):
    r = y[:3]
    v = y[3:6]
    mass_kg = float(y[6])
    acceleration = (
        acceleration_j2_eci(r, earth)
        + atmospheric_drag_acceleration_eci(
            r,
            v,
            mass_kg,
            environment_config,
            vehicle_role="chaser",
            earth=earth,
        )
    )

    if steering == "velocity_retrograde":
        v_norm = np.linalg.norm(v)
        if v_norm <= 0.0:
            thrust_direction = np.zeros(3)
        else:
            thrust_direction = -v / v_norm
    else:
        thrust_direction = fixed_direction

    thrust_acceleration = (thrust_N / mass_kg) / 1000.0 * thrust_direction
    mass_flow_kg_s = -thrust_N / (isp_s * G0_M_S2)
    return np.concatenate([v, acceleration + thrust_acceleration, [mass_flow_kg_s]])


def entry_altitude_event_factory(earth, entry_interface_altitude_km):
    def entry_altitude_event(t, y):
        return np.linalg.norm(y[:3]) - earth.radius - entry_interface_altitude_km

    entry_altitude_event.terminal = True
    entry_altitude_event.direction = -1
    return entry_altitude_event


def summarize_entry_result(result, entry_state, entry_time_s, earth):
    result.update({
        "coast_time_s": float(entry_time_s - result.get("burn_time_s", 0.0)),
        "total_time_s": float(entry_time_s),
        "entry_altitude_km": float(np.linalg.norm(entry_state[:3]) - earth.radius),
        "entry_fpa_deg": flight_path_angle_deg(entry_state[:3], entry_state[3:6]),
        "entry_speed_m_s": float(np.linalg.norm(entry_state[3:6]) * 1000.0),
        "entry_radius_km": float(np.linalg.norm(entry_state[:3])),
    })
    return result


def propagate_drag_deorbit(
    delta_v_m_s,
    start_altitude_km,
    entry_interface_altitude_km,
    environment_config,
    mass_kg,
    max_coast_time_s,
    max_step_s,
    rtol,
    atol,
    burn_model,
    burn_steering,
    finite_burn_thrust_N,
    finite_burn_isp_s,
    finite_burn_max_step_s,
    start_argument_deg=0.0,
    earth=EarthJ2(),
):
    state0 = circular_polar_state(
        earth.radius + start_altitude_km,
        np.radians(start_argument_deg),
        earth=earth,
    )
    v_hat = state0[3:6] / np.linalg.norm(state0[3:6])
    fixed_retrograde_direction = -v_hat
    burn_model = normalize_burn_model(burn_model)
    burn_steering = normalize_burn_steering(burn_steering)
    post_burn_mass_kg, propellant_used_kg, burn_duration_s = finite_burn_mass_and_duration(
        delta_v_m_s=delta_v_m_s,
        initial_mass_kg=mass_kg,
        thrust_N=finite_burn_thrust_N,
        isp_s=finite_burn_isp_s,
    )

    result = {
        "burn_model": burn_model.upper(),
        "burn_steering": burn_steering,
        "delta_v_m_s": float(delta_v_m_s),
        "commanded_delta_v_m_s": float(delta_v_m_s),
        "initial_mass_kg": float(mass_kg),
        "predicted_post_burn_mass_kg": float(post_burn_mass_kg),
        "predicted_propellant_used_kg": float(propellant_used_kg),
        "predicted_burn_duration_s": float(0.0 if burn_model == "impulsive" else burn_duration_s),
        "delivered_delta_v_m_s": float(delta_v_m_s),
        "initial_burn_direction_eci": to_jsonable(fixed_retrograde_direction),
    }

    entry_altitude_event = entry_altitude_event_factory(earth, entry_interface_altitude_km)

    if burn_model == "impulsive":
        state_after_burn = state0.copy()
        state_after_burn[3:6] = state_after_burn[3:6] + (delta_v_m_s / 1000.0) * fixed_retrograde_direction
        coast_mass_kg = post_burn_mass_kg
        result.update({
            "burn_completed": True,
            "burn_time_s": 0.0,
            "burn_solver_success": True,
            "burn_samples": 1,
            "post_burn_altitude_km": float(np.linalg.norm(state_after_burn[:3]) - earth.radius),
            "post_burn_fpa_deg": flight_path_angle_deg(state_after_burn[:3], state_after_burn[3:6]),
            "post_burn_speed_m_s": float(np.linalg.norm(state_after_burn[3:6]) * 1000.0),
        })
    else:
        burn_state0 = np.concatenate([state0, [mass_kg]])
        sol_burn = solve_ivp(
            lambda t, y: finite_burn_ode(
                t,
                y,
                earth,
                environment_config,
                finite_burn_thrust_N,
                finite_burn_isp_s,
                burn_steering,
                fixed_retrograde_direction,
            ),
            (0.0, burn_duration_s),
            burn_state0,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            max_step=finite_burn_max_step_s,
            events=entry_altitude_event,
        )

        if sol_burn.t_events[0].size > 0:
            entry_state_aug = sol_burn.y_events[0][0]
            actual_burn_time_s = float(sol_burn.t_events[0][0])
            final_mass_kg = float(entry_state_aug[6])
            delivered_delta_v_m_s = finite_burn_isp_s * G0_M_S2 * np.log(mass_kg / final_mass_kg)
            result.update({
                "reached_interface": True,
                "solver_success": bool(sol_burn.success),
                "solver_message": str(sol_burn.message),
                "samples": int(sol_burn.t.size),
                "burn_completed": False,
                "burn_time_s": actual_burn_time_s,
                "burn_solver_success": bool(sol_burn.success),
                "burn_samples": int(sol_burn.t.size),
                "predicted_post_burn_mass_kg": final_mass_kg,
                "predicted_propellant_used_kg": float(mass_kg - final_mass_kg),
                "delivered_delta_v_m_s": float(delivered_delta_v_m_s),
            })
            return summarize_entry_result(result, entry_state_aug[:6], actual_burn_time_s, earth)

        if not sol_burn.success:
            result.update({
                "reached_interface": False,
                "solver_success": False,
                "solver_message": str(sol_burn.message),
                "samples": int(sol_burn.t.size),
                "burn_completed": False,
                "burn_time_s": float(sol_burn.t[-1]),
                "burn_solver_success": False,
                "burn_samples": int(sol_burn.t.size),
            })
            return result

        burn_state_final = sol_burn.y[:, -1]
        state_after_burn = burn_state_final[:6]
        coast_mass_kg = float(burn_state_final[6])
        delivered_delta_v_m_s = finite_burn_isp_s * G0_M_S2 * np.log(mass_kg / coast_mass_kg)
        result.update({
            "burn_completed": True,
            "burn_time_s": float(sol_burn.t[-1]),
            "burn_solver_success": bool(sol_burn.success),
            "burn_samples": int(sol_burn.t.size),
            "predicted_post_burn_mass_kg": coast_mass_kg,
            "predicted_propellant_used_kg": float(mass_kg - coast_mass_kg),
            "delivered_delta_v_m_s": float(delivered_delta_v_m_s),
            "post_burn_altitude_km": float(np.linalg.norm(state_after_burn[:3]) - earth.radius),
            "post_burn_fpa_deg": flight_path_angle_deg(state_after_burn[:3], state_after_burn[3:6]),
            "post_burn_speed_m_s": float(np.linalg.norm(state_after_burn[3:6]) * 1000.0),
        })

    sol = solve_ivp(
        lambda t, y: single_body_ode(t, y, earth, environment_config, coast_mass_kg),
        (0.0, max_coast_time_s),
        state_after_burn,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=max_step_s,
        events=entry_altitude_event,
    )

    result.update({
        "reached_interface": bool(sol.t_events[0].size > 0),
        "solver_success": bool(sol.success),
        "solver_message": str(sol.message),
        "samples": int(sol.t.size),
    })

    if not result["reached_interface"]:
        final_state = sol.y[:, -1]
        result.update({
            "coast_time_s": float(sol.t[-1]),
            "total_time_s": float(result["burn_time_s"] + sol.t[-1]),
            "final_altitude_km": float(np.linalg.norm(final_state[:3]) - earth.radius),
            "final_fpa_deg": flight_path_angle_deg(final_state[:3], final_state[3:6]),
            "final_speed_m_s": float(np.linalg.norm(final_state[3:6]) * 1000.0),
        })
        return result

    entry_state = sol.y_events[0][0]
    return summarize_entry_result(result, entry_state, result["burn_time_s"] + sol.t_events[0][0], earth)


def solve_drag_deorbit(args):
    environment_config = build_environment_config(args)
    cache = {}

    def evaluate(delta_v_m_s):
        key = round(float(delta_v_m_s), 6)
        if key not in cache:
            cache[key] = propagate_drag_deorbit(
                delta_v_m_s=key,
                start_altitude_km=args.start_altitude_km,
                entry_interface_altitude_km=args.entry_interface_altitude_km,
                environment_config=environment_config,
                mass_kg=args.chaser_mass_kg,
                max_coast_time_s=args.max_coast_time_s,
                max_step_s=args.max_step_s,
                rtol=args.rtol,
                atol=args.atol,
                burn_model=args.normalized_burn_model,
                burn_steering=args.normalized_burn_steering,
                finite_burn_thrust_N=args.finite_burn_thrust_N,
                finite_burn_isp_s=args.finite_burn_isp_s,
                finite_burn_max_step_s=args.finite_burn_max_step_s,
                start_argument_deg=args.start_argument_deg,
            )
        return cache[key]

    target_fpa = abs(float(args.flight_path_angle_deg))

    def fpa_error(delta_v_m_s):
        result = evaluate(delta_v_m_s)
        if not result["reached_interface"]:
            return np.nan
        return -result["entry_fpa_deg"] - target_fpa

    fixed_delta_v = args.delta_v_m_s is not None
    if fixed_delta_v:
        best_delta_v = float(args.delta_v_m_s)
        best_result = evaluate(best_delta_v)
        method = "fixed_delta_v"
    else:
        grid = np.linspace(args.min_delta_v_m_s, args.max_delta_v_m_s, args.search_samples)
        valid = []
        for delta_v in grid:
            result = evaluate(delta_v)
            if result["reached_interface"]:
                valid.append((float(delta_v), fpa_error(delta_v)))

        if not valid:
            raise RuntimeError(
                "No delta-V in the search interval reached the entry interface. "
                "Increase --max-delta-v-m-s or --max-coast-time-s."
            )

        bracket = None
        for (dv_a, err_a), (dv_b, err_b) in zip(valid[:-1], valid[1:]):
            if err_a == 0.0:
                bracket = (dv_a, dv_a)
                break
            if np.sign(err_a) != np.sign(err_b):
                bracket = (dv_a, dv_b)
                break

        if bracket is not None and bracket[0] != bracket[1]:
            root = root_scalar(lambda dv: fpa_error(dv), bracket=bracket, xtol=args.delta_v_tol_m_s)
            best_delta_v = float(root.root)
            method = "target_fpa_root"
        elif bracket is not None:
            best_delta_v = float(bracket[0])
            method = "target_fpa_grid_exact"
        else:
            min_valid_delta_v = min(dv for dv, _ in valid)

            def objective(delta_v):
                result = evaluate(delta_v)
                if not result["reached_interface"]:
                    return 1e6 + abs(delta_v - args.max_delta_v_m_s)
                return abs(-result["entry_fpa_deg"] - target_fpa)

            opt = minimize_scalar(
                objective,
                bounds=(min_valid_delta_v, args.max_delta_v_m_s),
                method="bounded",
                options={"xatol": args.delta_v_tol_m_s},
            )
            best_delta_v = float(opt.x)
            method = "target_fpa_min_abs_error"

        best_result = evaluate(best_delta_v)

    rule_result = evaluate(args.rule_of_thumb_delta_v_m_s)
    return environment_config, method, best_result, rule_result, len(cache)


def relative_posix_path(path, base_dir):
    try:
        return path.resolve().relative_to(base_dir.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def attach_archive_metadata(config):
    created_at = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    settings = {
        "schema_version": config["schema_version"],
        "source": config["source"],
        "scenario": config["scenario"],
        "maneuver": config["maneuver"],
        "phase3": {
            "entry_interface_altitude_km": config["phase3"]["entry_interface_altitude_km"],
            "flight_path_angle_deg": config["phase3"]["flight_path_angle_deg"],
        },
        "environment": config["environment"],
    }
    result = {
        "settings": settings,
        "drag_deorbit": config["phase3"]["drag_deorbit"],
    }
    settings_hash = short_hash(settings, 16)
    result_hash = short_hash(result, 16)
    case_id = (
        f"drag_deorbit_h{safe_slug(config['scenario']['h_start_km'])}"
        f"_ei{safe_slug(config['phase3']['entry_interface_altitude_km'])}"
        f"_fpa{safe_slug(config['phase3']['flight_path_angle_deg'])}"
        f"_{settings_hash[:8]}_{result_hash[:8]}_"
        f"{created_at.replace('-', '').replace(':', '').replace('T', '_').replace('Z', 'Z')}"
    )
    config["archive"] = {
        "schema_version": 1,
        "created_at": created_at,
        "case_id": case_id,
        "settings_hash": settings_hash,
        "result_hash": result_hash,
        "settings_signature": settings,
    }
    return config


def update_index(index_path, config, archive_path, latest_path):
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

    archive = config["archive"]
    drag = config["phase3"]["drag_deorbit"]
    entry = {
        "case_id": archive["case_id"],
        "settings_hash": archive["settings_hash"],
        "result_hash": archive["result_hash"],
        "created_at": archive["created_at"],
        "path": relative_posix_path(archive_path, config_dir),
        "latest_alias_path": relative_posix_path(latest_path, config_dir),
        "scenario": config["scenario"],
        "phase3": {
            "mode": config["phase3"]["mode"],
            "entry_interface_altitude_km": config["phase3"]["entry_interface_altitude_km"],
            "flight_path_angle_deg": config["phase3"]["flight_path_angle_deg"],
            "burn_model": drag["burn_model"],
            "burn_steering": drag["burn_steering"],
            "drag_deorbit_delta_v_m_s": drag["delta_v_m_s"],
            "predicted_burn_duration_s": drag["predicted_burn_duration_s"],
            "predicted_entry_fpa_deg": drag["predicted_entry_fpa_deg"],
            "predicted_coast_time_s": drag["predicted_coast_time_s"],
        },
    }

    cases = [case for case in cases if case.get("case_id") != entry["case_id"]]
    cases.append(entry)
    cases.sort(key=lambda case: str(case.get("created_at", "")))
    index = {
        "schema_version": 1,
        "source": "DragDeorbitDesigner.py",
        "updated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "latest_case_id": entry["case_id"],
        "cases": cases,
    }
    index_path.write_text(json.dumps(index, indent=2, sort_keys=True), encoding="utf-8")


def write_outputs(config, output_path):
    output_path = Path(output_path)
    config_dir = output_path.parent
    archive_dir = config_dir / "drag_deorbit_runs"
    archive_dir.mkdir(parents=True, exist_ok=True)

    config = attach_archive_metadata(to_jsonable(config))
    archive_path = archive_dir / f"{config['archive']['case_id']}.json"
    archive_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    update_index(config_dir / "drag_deorbit_solution_index.json", config, archive_path, output_path)
    return output_path, archive_path


def build_matlab_config(args, environment_config, method, best_result, rule_result, evaluations):
    target_fpa = abs(float(args.flight_path_angle_deg))
    fpa_error = None
    if best_result["reached_interface"]:
        fpa_error = -best_result["entry_fpa_deg"] - target_fpa

    return {
        "schema_version": 1,
        "source": "DragDeorbitDesigner.py",
        "scenario": {
            "h_start_km": args.start_altitude_km,
            "start_argument_deg": args.start_argument_deg,
        },
        "environment": environment_config,
        "maneuver": {
            "burn_model": args.normalized_burn_model.upper(),
            "burn_steering": args.normalized_burn_steering,
            "initial_mass_kg": args.chaser_mass_kg,
            "finite_burn_thrust_N": args.finite_burn_thrust_N,
            "finite_burn_isp_s": args.finite_burn_isp_s,
            "finite_burn_dt_s": args.finite_burn_dt_s,
            "finite_burn_max_step_s": args.finite_burn_max_step_s,
        },
        "phase3": {
            "mode": "HOHMANN",
            "entry_interface_altitude_km": args.entry_interface_altitude_km,
            "flight_path_angle_deg": args.flight_path_angle_deg,
            "drag_deorbit": {
                "enabled": True,
                "design_method": method,
                "burn_model": best_result["burn_model"],
                "burn_direction": "retrograde",
                "burn_steering": best_result["burn_steering"],
                "delta_v_m_s": best_result["delta_v_m_s"],
                "commanded_delta_v_m_s": best_result["commanded_delta_v_m_s"],
                "delivered_delta_v_m_s": best_result["delivered_delta_v_m_s"],
                "predicted_burn_duration_s": best_result["predicted_burn_duration_s"],
                "predicted_burn_time_s": best_result["burn_time_s"],
                "burn_completed": best_result["burn_completed"],
                "initial_mass_kg": best_result["initial_mass_kg"],
                "predicted_post_burn_mass_kg": best_result["predicted_post_burn_mass_kg"],
                "predicted_propellant_used_kg": best_result["predicted_propellant_used_kg"],
                "finite_burn_thrust_N": args.finite_burn_thrust_N,
                "finite_burn_isp_s": args.finite_burn_isp_s,
                "finite_burn_dt_s": args.finite_burn_dt_s,
                "finite_burn_max_step_s": args.finite_burn_max_step_s,
                "post_burn_altitude_km": best_result.get("post_burn_altitude_km"),
                "post_burn_fpa_deg": best_result.get("post_burn_fpa_deg"),
                "post_burn_speed_m_s": best_result.get("post_burn_speed_m_s"),
                "predicted_coast_time_s": best_result.get("coast_time_s"),
                "predicted_total_time_s": best_result.get("total_time_s"),
                "predicted_entry_altitude_km": best_result.get("entry_altitude_km"),
                "predicted_entry_fpa_deg": best_result.get("entry_fpa_deg"),
                "predicted_entry_speed_m_s": best_result.get("entry_speed_m_s"),
                "fpa_error_deg": fpa_error,
                "max_coast_time_s": args.max_coast_time_s,
                "search_min_delta_v_m_s": args.min_delta_v_m_s,
                "search_max_delta_v_m_s": args.max_delta_v_m_s,
                "rule_of_thumb_delta_v_m_s": args.rule_of_thumb_delta_v_m_s,
                "rule_of_thumb_result": rule_result,
            },
        },
        "optimizer": {
            "success": bool(best_result["reached_interface"]),
            "method": method,
            "evaluations": evaluations,
            "target_fpa_deg": target_fpa,
            "fpa_error_deg": fpa_error,
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="Design a drag-aware single-retrograde-burn LEO deorbit for MATLAB Phase 3."
    )
    parser.add_argument("--matlab-config-out", default="configs/latest_drag_deorbit_solution.json")
    parser.add_argument("--start-altitude-km", type=float, default=500.0)
    parser.add_argument("--entry-interface-altitude-km", type=float, default=120.0)
    parser.add_argument("--flight-path-angle-deg", type=float, default=4.0)
    parser.add_argument("--start-argument-deg", type=float, default=0.0)
    parser.add_argument("--burn-model", default="impulsive",
                        help="'impulsive' applies an instantaneous velocity decrement. "
                             "'finite_burn' and 'continuous' integrate an optional finite-duration burn.")
    parser.add_argument("--burn-steering", default="velocity_retrograde",
                        help="Finite-burn steering law: 'velocity_retrograde' or 'fixed_initial_retrograde'.")
    parser.add_argument("--finite-burn-thrust-n", "--finite-burn-thrust-N",
                        dest="finite_burn_thrust_N", type=float, default=300.0)
    parser.add_argument("--finite-burn-isp-s", dest="finite_burn_isp_s", type=float, default=200.0)
    parser.add_argument("--finite-burn-dt-s", dest="finite_burn_dt_s", type=float, default=1.0,
                        help="MATLAB RK4 execution step saved to JSON for finite deorbit burns.")
    parser.add_argument("--finite-burn-max-step-s", dest="finite_burn_max_step_s", type=float, default=2.0,
                        help="Python solve_ivp max step during the finite burn.")
    parser.add_argument("--delta-v-m-s", type=float, default=None)
    parser.add_argument("--min-delta-v-m-s", type=float, default=40.0)
    parser.add_argument("--max-delta-v-m-s", type=float, default=320.0)
    parser.add_argument("--delta-v-tol-m-s", type=float, default=0.05)
    parser.add_argument("--search-samples", type=int, default=57)
    parser.add_argument("--rule-of-thumb-delta-v-m-s", type=float, default=100.0)
    parser.add_argument("--max-coast-time-s", type=float, default=6 * 3600.0)
    parser.add_argument("--max-step-s", type=float, default=10.0)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-11)
    parser.add_argument("--chaser-mass-kg", type=float, default=2000.0)
    parser.add_argument("--target-mass-kg", type=float, default=2000.0)
    parser.add_argument("--chaser-cd", type=float, default=2.2)
    parser.add_argument("--chaser-area-m2", type=float, default=4.0)
    parser.add_argument("--target-cd", type=float, default=2.2)
    parser.add_argument("--target-area-m2", type=float, default=4.0)
    parser.add_argument("--earth-rotation-rad-s", type=float, default=7.2921159e-5)
    parser.add_argument("--co-rotate-atmosphere", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    args.normalized_burn_model = normalize_burn_model(args.burn_model)
    args.normalized_burn_steering = normalize_burn_steering(args.burn_steering)

    environment_config, method, best_result, rule_result, evaluations = solve_drag_deorbit(args)
    config = build_matlab_config(args, environment_config, method, best_result, rule_result, evaluations)
    latest_path, archive_path = write_outputs(config, args.matlab_config_out)

    print("========== Drag-Aware Deorbit Design ==========")
    print(f"Start altitude          : {args.start_altitude_km:.3f} km")
    print(f"Entry interface         : {args.entry_interface_altitude_km:.3f} km")
    print(f"Target FPA magnitude    : {abs(args.flight_path_angle_deg):.3f} deg")
    print(f"Burn model              : {args.normalized_burn_model.upper()} ({args.normalized_burn_steering})")
    print(f"Design method           : {method}")
    print(f"Retrograde delta-V      : {best_result['delta_v_m_s']:.3f} m/s")
    if args.normalized_burn_model == "finite_burn":
        print(f"Finite burn duration    : {best_result['predicted_burn_duration_s'] / 60.0:.3f} min")
        print(f"Finite burn propellant  : {best_result['predicted_propellant_used_kg']:.3f} kg")
    if best_result["reached_interface"]:
        print(f"Predicted interface time: {best_result['total_time_s'] / 60.0:.3f} min")
        print(f"Predicted interface FPA : {best_result['entry_fpa_deg']:.4f} deg")
        print(f"Predicted interface v   : {best_result['entry_speed_m_s']:.3f} m/s")
    else:
        print("Predicted interface     : not reached within max coast time")
    print(f"100 m/s check reached   : {rule_result['reached_interface']}")
    if rule_result["reached_interface"]:
        print(f"100 m/s check FPA       : {rule_result['entry_fpa_deg']:.4f} deg")
    print(f"Saved latest JSON       : {latest_path}")
    print(f"Saved archived JSON     : {archive_path}")


if __name__ == "__main__":
    main()
