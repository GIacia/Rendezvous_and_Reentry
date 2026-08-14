# Unmanned Rendezvous Mission Simulator

MATLAB-based end-to-end rendezvous mission simulator for an unmanned chaser spacecraft operating in Earth orbit. The current codebase combines:

- **Phase 1**: 3-DOF phasing / homing with impulsive maneuvers and J2-aware orbital propagation
- **Phase 2**: LVLH waypoint-impulse proximity operations with cycloidal drift, R-bar hops, nonlinear J2 propagation, and mass bookkeeping
- **Phase 3**: selectable 3-DOF de-orbit / re-entry modes, including direct FPA-targeted descent and a configurable parking-orbit / R-bar-aligned re-entry setup
- **Phase 4**: atmospheric entry from a configurable entry interface, with re-entry vehicle shape selection, heat-rate diagnostics, and chaser-to-entry-vehicle line-of-sight checks

The project is designed as a mission-level simulation framework rather than a single guidance-law demo. It is useful for studying how orbit transfer logic, target-relative geometry, J2 perturbation, mass depletion, proximity-operation sequencing, and entry diagnostics interact in one connected workflow.

## Current Scope

This repository currently models:

- Earth central gravity + J2 perturbation
- Optional atmospheric drag using an ISA76-style standard-atmosphere density model
- A target spacecraft and a chaser spacecraft
- Hohmann-based phase planning with target co-propagation
- A J2-aware wait-time search before departure so the transfer arrives closer to a desired LVLH capture point
- Custom phased maneuver logic driven by externally tuned phase angle, delta-V, and gamma parameters, with impulsive execution by default and finite-burn execution kept as an explicit study option
- LVLH waypoint-impulse proximity operations with cleanup, hold trims, cycloidal drift, R-bar hops, braking impulses, and mass depletion
- selectable Phase 3 re-entry logic through `Mission_Run_Config.m`
- Simple thrust uncertainty / noise injection in selected phasing modes
- Mission-level delta-V and propellant budget tracking
- Visualization of trajectory, altitude history, proximity trajectory, mass depletion, atmospheric-entry heat flux, dynamic pressure, g-load, and line-of-sight margin

## Repository Structure

```text
.
|-- Main_Mission_Simulator.m      % Main entry point for the full mission
|-- Mission_Config.m              % Physical constants, vehicle data, mission parameters
|-- Mission_Run_Config.m          % User-facing run control panel
|-- Phasing_Propagator.m          % 3-DOF phasing, Hohmann, custom impulse logic, preliminary Multi-Hohmann branch
|-- Env_EOM.m                     % Environment and 3-DOF/6-DOF equations of motion
|-- Atmospheric_Drag_Acceleration.m % Shared MATLAB atmospheric-drag model
|-- Standard_Atmosphere_Density.m % Shared ISA76-style atmosphere helper
|-- Reentry_Propagator.m          % Atmospheric-entry propagation and LOS diagnostics
|-- Run_Paper_Reproduction_Suite.m % Zhang/Saito standalone audit entry point
|-- +reentry_core/                % Shared reusable atmospheric-entry physics
|-- +paperstudies/
|   |-- +zhang/                   % Zhang equations, conditions, provenance, surrogate adapter
|   `-- +saito/                   % Saito tables, equations, grids, surrogate adapter
|-- J2PolarHohmann.py             % Python J2 polar Hohmann propagation study
|-- J2PolarHohmannShooting.py     % Python shooting / constrained optimization helper
|-- DragDeorbitDesigner.py        % Python drag-aware Phase 3 deorbit design helper
|-- docs/
|   |-- ARCHITECTURE.md
|   |-- ASSUMPTIONS_AND_LIMITATIONS.md
|   |-- PAPER_REPRODUCTION_FRAMEWORK.md
|   `-- REENTRY_PAPER_MODELS.md
|-- validation/
|   |-- Check_Reentry_Code.m
|   |-- Run_All_Reentry_Validations.m
|   |-- Validate_Reentry_Core_Equivalence.m
|   |-- Validate_Reentry_Propagator.m
|   `-- Validate_Paper_Reproduction.m
|-- examples/
|   `-- README.md
|-- results/
|   `-- README.md
|-- legacy/                       % Older MATLAB versions, including GNC_Controller.m
|-- requirements.txt              % Python helper dependencies
`-- .gitignore
```

The active MATLAB `.m` files are intentionally kept at the repository root so that the existing MATLAB execution flow does not break.

## How to Run

1. Open the repository folder in MATLAB.
2. Make sure the current folder is the repository root.
3. Edit `Mission_Run_Config.m` if you want to change the scenario, burn model,
   re-entry mode, parking altitude, entry vehicle, or solver tolerances.
4. Run:

```matlab
Main_Mission_Simulator
```

`Mission_Run_Config.m` is the recommended control surface for normal use. It
keeps the commonly adjusted settings in one place:

- Python optimizer JSON selection: `AUTO`, `NONE`, `FILE`, `CASE_ID`, or `HASH`
- scenario altitudes and initial phase geometry
- maneuver model: `IMPULSIVE` or `FINITE_BURN`
- Phase 1 phasing parameters and capture tolerances
- Phase 2 proximity-operation waypoints and timing
- Phase 3 mode, parking altitude, entry-interface altitude, and FPA
- Phase 4 re-entry vehicle shape and atmosphere-entry settings
- optional orbital atmospheric drag

Common examples:

```matlab
% In Mission_Run_Config.m
run.maneuver.burn_model = "IMPULSIVE";
run.phase3.mode = "R_BAR_200_FPA";
run.phase3.parking_altitude_km = 200;
run.phase3.flight_path_angle_deg = 4;
run.reentry.shape = "TPS_MIN";
run.environment.atmospheric_drag.enabled = true;
run.environment.atmospheric_drag.apply_from_phase = "PHASE3";
```

Paper-derived reduced-order atmospheric-entry modes can be selected directly:

```matlab
% RLV polynomial aerodynamics, speed-scheduled AoA, and target-relay
% antenna/path-constraint diagnostics.
run.reentry.vehicle_mode = "SPACEPLANE";
run.reentry.communication.antenna_mount = "AFT"; % or PAPER_TOP
run.reentry.communication.relay_mode = "MISSION_TARGET_DYNAMIC_ORBIT";
% Optional paper relay reference:
% run.reentry.communication.relay_mode = "PAPER_TDRS_STATIC_EARTH_FIXED";
% run.reentry.communication.earth_fixed_to_eci_angle_at_mission_epoch_deg = 0;
% Optional only after a link-budget-derived range is available:
run.reentry.communication.max_range_m = inf;

% 60 kg capsule separated at the 120 km entry interface, with the paper
% reference FPA and parachute-speed terminal event.
run.reentry.vehicle_mode = "CAPSULE";
run.reentry.capsule.separation_mode = "ENTRY_INTERFACE";
run.reentry.capsule.add_to_chaser_initial_mass = true;
run.reentry.capsule.altitude_termination_enabled = false; % use 240 m/s event
```

See `docs/REENTRY_PAPER_MODELS.md` for the implemented equations, paper
values, uncertainty scales, and the boundary between diagnostics and guidance.

### Standalone paper reproduction studies

The integrated mission modes above are reduced mission-level models. Separate
Zhang and Saito study packages preserve the published conditions without
silently mixing them with the normal mission scenario:

```matlab
% Fast deterministic audit; no optimization or forward propagation.
report = Run_Paper_Reproduction_Suite();

% Optional forward runs through the same physical kernel used by the mission.
% Missing paper inputs remain explicitly labeled surrogate assumptions.
report = Run_Paper_Reproduction_Suite(struct('forward', true));
```

The paper suite runs selected source-anchor, equation, state-construction,
coordinate-geometry, and uncertainty-grid regressions;
`Run_All_Reentry_Validations` additionally checks shared-kernel consistency
and both adapters. The tables remain traceable and machine-readable, but the
automated suite is not independent proof of every transcribed cell. These checks
do not claim full numerical reproduction of either paper: Zhang omits material
vehicle/OCP data, while Saito omits the proprietary aerodynamic database and
several guidance/controller parameters. See
`docs/PAPER_REPRODUCTION_FRAMEWORK.md` for the evidence matrix, exact blockers,
commands, and interpretation rules.

```matlab
addpath validation
Run_All_Reentry_Validations
```

The entry attitude inputs are intentionally simple in the current model.
SPACEPLANE uses its configured speed-scheduled AoA unless
`run.reentry.aoa_deg` supplies a constant override; CAPSULE uses its constant
trim surrogate. `run.reentry.bank_angle_deg` is held constant in both modes.
The simulator does not yet propagate vehicle attitude, solve trim, or apply
closed-loop bank guidance for the separated re-entry vehicle.

Leave Phase 1 optimizer outputs such as `run.phase1.phase_angle_deg`,
`run.phase1.delta_v_m_s`, and `run.phase1.gamma_deg` as `[]` when you want
MATLAB to use the selected Python optimizer JSON values. Fill them in only when
you want to manually override the optimizer solution.

Python shooting results can be exported directly into a MATLAB-readable mission config:

```bash
python J2PolarHohmannShooting.py --no-plot --matlab-config-out configs/latest_python_solution.json
```

Every Python export now writes three things:

- `configs/latest_python_solution.json`: latest alias for backward compatibility
- `configs/python_runs/<case_id>.json`: archived result that is never overwritten by a different case
- `configs/python_solution_index.json`: searchable index with `case_id`, `settings_hash`, burn model, scenario, drag settings, and optimizer summary

Atmospheric drag is off by default. To optimize with the ISA76 drag option and export the same environment settings to MATLAB:

```bash
python J2PolarHohmannShooting.py --no-plot --atmospheric-drag isa76 --matlab-config-out configs/latest_python_solution.json
```

The default run-control scope is
`run.environment.atmospheric_drag.apply_from_phase = "PHASE3"`. With that
temporary setting, Phase 1/2 remain gravity + J2 only so the existing drag-off
Phase 1 optimizer archives stay usable, and orbital drag is enabled after
berthing for Phase 3 and the orbiting relay chaser during Phase 4. Set
`apply_from_phase = "PHASE1"` only when you also have a matching drag-on Phase 1
optimizer JSON. The separated re-entry vehicle always uses the atmospheric
entry model in `Reentry_Propagator.m`; this orbital-drag switch does not disable
Phase 4 entry drag/lift.

When orbital drag is enabled and Phase 3 is `HOHMANN`, MATLAB uses a separate
Python-generated drag-aware deorbit design instead of the legacy dragless conic
FPA injection. Generate or refresh that design with:

```bash
python DragDeorbitDesigner.py --start-altitude-km 500 --entry-interface-altitude-km 120 --flight-path-angle-deg 4 --matlab-config-out configs/latest_drag_deorbit_solution.json
```

The generated file is loaded through:

```matlab
% In Mission_Run_Config.m
run.phase3.drag_deorbit_design.mode = "AUTO";
run.phase3.drag_deorbit_design.file = "configs/latest_drag_deorbit_solution.json";
```

By default this design helper writes an impulsive retrograde deorbit design
followed by J2 + atmospheric-drag propagation to the configured entry
interface. `--burn-model finite_burn` remains available only as an explicit
experiment; it is not the recommended default and is not a full low-thrust
trajectory optimizer. `--burn-model continuous` is accepted only as an alias
for that finite-burn experiment.

The common "about 100 m/s from LEO" rule of thumb is valid for lower start
altitudes and shallower entry-interface angles, but it is not universal. For
the repository's current 500 km / 120 km / 4 deg default, the generated
impulsive design is about 275 m/s; the same JSON records whether a 100 m/s
check reaches the interface. A low-thrust finite burn such as 300 N can last
long enough that the vehicle reaches the interface before the commanded burn
completes, so it is not used as the default deorbit method.

Then load that exact latest config by setting the run control file:

```matlab
% In Mission_Run_Config.m
run.python_config.mode = "FILE";
run.python_config.file = "configs/latest_python_solution.json";
Main_Mission_Simulator
```

Or let MATLAB select a matching archived Python result automatically:

```matlab
% In Mission_Run_Config.m
run.python_config.mode = "AUTO";
run.maneuver.burn_model = "IMPULSIVE";
Main_Mission_Simulator
```

The automatic selector reads `configs/python_solution_index.json` and chooses
the best archived case for the current burn model, scenario altitude pair, and
atmospheric-drag setting. The burn model is treated as a hard filter; scenario
and drag settings are scored, so use `FILE`, `CASE_ID`, or `HASH` when exact
reproduction matters. You can also pin an exact archived run:

```matlab
% In Mission_Run_Config.m
run.python_config.mode = "CASE_ID";
run.python_config.case_id = "impulsive_drag-off_h300.0-500.0_phase90.0_31c1e45c_410e0069_20260630_134559Z";
Main_Mission_Simulator
```

The JSON file overrides only the fields it contains. `Mission_Run_Config.m` is
then applied again as the final user-facing override, while `Mission_Config.m`
remains the source for physical constants and vehicle shape tables.

For a MATLAB-only run that ignores Python optimizer JSON:

```matlab
% In Mission_Run_Config.m
run.python_config.mode = "NONE";
Main_Mission_Simulator
```

The active MATLAB path is intended to use standard MATLAB numerical functionality. The Python helper scripts require the packages listed in `requirements.txt`.

For `R_BAR_200_FPA`, the final parking-orbit to entry-interface injection fuel
update is controlled independently. The default is `false`; to include that burn
in the propellant and remaining-mass budget:

```matlab
run.phase3.charge_final_reentry_fuel = true; % research default
Main_Mission_Simulator
```

Legacy environment-variable overrides are still supported for batch scripts:
`RENDEZVOUS_CONFIG_JSON`, `RENDEZVOUS_CONFIG_CASE_ID`,
`RENDEZVOUS_CONFIG_HASH`, `RENDEZVOUS_BURN_MODEL`,
`RENDEZVOUS_ATMOSPHERIC_DRAG`, `RENDEZVOUS_REENTRY_SHAPE`,
`RENDEZVOUS_PHASE3_MODE`, and `RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL`. Set
`run.runtime.allow_environment_overrides = false` to make `Mission_Run_Config.m`
the only run-control source.

## Simulation Flow

### Phase 1: Phasing & Homing

`Phasing_Propagator.m` performs the pre-rendezvous orbital transfer.

The active root script currently uses `CUSTOM_IMPULSE` mode by default, where
the departure phase angle, burn magnitude, and burn flight-path tilt are supplied
through the selected Python optimizer JSON or through `Mission_Run_Config.m`.
The maneuver execution model is selected with `run.maneuver.burn_model`:
`"IMPULSIVE"` is the default instantaneous delta-V model. `"FINITE_BURN"`
applies the same delta-V direction over a finite-duration firing using the
configured thrust, Isp, and mass depletion, but it should be selected
explicitly for finite-burn sensitivity studies.

In `HOHMANN` mode, the logic is:

1. Co-propagate the chaser and target under J2.
2. Search for a departure wait time that best aligns the transfer arrival with a desired LVLH capture point.
3. Execute the first Hohmann impulse.
4. Numerically propagate the transfer arc.
5. Circularize at the destination altitude.

### Phase 2: Proximity Operations

The current root-level mission script performs proximity operations directly in `Main_Mission_Simulator.m`.

The active Phase 2 implementation includes:

- LVLH relative position / velocity computation
- optional cleanup to the S2 -V-bar hold point
- S2 hold trim to cancel residual relative velocity
- one S2-to-S3 V-bar impulse followed by natural cycloidal drift
- S3 detection at the first V-bar crossing
- CW/Hill waypoint targeting for S3-to-S4 R-bar hops
- nonlinear free propagation with `Env_EOM.m`
- braking impulses and propellant bookkeeping at hold points

Older force-based `GNC_Controller.m` implementations are preserved under `legacy/`, but they are not called by the current root-level mission script.

### Phase 3: Re-entry / Descent

Phase 3 is selected with `run.phase3.mode` in `Mission_Run_Config.m`.

Current modes:

- `HOHMANN`: performs a direct FPA-targeted de-orbit injection and stops at `run.phase3.entry_interface_altitude_km`. It no longer performs a nonphysical circularization after crossing the interface.
- `R_BAR_200_FPA`: first lowers from the station orbit region to `run.phase3.parking_altitude_km`, waits until the vehicle is below the target on the target R-bar, then performs a final injection to `run.phase3.entry_interface_altitude_km` using the configured FPA geometry. The final injection delta-V and propellant are included by default; disabling `run.phase3.charge_final_reentry_fuel` is an explicitly nonphysical mass-ledger sensitivity option.

Both modes still use the same post-run check of the actual flight-path angle near the configured entry-interface altitude.

If orbital atmospheric drag is enabled for Phase 3, `HOHMANN` switches to the
drag-aware single-retrograde-burn design loaded from
`run.phase3.drag_deorbit_design.file`. Set
`run.phase3.drag_deorbit_design.mode = "OFF"` to force the older dragless conic
calculation even with orbital drag enabled.

For explicitly generated finite-burn deorbit JSONs, MATLAB recomputes burn
duration and propellant from the actual Phase 3 mass after
rendezvous/proximity operations. The Python design mass is still recorded in
the JSON for traceability.

### Phase 4: Atmospheric Entry

After Phase 3 reaches the configured entry interface, `Reentry_Propagator.m` propagates a separated re-entry vehicle through a co-rotating ISA76 atmosphere. The translational state remains ECI, but drag and lift use atmosphere-relative velocity, which is equivalent to an ECEF-relative aerodynamic velocity model.

SPACEPLANE resolves AoA from the configured speed schedule at every dynamics
evaluation unless `run.reentry.aoa_deg` supplies a constant override. CAPSULE
uses a constant trim surrogate. `run.reentry.bank_angle_deg` rotates the lift
direction about the atmosphere-relative velocity vector and remains open-loop
constant. The flight-path angle is not commanded during Phase 4; it is computed
from the propagated position and velocity state and evolves naturally under
gravity, drag, and lift.

The pre-Phase-3 chaser state is propagated separately as an orbiting relay. The code checks whether the re-entry vehicle and relay chaser maintain geometric line of sight by testing Earth-limb clearance along the connecting segment. It also plots heat flux, total heat load, dynamic pressure, aero g-load, and LOS clearance/elevation.

## Main Files

### `Main_Mission_Simulator.m`

Orchestrates the complete mission:

- initializes spacecraft states
- runs Phase 1 / 2 / 3 in sequence
- collects delta-V and propellant budgets
- generates mission plots

### `Mission_Config.m`

Contains mission constants and tuning parameters such as:

- Earth constants (`mu`, `Re`, `J2`, `g0`)
- vehicle mass and inertia
- mission altitudes
- propulsion system parameters
- optional atmospheric-drag settings: model, Cd, reference area, and atmosphere co-rotation
- re-entry vehicle shape definitions and entry settings
- sensor / thrust uncertainty settings
- nominal capture-point settings

### `Phasing_Propagator.m`

Handles impulsive orbital transfer logic for:

- Hohmann transfer
- custom impulse phase targeting
- target co-propagation
- history logging of chaser/target states
- J2-aware wait-time targeting for capture-point arrival

### `Env_EOM.m`

Defines the equations of motion used for translational and rotational propagation.

### Python Helpers

`J2PolarHohmann.py` and `J2PolarHohmannShooting.py` provide a SciPy-based side workflow for studying and tuning J2 polar Hohmann-like rendezvous parameters. They are not required for a normal MATLAB run, but they explain where the active `CUSTOM_IMPULSE` values can come from. The Python workflow can also export atmospheric-drag settings into the MATLAB JSON config.

## Model Assumptions

This repository is best interpreted as a research / educational simulator rather than a flight-grade GNC tool.

Important assumptions include:

- Earth gravity + J2, with optional atmospheric drag
- no SRP, third-body gravity, or full Earth-fixed frame dynamics beyond the simple co-rotating atmosphere used by the drag and entry models
- drag-aware `HOHMANN` deorbit design is a single retrograde burn followed by numerical drag propagation, not a full entry-guidance, landing-target, or high-dimensional trajectory optimizer
- SPACEPLANE uses the Zhang polynomial `CL/CD` with project-assumed area and mass; CAPSULE uses constant `Cd/L/D` because the source aerodynamic database is unavailable
- SPACEPLANE AoA is speed-scheduled and CAPSULE trim is constant; bank is open-loop constant and no attitude propagation, trim solve, or bank guidance law is included yet
- simplified impulsive burns in the phasing and proximity stages
- simplified thrust and sensor error models
- simplified capture / berthing logic
- no high-fidelity actuator, navigation filter, abort logic, or docking contact model

A more detailed discussion is provided in `docs/ASSUMPTIONS_AND_LIMITATIONS.md`.

For a Korean technical report with theory, usage, limitations, and references, see `docs/RENDEZVOUS_REENTRY_REPORT_KR.md`.

## Current Known Limitations

At the current development stage, the code is useful for mission logic prototyping and relative-motion debugging, but not yet for high-fidelity mission reconstruction.

Examples of current limitations:

- Hohmann targeting is still a best-fit capture-point search, not a strict boundary-value solution.
- The transfer is not yet optimized for exact arrival relative velocity.
- The active proximity-operations model is waypoint-impulsive and not yet tied to a realistic navigation-estimation pipeline.
- Orbit initialization and inclination usage should be reviewed for strict physical consistency.
- Mission phases are logically connected, but still need refinement for a fully validated end-to-end scenario.

## Suggested Development Roadmap

1. Exact capture-point targeting using shooting / differential correction
2. Relative-velocity matching at capture
3. Navigation filter for noisy state estimation
4. More realistic perturbation set: drag, SRP, higher-order gravity if needed
5. Monte Carlo campaign support
6. Automated regression test scripts
7. Batch mission-case runner for parameter sweeps

## Citation and Research Attribution

If you use this simulator in research, cite the repository using
[`CITATION.cff`](CITATION.cff) and cite the Zhang and/or Saito source article
for every paper-derived model component used in the analysis. Detailed reuse
boundaries, article notices, DOI links, provenance classes, and the
non-endorsement statement are recorded in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

The source-paper PDFs, article figures, substantial article text, proprietary
data, and unpublished inputs are not redistributed by this repository. Local
reference PDFs under `tmp/pdfs` are intentionally excluded from Git tracking.

## License

Original source code and documentation in this repository are licensed under
the [BSD 3-Clause License](LICENSE), unless a file explicitly states otherwise.
The license does not apply to the cited papers or other third-party material;
see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
