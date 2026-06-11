# Unmanned Rendezvous Mission Simulator

MATLAB-based end-to-end rendezvous mission simulator for an unmanned chaser spacecraft operating in Earth orbit. The current codebase combines:

- **Phase 1**: 3-DOF phasing / homing with impulsive maneuvers and J2-aware orbital propagation
- **Phase 2**: LVLH waypoint-impulse proximity operations with cycloidal drift, R-bar hops, nonlinear J2 propagation, and mass bookkeeping
- **Phase 3**: selectable 3-DOF de-orbit / re-entry modes, including direct FPA-targeted descent and a 200 km R-bar-aligned re-entry setup

The project is designed as a mission-level simulation framework rather than a single guidance-law demo. It is useful for studying how orbit transfer logic, target-relative geometry, J2 perturbation, mass depletion, and proximity-operation sequencing interact in one continuous workflow.

## Current Scope

This repository currently models:

- Earth central gravity + J2 perturbation
- A target spacecraft and a chaser spacecraft
- Hohmann-based phase planning with target co-propagation
- A J2-aware wait-time search before departure so the transfer arrives closer to a desired LVLH capture point
- Custom impulse phasing driven by externally tuned phase angle, delta-V, and gamma parameters
- LVLH waypoint-impulse proximity operations with cleanup, hold trims, cycloidal drift, R-bar hops, braking impulses, and mass depletion
- selectable Phase 3 re-entry logic through `reentry_mode`
- Simple thrust uncertainty / noise injection in selected phasing modes
- Mission-level delta-V and propellant budget tracking
- Visualization of trajectory, altitude history, proximity trajectory, and mass depletion

## Repository Structure

```text
.
|-- Main_Mission_Simulator.m      % Main entry point for the full mission
|-- Mission_Config.m              % Physical constants, vehicle data, mission parameters
|-- Phasing_Propagator.m          % 3-DOF phasing / Hohmann / Lambert / custom impulse logic
|-- Env_EOM.m                     % Environment and 3-DOF/6-DOF equations of motion
|-- J2PolarHohmann.py             % Python J2 polar Hohmann propagation study
|-- J2PolarHohmannShooting.py     % Python shooting / constrained optimization helper
|-- docs/
|   |-- ARCHITECTURE.md
|   `-- ASSUMPTIONS_AND_LIMITATIONS.md
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
3. Run:

```matlab
Main_Mission_Simulator
```

The active MATLAB path is intended to use standard MATLAB numerical functionality. The Python helper scripts require the packages listed in `requirements.txt`.

Phase 3 is selected by `reentry_mode` in `Main_Mission_Simulator.m`. You can override it without editing the file:

```matlab
setenv('RENDEZVOUS_PHASE3_MODE','R_BAR_200_FPA')
Main_Mission_Simulator
```

For `R_BAR_200_FPA`, the final 200 km to 120 km injection fuel update is controlled independently. The default in the script is `off`; to include that burn in the propellant and remaining-mass budget, run:

```matlab
setenv('RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL','on')
Main_Mission_Simulator
```

## Simulation Flow

### Phase 1: Phasing & Homing

`Phasing_Propagator.m` performs the pre-rendezvous orbital transfer.

The active root script currently uses `CUSTOM_IMPULSE` mode, where the departure phase angle, burn magnitude, and burn flight-path tilt are supplied through `custom_params`. Those values can be generated or refined using the Python shooting helpers.

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

Older continuous-force `GNC_Controller.m` implementations are preserved under `legacy/`, but they are not called by the current root-level mission script.

### Phase 3: Re-entry / Descent

Phase 3 is selected with `reentry_mode` in `Main_Mission_Simulator.m`.

Current modes:

- `HOHMANN`: preserves the previous direct FPA-targeted Hohmann-style descent. It computes the radius needed to pass near 120 km at the configured 4 deg flight-path angle, then sends that target radius to `Phasing_Propagator.m`.
- `R_BAR_200_FPA`: first lowers from the station orbit region to a 200 km parking orbit, waits until the vehicle is below the target on the target R-bar, then performs a final 200 km to 120 km injection using the same FPA geometry. The final injection delta-V is always reported; its propellant can be included or excluded with `charge_final_reentry_fuel` or `RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL`.

Both modes still use the same post-run check of the actual flight-path angle near 120 km altitude.

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
- sensor / thrust uncertainty settings
- nominal capture-point settings

### `Phasing_Propagator.m`

Handles impulsive orbital transfer logic for:

- Hohmann transfer
- Lambert transfer
- custom impulse phase targeting
- target co-propagation
- history logging of chaser/target states
- J2-aware wait-time targeting for capture-point arrival

### `Env_EOM.m`

Defines the equations of motion used for translational and rotational propagation.

### Python Helpers

`J2PolarHohmann.py` and `J2PolarHohmannShooting.py` provide a SciPy-based side workflow for studying and tuning J2 polar Hohmann-like rendezvous parameters. They are not required for a normal MATLAB run, but they explain where the active `CUSTOM_IMPULSE` values can come from.

## Model Assumptions

This repository is best interpreted as a research / educational simulator rather than a flight-grade GNC tool.

Important assumptions include:

- Earth gravity + J2 only
- no atmospheric drag, SRP, third-body gravity, or full Earth-fixed frame effects in the orbital propagator
- simplified impulsive burns in the phasing and proximity stages
- simplified thrust and sensor error models
- simplified capture / berthing logic
- no high-fidelity actuator, navigation filter, abort logic, or docking contact model

A more detailed discussion is provided in `docs/ASSUMPTIONS_AND_LIMITATIONS.md`.

## Current Known Limitations

At the current development stage, the code is useful for mission logic prototyping and relative-motion debugging, but not yet for high-fidelity mission reconstruction.

Examples of current limitations:

- Hohmann targeting is still a best-fit capture-point search, not a strict boundary-value solution.
- The transfer is not yet optimized for exact arrival relative velocity.
- The active proximity-operations model is waypoint-impulsive and not yet tied to a realistic navigation-estimation pipeline.
- Orbit initialization and inclination usage should be reviewed for strict physical consistency.
- Mission phases are logically connected, but still need refinement for a fully validated end-to-end scenario.

## Suggested Development Roadmap

1. Exact capture-point targeting using Lambert / shooting / differential correction
2. Relative-velocity matching at capture
3. Navigation filter for noisy state estimation
4. More realistic perturbation set: drag, SRP, higher-order gravity if needed
5. Monte Carlo campaign support
6. Automated regression test scripts
7. Batch mission-case runner for parameter sweeps

## License

No license file is included yet. Before making the repository public, choose a license explicitly, such as MIT, BSD-3-Clause, or GPL-3.0 depending on the intended reuse model.
