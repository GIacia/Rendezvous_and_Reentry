# Unmanned Rendezvous Mission Simulator

MATLAB-based end-to-end rendezvous mission simulator for an unmanned chaser spacecraft operating in Earth orbit. The current codebase combines:

- **Phase 1**: 3-DOF phasing / homing with impulsive maneuvers and J2-aware orbital propagation
- **Phase 2**: LVLH waypoint-impulse proximity operations with cycloidal drift, R-bar hops, nonlinear J2 propagation, and mass bookkeeping
- **Phase 3**: selectable 3-DOF de-orbit / re-entry modes, including direct FPA-targeted descent and a 200 km R-bar-aligned re-entry setup
- **Phase 4**: atmospheric entry from the 120 km interface, with re-entry vehicle shape selection, heat-rate diagnostics, and chaser-to-entry-vehicle line-of-sight checks

The project is designed as a mission-level simulation framework rather than a single guidance-law demo. It is useful for studying how orbit transfer logic, target-relative geometry, J2 perturbation, mass depletion, proximity-operation sequencing, and entry diagnostics interact in one connected workflow.

## Current Scope

This repository currently models:

- Earth central gravity + J2 perturbation
- Optional atmospheric drag using an ISA76-style standard-atmosphere density model
- A target spacecraft and a chaser spacecraft
- Hohmann-based phase planning with target co-propagation
- A J2-aware wait-time search before departure so the transfer arrives closer to a desired LVLH capture point
- Multi-Hohmann phasing option for splitting large maneuvers across multiple thermally limited burns
- Custom phased maneuver logic driven by externally tuned phase angle, delta-V, and gamma parameters, with selectable impulsive or finite-burn execution
- LVLH waypoint-impulse proximity operations with cleanup, hold trims, cycloidal drift, R-bar hops, braking impulses, and mass depletion
- selectable Phase 3 re-entry logic through `reentry_mode`
- Simple thrust uncertainty / noise injection in selected phasing modes
- Mission-level delta-V and propellant budget tracking
- Visualization of trajectory, altitude history, proximity trajectory, mass depletion, atmospheric-entry heat flux, dynamic pressure, g-load, and line-of-sight margin

## Repository Structure

```text
.
|-- Main_Mission_Simulator.m      % Main entry point for the full mission
|-- Mission_Config.m              % Physical constants, vehicle data, mission parameters
|-- Phasing_Propagator.m          % 3-DOF phasing, Hohmann, Multi-Hohmann, custom impulse logic
|-- Env_EOM.m                     % Environment and 3-DOF/6-DOF equations of motion
|-- Atmospheric_Drag_Acceleration.m % Shared MATLAB atmospheric-drag model
|-- Standard_Atmosphere_Density.m % Shared ISA76-style atmosphere helper
|-- Reentry_Propagator.m          % 120 km atmospheric-entry propagation and LOS diagnostics
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

Python shooting results can be exported directly into a MATLAB-readable mission config:

```bash
python J2PolarHohmannShooting.py --no-plot --matlab-config-out configs/latest_python_solution.json
```

Atmospheric drag is off by default. To optimize with the ISA76 drag option and export the same environment settings to MATLAB:

```bash
python J2PolarHohmannShooting.py --no-plot --atmospheric-drag isa76 --matlab-config-out configs/latest_python_solution.json
```

Then load that config in MATLAB:

```matlab
setenv('RENDEZVOUS_CONFIG_JSON','configs/latest_python_solution.json')
Main_Mission_Simulator
```

The JSON file overrides only the fields it contains, so `Mission_Config.m` remains the default source of truth.

For a quick MATLAB-only drag run without editing files:

```matlab
setenv('RENDEZVOUS_ATMOSPHERIC_DRAG','ISA76')
Main_Mission_Simulator
```

Clear it with `setenv('RENDEZVOUS_ATMOSPHERIC_DRAG','OFF')` or by clearing the environment variable.

The atmospheric entry vehicle shape is selected independently. Available values are `COMPROMISE`, `HEATLOAD_MIN`, `PAYLOAD_MAX`, and `TPS_MIN`:

```matlab
setenv('RENDEZVOUS_REENTRY_SHAPE','TPS_MIN')
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

The active root script currently uses `CUSTOM_IMPULSE` mode, where the departure phase angle, burn magnitude, and burn flight-path tilt are supplied through `custom_params`. The maneuver execution model is selected separately with `custom_params.burn_model`: `"IMPULSIVE"` preserves the legacy instantaneous delta-V behavior, while `"FINITE_BURN"` applies the same delta-V direction over a short high-thrust burn using the configured thrust, Isp, and mass depletion. Those values can be generated or refined using the Python shooting helpers.

To switch the MATLAB custom maneuver to finite-burn mode without editing the script:

```matlab
setenv('RENDEZVOUS_BURN_MODEL','FINITE_BURN')
Main_Mission_Simulator
```

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

Phase 3 is selected with `reentry_mode` in `Main_Mission_Simulator.m`.

Current modes:

- `HOHMANN`: performs a direct FPA-targeted de-orbit injection and stops at the 120 km entry interface. It no longer performs a nonphysical circularization after crossing the interface.
- `R_BAR_200_FPA`: first lowers from the station orbit region to a 200 km parking orbit, waits until the vehicle is below the target on the target R-bar, then performs a final 200 km to 120 km injection using the same FPA geometry. The final injection delta-V is always reported; its propellant can be included or excluded with `charge_final_reentry_fuel` or `RENDEZVOUS_CHARGE_FINAL_REENTRY_FUEL`.

Both modes still use the same post-run check of the actual flight-path angle near 120 km altitude.

### Phase 4: Atmospheric Entry

After Phase 3 reaches the 120 km interface, `Reentry_Propagator.m` propagates a separated re-entry vehicle through a co-rotating ISA76 atmosphere. The translational state remains ECI, but drag and lift use atmosphere-relative velocity, which is equivalent to an ECEF-relative aerodynamic velocity model.

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
- Multi-Hohmann transfer
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
- re-entry vehicle Cd is currently an assumed constant; only L/D trends and geometry are taken from the supplied shape PPT
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

## License

No license file is included yet. Before making the repository public, choose a license explicitly, such as MIT, BSD-3-Clause, or GPL-3.0 depending on the intended reuse model.
