# Architecture Overview

## 1. Top-Level Execution

The simulation starts in `Main_Mission_Simulator.m`.

Its responsibilities are:

1. load mission constants from `Mission_Config.m`
2. initialize chaser and target states
3. run the mission in four major phases
4. accumulate delta-V / fuel usage
5. visualize the mission history

## 2. Module Responsibilities

### `Mission_Config.m`

Provides a single configuration structure `sys` containing:

- environmental constants
- spacecraft mass and inertia
- propulsion parameters
- GNC-related gains retained for compatibility and future controller work
- uncertainty settings
- mission altitude targets
- capture-point settings

This file is the main place to edit scenario assumptions.

### `Phasing_Propagator.m`

Responsible for mission phases dominated by orbital transfer rather than close-range feedback control.

It currently supports:

- impulsive Hohmann departure and circularization
- custom phase targeting with impulsive execution by default and finite-burn execution as an explicit study option
- preliminary Multi-Hohmann transfer code retained for future burn-duration / thermal-constraint studies, but not part of the default workflow
- target co-propagation during transfer
- J2-aware wait-time scanning before departure
- history logging for both chaser and target

Internally, the Hohmann mode behaves more like a target-aware transfer planner than a simple altitude-raising maneuver.

### Phase 2 Helpers in `Main_Mission_Simulator.m`

The active root-level Phase 2 implementation is embedded in `Main_Mission_Simulator.m` as local helper functions instead of using the legacy `GNC_Controller.m`.

Major internal steps:

1. compute target-centered LVLH frame
2. compute relative position and velocity with LVLH frame rotation
3. optionally clean up the Phase 1 terminal state to the S2 -V-bar hold point
4. trim residual relative velocity at S2
5. apply a single V-bar impulse and coast through a cycloidal S2-to-S3 drift
6. detect S3 at the first V-bar crossing
7. perform CW/Hill targeted R-bar hops from S3 to S4
8. apply braking impulses at hold points
9. propagate the nonlinear chaser/target dynamics with `Env_EOM.m`

Older force-based `GNC_Controller.m` files are preserved in `legacy/`, but they are not part of the active root-level execution path.

### `Env_EOM.m`

Provides the physical propagation model.

Depending on mode, it supports:

- translational 3-DOF propagation under gravity + J2
- translational + rotational 6-DOF propagation
- thrust-driven mass depletion

### `Reentry_Propagator.m` and `+reentry_core`

`Reentry_Propagator.m` is the atmospheric-entry orchestration layer. It owns
vehicle/scenario selection, relay propagation, the Zhang blackout-zone latch,
history logging, path-constraint aggregation, and terminal-event policy.

The numerically reusable part is in the `+reentry_core` MATLAB package:

- seven-state `[r_ECI; v_ECI; mass]` dynamics and RK4 stepping
- central/J2 gravity and co-rotating-atmosphere relative velocity
- spaceplane/capsule aerodynamic evaluation and command resolution
- altitude and speed event refinement
- FPA, Earth-limb LOS, RAAP, and antenna geometry

Keeping the physical kernel independent of mission and paper policy lets the
integrated mission and standalone studies exercise one implementation. The
public `Reentry_Propagator` signature and history/summary interface remain
unchanged.

### `+paperstudies`

The `+paperstudies/+zhang` and `+paperstudies/+saito` packages are standalone
paper-condition harnesses. They store published tables and equations,
machine-readable provenance and blockers, deterministic self-tests, and
optional surrogate adapters to the shared physical kernel. They do not change
the normal mission configuration and do not present omitted paper inputs as
known quantities.

`Run_Paper_Reproduction_Suite.m` is the common entry point. See
`docs/PAPER_REPRODUCTION_FRAMEWORK.md` for the evidence levels and the boundary
between structural reproduction and full numerical trajectory reproduction.

### Python Helpers

`J2PolarHohmann.py` and `J2PolarHohmannShooting.py` form a separate SciPy workflow for studying J2 polar Hohmann-like rendezvous behavior and tuning the active `CUSTOM_IMPULSE` parameters used by the MATLAB script. The shooting script can export a MATLAB-readable JSON configuration.

## 3. Data Flow

### Phase 1

`Main_Mission_Simulator`
to `Phasing_Propagator`
to updated chaser state, updated target state, mission history, delta-V, fuel

### Phase 2

`Main_Mission_Simulator`
to local LVLH/CW/impulse helper functions
to `Env_EOM`
to updated chaser and target states, proximity history, delta-V, fuel

### Phase 3

`Main_Mission_Simulator`
to Phase 3 local helpers / `Phasing_Propagator`
to 120 km entry-interface state and budget update

### Phase 4

`Main_Mission_Simulator`
to `Reentry_Propagator`
to atmospheric-entry history, heating/dynamic-pressure/g-load diagnostics, and chaser-to-entry-vehicle LOS summary

### Standalone paper studies

`Run_Paper_Reproduction_Suite`
to `paperstudies.zhang.run` and `paperstudies.saito.run`
to published-condition regression reports
to optional `Reentry_Propagator` surrogate forward histories

## 4. Frames Used

The project uses both:

- ECI-like inertial coordinates for orbital propagation
- LVLH coordinates for proximity guidance and interpretation of relative motion
- atmosphere-relative velocity with a co-rotating atmosphere for drag and entry aerodynamics

Important practical note: the physical meaning of the final approach depends strongly on how the LVLH frame is built from the target state and how relative velocity is defined with frame rotation included.

## 5. Why the Repository Is Structured This Way

The current structure reflects a useful separation of concerns:

- orbital mission design logic
- close-range waypoint / impulse sequencing
- equations of motion / physics model
- scenario configuration
- mission orchestration and plotting
- Python-side parameter studies

That separation makes future upgrades easier, such as:

- replacing current capture targeting with a validated boundary-value / shooting correction method
- restoring or replacing the proximity controller with CW / HCW / MPC / LQR logic
- introducing estimation and navigation filters
- adding batch Monte Carlo runners
