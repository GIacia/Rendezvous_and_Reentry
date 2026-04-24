# Unmanned Rendezvous Mission 6-DOF Simulator

MATLAB-based end-to-end rendezvous mission simulator for an unmanned chaser spacecraft operating in Earth orbit. The current codebase combines:

- **Phase 1**: 3-DOF phasing / homing with impulsive maneuvers and **J2-aware orbital propagation**
- **Phase 2**: 6-DOF proximity operations with **LVLH relative navigation**, **R-bar approach guidance**, and **attitude dynamics**
- **Phase 3**: 3-DOF de-orbit / re-entry descent to a lower orbit

The project is designed as a mission-level simulation framework rather than a single guidance law demo. It is especially useful for studying how **orbit transfer logic**, **target-relative geometry**, **J2 perturbation**, **mass depletion**, and **proximity-operation control** interact in one continuous workflow.

---

## Current Scope

This repository currently models:

- Earth central gravity + **J2 perturbation**
- A target spacecraft and a chaser spacecraft
- Hohmann-based phase planning with **target co-propagation**
- A **J2-aware wait-time search** before departure so the transfer arrives closer to a desired LVLH capture point
- 6-DOF chaser translational + rotational dynamics during proximity operations
- Simple thrust uncertainty / noise injection
- Mission-level **delta-V and propellant budget tracking**
- Visualization of trajectory, altitude history, and mass depletion

---

## Repository Structure

```text
.
├── Main_Mission_Simulator.m      % Main entry point for the full mission
├── Mission_Config.m              % Physical constants, vehicle data, mission parameters
├── Phasing_Propagator.m          % 3-DOF phasing / Hohmann / Lambert propagation logic
├── GNC_Controller.m              % LVLH relative guidance and control for proximity ops
├── Env_EOM.m                     % Environment and 6-DOF equations of motion
├── docs/
│   ├── ARCHITECTURE.md
│   └── ASSUMPTIONS_AND_LIMITATIONS.md
├── results/
│   └── README.md
├── examples/
│   └── README.md
└── .gitignore
```

The current `.m` files are intentionally kept at the repository root so that the existing MATLAB execution flow does **not** break.

---

## How to Run

1. Open the repository folder in MATLAB.
2. Make sure the current folder is the repository root.
3. Run:

```matlab
Main_Mission_Simulator
```

No external toolbox dependencies are intentionally assumed beyond standard MATLAB numerical functionality.

---

## Simulation Flow

### Phase 1 — Phasing & Homing

`Phasing_Propagator.m` performs the pre-rendezvous orbital transfer.

In the current Hohmann mode, the logic is:

1. Co-propagate the **chaser** and **target** under J2
2. Search for a departure wait time that best aligns the transfer arrival with a desired LVLH capture point
3. Execute the first Hohmann impulse
4. Numerically propagate the transfer arc
5. Circularize at the destination altitude

This is more physically meaningful than immediately burning at `t = 0`, because it includes target motion during the waiting phase.

### Phase 2 — Proximity Operations

`GNC_Controller.m` computes the commanded force/torque for R-bar closing and final berthing in the target-centered LVLH frame.

The controller currently includes:

- LVLH relative position / velocity computation
- Environmental feedforward compensation
- PD-like translational guidance
- Mode switching between `R-BAR_CLOSING` and `FINAL_BERTHING`
- Simple body torque damping for attitude stabilization

### Phase 3 — Re-entry / Descent

The same phasing propagator is reused to lower the orbit from the target-altitude region toward the re-entry interface altitude.

---

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

- **Hohmann transfer**
- **Lambert transfer**
- target co-propagation
- history logging of chaser/target states
- J2-aware wait-time targeting for capture-point arrival

### `GNC_Controller.m`
Implements the proximity-operations guidance and control logic in LVLH.

### `Env_EOM.m`
Defines the equations of motion used for translational and rotational propagation.

---

## Model Assumptions

This repository is best interpreted as a **research / educational simulator** rather than a flight-grade GNC tool.

Important assumptions include:

- Earth gravity + J2 only
- no atmospheric drag, SRP, third-body gravity, or Earth rotation effects in the orbital propagator
- simplified impulsive burns in the phasing stage
- simplified thrust and sensor error models
- simplified capture / berthing logic
- no high-fidelity actuator or navigation filter model

A more detailed discussion is provided in [`docs/ASSUMPTIONS_AND_LIMITATIONS.md`](docs/ASSUMPTIONS_AND_LIMITATIONS.md).

---

## Current Known Limitations

At the current development stage, the code is useful for **mission logic prototyping** and **relative-motion debugging**, but not yet for high-fidelity mission reconstruction.

Examples of current limitations:

- Hohmann targeting is still a **best-fit capture-point search**, not a strict boundary-value solution
- the transfer is not yet optimized for exact arrival relative velocity
- the controller is not yet tied to a realistic navigation-estimation pipeline
- orbit initialization and inclination usage should be reviewed for strict physical consistency
- mission phases are logically connected, but still need refinement for a fully validated end-to-end scenario

---

## Suggested Development Roadmap

1. **Exact capture-point targeting** using Lambert / shooting / differential correction
2. **Relative-velocity matching** at capture
3. **Navigation filter** for noisy state estimation
4. **More realistic perturbation set**: drag, SRP, higher-order gravity if needed
5. **Monte Carlo campaign support**
6. **Automated regression test scripts**
7. **Batch mission-case runner** for parameter sweeps

---

## Recommended GitHub Usage

Good commits for this repository are usually organized by one physical topic at a time, for example:

- `fix phase-1 target co-propagation`
- `add J2-aware wait-time search`
- `refactor LVLH relative-state logging`
- `tune final berthing controller`

This makes the evolution of the mission logic much easier to trace.

---

## License

No license file is included yet.

Before making the repository public, it is recommended to choose a license explicitly, for example:

- **MIT License** for maximum reuse
- **BSD-3-Clause** for permissive academic/project use
- **GPL-3.0** if derivative open-source sharing is important

---

## Author Notes

This repository is being actively shaped toward an **end-to-end unmanned rendezvous mission simulator**. The code already contains several physically meaningful building blocks, and the current documentation is written to make that development path easier to understand and continue.
