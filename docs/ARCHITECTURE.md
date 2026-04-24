# Architecture Overview

## 1. Top-Level Execution

The simulation starts in `Main_Mission_Simulator.m`.

Its responsibilities are:

1. load all mission constants from `Mission_Config.m`
2. initialize chaser and target states
3. run the mission in three major phases
4. accumulate delta-V / fuel usage
5. visualize the mission history

---

## 2. Module Responsibilities

### `Mission_Config.m`
Provides a single configuration structure `sys` containing:

- environmental constants
- spacecraft mass and inertia
- propulsion parameters
- GNC gains
- uncertainty settings
- mission altitude targets
- capture-point settings

This file is the main place to edit scenario assumptions.

### `Phasing_Propagator.m`
Responsible for mission phases that are dominated by orbital transfer rather than close-range feedback control.

It currently supports:

- impulsive Hohmann departure and circularization
- Lambert-style transfer option
- target co-propagation during transfer
- J2-aware wait-time scanning before departure
- history logging for both chaser and target

Internally, the Hohmann mode now behaves more like a **target-aware transfer planner** than a simple altitude-raising maneuver.

### `GNC_Controller.m`
Responsible for close-range relative-motion control.

Major internal steps:

1. compute target-centered LVLH frame
2. compute relative position and velocity
3. estimate environmental acceleration terms in LVLH
4. generate guidance commands for R-bar closing / final berthing
5. map LVLH force demand to body-frame thrust command
6. damp angular rates using attitude control torque

### `Env_EOM.m`
Provides the physical propagation model.

Depending on mode, it supports:

- translational 3-DOF propagation under gravity + J2
- translational + rotational 6-DOF propagation
- thrust-driven mass depletion

---

## 3. Data Flow

### Phase 1
`Main_Mission_Simulator`
→ `Phasing_Propagator`
→ updated chaser state, updated target state, mission history, delta-V, fuel

### Phase 2
`Main_Mission_Simulator`
→ `GNC_Controller`
→ commanded force / torque
→ `Env_EOM`
→ updated chaser and target states

### Phase 3
`Main_Mission_Simulator`
→ `Phasing_Propagator`
→ descent / de-orbit style propagation and budget update

---

## 4. Frames Used

The project uses both:

- **ECI-like inertial coordinates** for orbital propagation
- **LVLH coordinates** for proximity guidance and interpretation of relative motion

Important practical note:
The physical meaning of the final approach depends strongly on how the LVLH frame is built from the target state and how relative velocity is defined with frame rotation included.

---

## 5. Why the Repository Is Structured This Way

The current structure reflects a useful separation of concerns:

- orbital mission design logic
- close-range guidance and control logic
- equations of motion / physics model
- scenario configuration
- mission orchestration and plotting

That separation makes future upgrades easier, such as:

- replacing Hohmann targeting with Lambert + correction
- swapping the controller for CW / HCW / MPC / LQR logic
- introducing estimation and navigation filters
- adding batch Monte Carlo runners
