# Assumptions and Limitations

This simulator is physically motivated, but it is **not** a flight-grade rendezvous analysis tool in its current form.

## 1. Orbital Environment

Currently included:

- central Earth gravity
- J2 perturbation
- optional atmospheric drag with an ISA76-style standard-atmosphere density model
- 120 km atmospheric-entry propagation with simplified drag/lift, heat-flux, and line-of-sight diagnostics

Currently neglected:

- solar radiation pressure
- third-body perturbations
- full Earth orientation / Earth-fixed frame effects beyond simple co-rotating atmosphere velocity for drag
- higher-order geopotential terms

This means the code is appropriate for **conceptual mission studies and controller debugging**, but not for high-precision orbit determination or operational conjunction analysis.

---

## 2. Maneuver Modeling

The phasing stage supports simplified impulsive maneuver logic and an optional finite-burn execution model for the custom phased maneuver.

This is useful for:

- first-order delta-V budgeting
- transfer geometry studies
- mission-sequence prototyping
- early finite-burn sensitivity checks using a fixed thrust/Isp model

But it is limited for:

- finite-burn attitude-coupled maneuver analysis
- exact burn reconstruction
- high-fidelity thruster execution modeling

The atmospheric-entry stage is also simplified:

- the re-entry vehicle keeps a fixed AoA and bank angle during entry
- Cd is assumed constant because the supplied PPT provides geometry and L/D trends, not full Cd/Cl tables
- L/D curves are approximate digitizations from slide images
- heat flux uses a Sutton-Graves stagnation-point estimate
- communication is checked as geometric Earth-limb line of sight only; plasma blackout physics is not modeled yet

---

## 3. Rendezvous Interpretation

The current Hohmann-based phase planner is **target-aware**, but it still works as a best-fit timing search rather than a strict exact rendezvous solver.

So the code can meaningfully answer questions like:

- when should the chaser depart to get closer to a desired capture geometry?
- how does J2-aware co-propagation affect transfer timing?
- how does the mission-level delta-V change?

But it should not yet be treated as an authoritative tool for:

- exact capture-state matching
- exact relative-velocity nulling at arrival
- certified rendezvous corridor analysis

---

## 4. Proximity Operations

The close-range controller uses LVLH relative motion, environmental feedforward, and PD-style logic.

This is useful for:

- controller prototyping
- relative-motion visualization
- debugging approach trajectories

But the following are still simplified or missing:

- state estimation filter
- measurement scheduling / asynchronous sensors
- realistic actuator dynamics and saturation logic
- fault detection / abort logic
- docking contact dynamics

---

## 5. State Consistency

A few modeling choices should be reviewed carefully when the project is upgraded:

- consistency between declared inclination and actual initialized orbital states
- exact interpretation of the capture point in LVLH vs inertial coordinates
- continuity between mission phases under realistic arrival constraints
- separation between debug resets and physically continuous propagation

---

## 6. Recommended Use Cases

Recommended:

- educational orbital mechanics projects
- early-stage rendezvous mission design logic
- J2-aware phasing experiments
- controller debugging and visualization
- mission delta-V / mass bookkeeping studies

Not recommended yet:

- operational mission planning
- precise delta-V inversion from real-world data
- collision risk assessment
- high-fidelity spacecraft certification work
