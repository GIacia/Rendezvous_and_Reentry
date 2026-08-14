# Assumptions and Limitations

This simulator is physically motivated, but it is **not** a flight-grade rendezvous analysis tool in its current form.

## 1. Orbital Environment

Currently included:

- central Earth gravity
- J2 perturbation
- optional atmospheric drag with an ISA76-style standard-atmosphere density model
- configurable entry-interface atmospheric propagation with simplified drag/lift, heat-flux, and line-of-sight diagnostics

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
- drag-aware LEO deorbit checks using either an impulsive or finite-duration single retrograde burn followed by numerical atmospheric-drag propagation

But it is limited for:

- finite-burn attitude-coupled maneuver analysis
- exact burn reconstruction
- high-fidelity thruster execution modeling
- full deorbit guidance, landing footprint targeting, or dispersions under uncertain atmosphere
- high-dimensional trajectory optimization; use this simulator for validation, sensitivity checks, and Monte Carlo after an optimizer has produced a candidate trajectory

The atmospheric-entry stage is also simplified:

- SPACEPLANE uses a paper-derived speed-scheduled AoA and polynomial `CL/CD`;
  CAPSULE uses a constant `Cd=1.3`, `L/D=0.25` surrogate because the paper's
  proprietary `CA/CN(Mach,AoA)` database and trim AoA are unavailable
- bank remains an open-loop constant; no trajectory-level optimizer,
  predictor-corrector guidance, or attitude/RCS controller is present
- built-in spaceplane area/nose radius and all capsule nose-radius data are
  project assumptions, so aerodynamic and aerothermal model-form uncertainty
  dominates fine numerical accuracy
- heat flux is a Sutton-Graves surrogate; it is not quantitatively comparable
  to either paper's published heat law, so paper heat compliance is reported
  as unevaluated
- communication checks RAAP cone, configured geometric range, and Earth
  occultation; it does not model RF link margin, antenna gain pattern, plasma
  attenuation, or service availability
- the analytic paper TDRS accumulates all mission-phase elapsed time, but the
  mission has no UTC epoch; its default Greenwich-meridian ECI angle is the
  explicit zero-degree project convention and must be replaced by a
  GMST-derived value for an epoch-tied RAAP/LOS analysis

The standalone Zhang/Saito study packages improve traceability but do not
remove those physical limitations. Their regression tests establish software
equivalence, source transcription, public-equation algebra, coordinate
round-trips, and uncertainty-grid construction. They are not external
validation of the assumed atmosphere, aerodynamics, heating, or guidance. Full
trajectory-level reproduction remains unavailable until Zhang's vehicle/OCP
inputs and Saito's proprietary aerodynamic and omitted guidance/controller
inputs are supplied. See `PAPER_REPRODUCTION_FRAMEWORK.md` for the exact claim
boundary.

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
