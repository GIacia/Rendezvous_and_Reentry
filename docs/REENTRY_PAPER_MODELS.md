# Paper-Based Re-entry Modes

The atmospheric-entry propagator supports two selectable reduced-order
vehicle modes:

```matlab
run.reentry.vehicle_mode = "SPACEPLANE";
% or
run.reentry.vehicle_mode = "CAPSULE";
```

These are deterministic trajectory and constraint-evaluation models. They do
not yet include the optimal-control or closed-loop bank guidance algorithms
from the source papers.

The paper-style 3-DOF propagation defaults to central spherical gravity:

```matlab
run.reentry.gravity_model = "CENTRAL_SPHERICAL";
% Optional consistency study with the orbital propagator:
run.reentry.gravity_model = "J2";
```

## SPACEPLANE

Source:

- R. Zhang et al., "Entry trajectory optimization considering blackout zone
  communication constraint," *Advances in Space Research* 77 (2026),
  11407-11417, DOI: 10.1016/j.asr.2025.11.062.

### Aerodynamics and AoA

The paper's conceptual RLV model is used:

```text
CL = -0.041065 + 0.016292 alpha + 0.0002602 alpha^2
CD =  0.080505 - 0.03026 CL + 0.86495 CL^2
```

`alpha` is in degrees. The prescribed velocity profile is:

- 15 deg at and below 2,000 m/s
- linear from 15 to 40 deg between 2,000 and 5,000 m/s
- 40 deg at and above 5,000 m/s

The selected built-in spaceplane shape still supplies its reference area and
nose radius.

### Path constraints

The following paper values are evaluated, but are not actively enforced by a
controller:

- dynamic pressure: 20 kPa
- aerodynamic load: 2.5 g
- heat-rate reference limit: 160 kW/m^2
- bank magnitude: 60 deg
- bank rate: 10 deg/s

### Relay antenna geometry

The relative angle of antenna pointing (RAAP) is evaluated as:

```text
W = r_relay - r_vehicle
beta = acos((A dot W_body) / (|A| |W_body|))
beta <= beam half-angle
```

The velocity-frame-to-body-frame transformation follows the paper's AoA and
bank rotations. The paper uses an upper-body boresight `[0;1;0]`. This project
defaults to the requested aft boresight `[-1;0;0]`:

```matlab
run.reentry.communication.antenna_mount = "AFT";
% Exact paper mounting:
run.reentry.communication.antenna_mount = "PAPER_TOP";
```

The beam half-angle is 45 deg. `X_target`, not a duplicate of the deorbiting
chaser, is propagated as the independently orbiting relay satellite.

Communication geometry is available only when all three tests pass:

- RAAP is inside the antenna half-cone
- slant range is inside the configured minimum/maximum range
- the Earth does not occult the vehicle-to-target line segment

The default maximum range is infinite because no defensible RF link budget is
available yet. A study-specific geometric range limit can be supplied:

```matlab
run.reentry.communication.min_range_m = 0;
run.reentry.communication.max_range_m = 2.0e6;
```

This remains a geometric constraint, not a received-power calculation. It
does not model frequency, antenna gain roll-off, polarization, plasma
attenuation, coding margin, or TDRS service availability.

The code records:

- RAAP for the prescribed bank angle
- the minimum RAAP attainable by any instantaneous bank angle within
  +/-60 deg
- cone, range, Earth-LOS, and combined communication-geometry flags

The second quantity is a geometric feasibility diagnostic, not a bank
controller. Maintaining the constraint throughout entry still requires a
bank-history optimizer or guidance law.

The paper applies the BZC constraint between 80 and 60 km. This project
defaults to continuous target tracking because the requested aft antenna
should keep the target in view. To reproduce the paper's altitude gating:

```matlab
run.reentry.communication.tracking_scope = "BLACKOUT_ONLY";
```

## CAPSULE

Source:

- T. Saito et al., "Guidance strategies for controlled Earth reentry of small
  spacecraft in low Earth orbit," *Acta Astronautica* 229 (2025), 684-697,
  DOI: 10.1016/j.actaastro.2024.12.054.

### Properties

The paper capsule has a mass of 150 kg. This implementation intentionally
uses the requested 60 kg while retaining the paper's published reference
area and drag coefficient:

```text
mass             = 60 kg
reference area   = 0.554 m^2
Cd               = 1.3
nominal L/D      = 0.25 (paper range 0.20-0.30)
beta             = 83.31 kg/m^2
```

For provenance, the configuration also stores the paper's 330 kg total
spacecraft mass and the satellite-module values `Sref = 0.640 m^2` and
`Cd = 2.2`. They do not overwrite this project's chaser geometry or mass.
The paper's hybrid-thruster/RCS deorbit attitude sequence is likewise not
substituted for the existing Phase 3 deorbit design.

The paper does not publish the Mach-10 trim angle or the axial/normal
coefficient table. The reduced model therefore uses constant nominal `L/D`
instead of inventing those data. The logged AoA is a zero-degree trim
surrogate and must not be interpreted as a validated trim solution.

The nose radius is also not stated in the paper. A 0.420 m surrogate (half
the published 0.840 m diameter) is used for Sutton-Graves diagnostics and is
explicitly an assumption.

### Entry and terminal conditions

When capsule mode and `use_paper_entry_conditions` are selected, Phase 3 uses:

- entry interface: 120 km
- target FPA magnitude: 1.16 deg

The paper reference entry speed, 7.97 km/s, is retained as a validation value
but is not forcibly written into the propagated state.

The following paper phase markers are evaluated:

- aeroassist activation when drag acceleration exceeds 0.20 g
- guidance reference cutoff at Mach 3
- analysis termination at 240 m/s, corresponding to parachute deployment

No parachute dynamics are modeled.

### Separation

The default assumes the chaser-capsule stack performs the deorbit maneuver
and the capsule separates at the entry interface with zero separation
impulse:

```matlab
run.reentry.capsule.separation_mode = "ENTRY_INTERFACE";
```

By default, `Chaser_Mass_Init` is treated as the base chaser wet mass and the
60 kg capsule is added to it for Phases 1-3. If an imported
`maneuver.initial_mass_kg` already describes the complete stack, disable the
addition explicitly:

```matlab
run.reentry.capsule.add_to_chaser_initial_mass = false;
```

The carrier's post-separation mass is recorded as
`stack mass - capsule mass`, but its atmospheric trajectory is explicitly
marked `NOT_PROPAGATED`. The model rejects a capsule mass larger than the
entry stack. Position and velocity are continuous across the nominal
separation because the separation impulse is fixed to zero.

The mission budget retains both `Remaining_Mass_kg` (the actively propagated
vehicle) and `Total_Accounted_Mass_kg` (capsule plus the unpropagated carrier
after separation), so separation does not look like an unexplained mass loss.

For a deliberately attached sensitivity case:

```matlab
run.reentry.capsule.separation_mode = "ATTACHED";
```

The attached case retains the capsule aerodynamic area and coefficients, so
it is a ballistic-coefficient sensitivity case rather than a validated
composite-stack aerodynamic model.

The comparison utility also reports a second, explicitly labeled stack
surrogate using the current orbital-drag chaser values (`Cd = 2.2`,
`A = 4.0 m^2`). This is useful for bracketing but is not a substitute for
composite-stack aerodynamics. Study-specific values can be passed directly to
`Reentry_Propagator` as `reference_area_m2`, `cd`, and `nominal_ld`
overrides. The utility's default attached mass is the initial wet stack
(2,060 kg), not the post-burn mass at a particular entry interface; pass an
EI mass as its first argument when that value is available.

Run the direct comparison with:

```matlab
addpath validation
comparison = Compare_Capsule_Separation_Sensitivity();
```

### Uncertainty cases

The explicit-guidance cases use +/-10 percent atmospheric-density and drag
coefficient biases. The broader predictor-corrector campaign uses density
scales 0.9-1.1 and drag-coefficient scales 0.8-1.2. The published capsule
family spans `L/D = 0.20-0.30`, or +/-20 percent about the nominal 0.25.
These can be applied without changing the nominal database:

```matlab
run.reentry.uncertainty.density_scale = 0.9; % 0.9, 1.0, 1.1
run.reentry.uncertainty.cd_scale = 0.8;      % 0.8 ... 1.2
run.reentry.uncertainty.ld_scale = 1.2;      % maps 0.25 to 0.30
```

The paper's 1.0 MW/m^2 heat-rate limit and 200 MJ/m^2 shallow-entry heat-load
reference are stored for diagnostics. The current heating calculation remains
a Sutton-Graves stagnation-point convective estimate.

## Important limitations

- No GPOPS-II or pseudospectral trajectory optimization is included.
- No predictor-corrector, explicit range guidance, terminal bank guidance, or
  RCS control is included.
- The spaceplane antenna constraint is evaluated but not enforced.
- Antenna range is geometric only; no RF or plasma link budget is included.
- The capsule axial/normal coefficient table is not available in the paper.
- The aft antenna is a requested modification; the exact paper antenna points
  along the upper body axis.
- Capsule separation currently supports entry-interface separation or an
  attached sensitivity case; arbitrary separation altitude and dual-body
  atmospheric propagation are not yet implemented.
- The capsule paper uses a 100 Hz RK4 propagation and neglects wind. This
  project defaults to a 0.5 s step and a co-rotating atmosphere for efficient
  preliminary studies. Use `run.reentry.dt_s = 0.01` and
  `run.environment.atmospheric_drag.co_rotate_atmosphere = false` for the
  closer paper-style numerical/velocity convention, then repeat a timestep
  convergence study.
- Geodetic WGS-84 altitude, rarefied-flow aerodynamics, plasma link budget,
  ablation, and TPS response remain future work.
