# Paper-Based Re-entry Modes

The atmospheric-entry propagator supports two selectable reduced-order
vehicle modes:

```matlab
run.reentry.vehicle_mode = "SPACEPLANE";
% or
run.reentry.vehicle_mode = "CAPSULE";
```

These are deterministic reduced-order trajectory and diagnostic models. They
are not paper-reproduction models and do not include the papers' optimal
control or closed-loop bank guidance.

The translational state is propagated in Cartesian ECI coordinates with
central spherical gravity and atmosphere-relative velocity
`v_rel = v_ECI - omega x r`. Under spherical-Earth/no-wind assumptions this
is physically equivalent to the rotating spherical-coordinate 3-DOF form
used in both papers, although it has not been cross-validated against the
published trajectories:

```matlab
run.reentry.gravity_model = "CENTRAL_SPHERICAL";
% Optional consistency study with the orbital propagator:
run.reentry.gravity_model = "J2";
```

The configuration and summary deliberately distinguish three levels:

- `paper reference`: a value or assumption transcribed from a source paper
- `implemented reduced model`: a compatible approximation used in propagation
- `surrogate diagnostic`: a useful number that is not quantitatively
  comparable with the paper result

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
nose radius. The paper does not publish either value in this article, so
combining its coefficient polynomials with a built-in shape assumes a
compatible coefficient/reference-area convention. This is an unresolved
model-form uncertainty, not a paper-derived vehicle definition.

### Path constraints

The following values come from Table 3:

- dynamic pressure: 20 kPa
- aerodynamic load: 2.5 g
- bank magnitude: 60 deg
- bank rate: 10 deg/s

Those four constraints are evaluated but not actively enforced. Heating is
handled differently. Zhang Eq. (32) uses
`Qdot = k_q*sqrt(rho)*V^3.15`, while the current code computes a
Sutton-Graves `sqrt(rho/Rn)*V^3` surrogate. The paper does not provide the
required `k_q`, and its own presentation is inconsistent: Table 3 gives
`160,000 W`, whereas Fig. 8(c) appears to show a boundary near
`1.6 MW/m^2`. Consequently the configured 160 kW value is retained only as
provenance; no paper heat pass/fail is claimed. The summary reports
`heat_constraint_evaluated = false` and leaves the aggregate status
unevaluated unless another implemented constraint is definitely violated.

The paper cases start at 80 km and terminate at 25 km/800 m/s, targeting
longitude 111.3 deg and latitude 42.2 deg. The three initial
longitude/latitude/speed/FPA/heading rows are stored in
`spaceplane.paper_initial_scenarios`. Normal mission runs instead begin at
the mission-generated 120 km interface and currently stop at 20 km; therefore
they do not reproduce Figs. 4-11.

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

The beam half-angle is 45 deg. Two explicitly different relay cases are
available:

```matlab
% User-requested generalization; this is the default.
run.reentry.communication.relay_mode = "MISSION_TARGET_DYNAMIC_ORBIT";

% Paper reference: equatorial GEO, Earth-fixed lon/lat = 77/0 deg,
% geocentric radius 42,164 km.
run.reentry.communication.relay_mode = "PAPER_TDRS_STATIC_EARTH_FIXED";
```

In the first case, `X_target`, not a duplicate of the deorbiting chaser, is
propagated as the independently orbiting relay. In the second case an
analytic Earth-fixed TDRS state is generated; the supplied target state is
ignored. The paper treats this relay as stationary relative to Earth over
the 800-2,000 s RLV flight.

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

The second quantity is an instantaneous, independent bank-only scan. It does
not enforce the 10 deg/s bank-rate limit, propagate the trajectory resulting
from each bank choice, or satisfy terminal conditions. It is therefore
reported as `instantaneous_bank_geometric_reachable`, not as guidance
feasibility. If `evaluate_bank_feasibility` is disabled, the scan histories
and summary status remain `NaN` rather than reinterpreting the prescribed
bank angle as an optimization result.

The paper starts at the nominal blackout upper boundary of 80 km, and
Eq. (34) applies BZC from that initial point until the first downward 60 km
crossing. It explicitly remains active during a short skip above 80 km.
Because the integrated mission trajectory begins at 120 km, `PAPER_BZC`
uses a state machine: phase 0 waits for the first downward entry through
80 km, phase 1 remains active through any subsequent skip, and phase 2
begins at the first downward 60 km exit and never reactivates. An initial
state already within 60–80 km starts in phase 1 because its prehistory is
unknown; an initial state below 60 km starts in terminal phase 2. This project
defaults to continuous target tracking for the requested aft antenna. The
paper gating is:

```matlab
run.reentry.communication.tracking_scope = "PAPER_BZC";
```

`"BLACKOUT_ONLY"` is retained as an alias. `"ALTITUDE_BAND"` gives the older
literal 60-80 km band, but that is not the paper rule. Antenna failure is
included in `overall_constraints_satisfied` and in the main-run warning.

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

The corresponding 150 kg paper value is `208.28 kg/m^2`; with the same
`Cd*A`, the requested 60 kg vehicle has 2.5 times the aerodynamic
acceleration at a common state. Published guidance performance for the
150 kg capsule therefore cannot validate the 60 kg trajectory.

For provenance, the configuration also stores the paper's 330 kg total
spacecraft mass and the satellite-module values `Sref = 0.640 m^2` and
`Cd = 2.2`. They do not overwrite this project's chaser geometry or mass.
The article introduction says 300 kg, while Sec. 2.1/Table 1 gives
180 + 150 = 330 kg; the configuration follows the detailed table.
The paper's hybrid-thruster/RCS deorbit attitude sequence is likewise not
substituted for the existing Phase 3 deorbit design. The source case uses a
450 km sun-synchronous orbit, a 172 N HTS, spin stabilization, a 210 s
firing, 540 s residual-thrust standby, de-spin, release-attitude setup, and
then capsule separation. The current Hohmann/generic finite-burn Phase 3 is a
separate mission model.

The paper does not publish the Mach-10 trim angle or the axial/normal
coefficient table. The reduced model therefore uses constant nominal `L/D`
instead of inventing those data. The logged AoA is a zero-degree trim
surrogate and must not be interpreted as a validated trim solution.

The nose radius is also not stated in the paper. A 0.420 m surrogate (half
the published 0.840 m diameter) is used for Sutton-Graves diagnostics and is
explicitly an assumption.

### Entry and terminal conditions

When capsule mode and `use_paper_entry_conditions` are selected, Phase 3
currently applies only:

- entry interface: 120 km
- target FPA magnitude: 1.16 deg

The scope is recorded as `ALTITUDE_AND_FPA_ONLY`. The paper's four Table 4
states are geodetic altitudes 121.01-121.15 km, speed 7.97 km/s, FPA
-1.16/-1.17 deg, and range-to-go 4,852-4,920 km. Their longitude, latitude,
azimuth, range-to-go, and target at 22.74 deg N/158.78 deg E are stored as
provenance but are not forced into the propagated mission state. The 120 km
spherical event must not be treated as numerically identical to those
geodetic states.

The following paper phase markers are evaluated:

- latched aeroassist activation at the first drag-acceleration crossing of
  0.20 g
- guidance reference cutoff at Mach 3
- analysis termination at 240 m/s, corresponding to parachute deployment

Capsule mode disables the generic 20 km termination so the paper's 240 m/s
criterion has priority. A run that does not reach it ends at `max_time`; no
parachute dynamics are modeled. A zero-altitude `GROUND_IMPACT` safety event
prevents propagation below the spherical surface. A study-specific altitude
stop can be re-enabled explicitly with:

```matlab
run.reentry.capsule.altitude_termination_enabled = true;
```

### Separation

The default assumes the chaser-capsule stack performs the deorbit maneuver
and the capsule separates at the entry interface with zero separation
impulse:

```matlab
run.reentry.capsule.separation_mode = "ENTRY_INTERFACE";
```

This is a project assumption, not the paper's separation event. The paper
separates after deorbit firing, standby, de-spin, and release-attitude setup,
before the capsule reaches the entry interface. A defensible arbitrary
separation epoch additionally requires composite stack and separated-carrier
`Cd*A`, which are not available.

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

The paper contains three related but distinct quantities:

- explicit Table 7 cases: density +/-10 percent and `Cd` +/-10 percent,
  varied one at a time for each range
- RPC Table 9 campaign: density 0.9-1.1 and `Cd` 0.8-1.2 in a full matrix
- capsule-family envelope: `L/D = 0.20-0.30`; separately, the
  recession comparison suggests about 10 percent change from a nominal shape

The family envelope is not itself a +/-20 percent uncertainty distribution.
The paper's prose sometimes says lift-to-drag uncertainty where its tables
are labeled drag-coefficient scale, so both independent study knobs are kept:

```matlab
run.reentry.uncertainty.density_scale = 0.9; % 0.9, 1.0, 1.1
run.reentry.uncertainty.cd_scale = 0.8;      % 0.8 ... 1.2
run.reentry.uncertainty.ld_scale = 1.1;      % recession-style sensitivity
% Family-envelope endpoint, not a stated probability bound:
run.reentry.uncertainty.ld_scale = 1.2;      % maps 0.25 to 0.30
```

The paper uses a Detra-Kemp-Riddell Eq. (28) with a `V^3.15` term,
sea-level normalization, and stagnation/hot-wall enthalpy terms. The current
calculation is a Sutton-Graves surrogate and cannot be validated against the
paper without the missing nose radius and a consistent thermochemical model.
The HSRC 1.0 MW/m^2 allowable heat rate and the roughly 200 MJ/m^2 observed
shallow-entry heat load are stored as references. The latter is not a stated
allowable limit. Both direct pass/fail comparisons are therefore marked
unevaluated; raw surrogate values remain available for relative studies.

## Important limitations

- No GPOPS-II or pseudospectral trajectory optimization is included.
- No predictor-corrector, explicit range guidance, terminal bank guidance, or
  RCS control is included.
- The spaceplane antenna constraint is evaluated and reported in the overall
  status but is not enforced. The paper's bank-rate-as-control two-phase OCP
  and `integral(abs(bank_rate))` objective are absent.
- Antenna range is geometric only; no RF or plasma link budget is included.
- The capsule axial/normal coefficient table is not available in the paper.
- The aft antenna is a requested modification; the exact paper antenna points
  along the upper body axis.
- Capsule separation currently supports entry-interface separation or an
  attached sensitivity case; arbitrary separation altitude and dual-body
  atmospheric propagation are not yet implemented.
- The capsule paper uses a 100 Hz RK4 propagation, spherical gravity, Earth
  rotation, and no upper-atmosphere wind. This project defaults to a 0.5 s
  step. Use `run.reentry.dt_s = 0.01` and retain
  `run.environment.atmospheric_drag.co_rotate_atmosphere = true` for the
  closer paper-style numerical/velocity convention, then perform a timestep
  convergence study. Setting it false incorrectly uses inertial velocity as
  aerodynamic velocity.
- The paper's `CA,CN(Mach,alpha)` database, Mach-10 trim AoA, fading-filter
  coefficient, predictor perturbation, cross-range deadband, and RCS control
  gains are not published. Exact capsule guidance replication is impossible
  from this article alone.
- Spaceplane ISA76 density, built-in reference area, and nose radius are
  project assumptions; the Zhang article does not provide them.
- Geodetic WGS-84 altitude, rarefied-flow aerodynamics, plasma link budget,
  ablation, and TPS response remain future work.
