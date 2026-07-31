# Codex Checkpoint — Re-entry Paper Audit and Reduced Models

Date: 2026-07-31  
Workspace: `C:\Users\user\Desktop\Rendezvous_and_Reentry`

## Current status

- The interrupted PDF reviews were restarted and completed.
- Spaceplane paper: all 11/11 pages rendered and visually inspected; key
  tables/figures were rechecked at high resolution.
- Capsule paper: all 14/14 pages rendered and visually inspected; equations
  and tables were cross-checked with text extraction.
- Both independent paper-review agents completed. The unused architecture
  agent was stopped before starting.
- No MATLAB process is running.
- Default run mode has been restored to `SPACEPLANE`.
- Temporary page renders, extracted text, and `.mat` sensitivity artifacts
  were deleted. The Git-tracked source copies remain:
  - `tmp/pdfs/spaceplane.pdf`
  - `tmp/pdfs/capsule.pdf`
- No file is staged or committed.
- Current `HEAD` remains `373d237` (`Checkpoint`).

## Working tree

Modified:

- `CODEX_CHECKPOINT.md`
- `Main_Mission_Simulator.m`
- `Mission_Config.m`
- `Mission_Run_Config.m`
- `README.md`
- `Reentry_Propagator.m`
- `docs/ASSUMPTIONS_AND_LIMITATIONS.md`
- `docs/REENTRY_PAPER_MODELS.md`
- `validation/Validate_Reentry_Propagator.m`

Untracked: none.

## Source papers and audit conclusions

### Zhang et al. — SPACEPLANE/RLV

Source: `tmp/pdfs/spaceplane.pdf`

The current mode is now explicitly described as:

> Paper-derived RLV aerodynamics and RAAP geometry in a forward propagator,
> with either a user-requested aft antenna/dynamic mission target or an
> optional paper-reference TDRS; no paper trajectory optimization or BZC
> enforcement.

Confirmed paper values:

- `CL = -0.041065 + 0.016292*alpha + 0.0002602*alpha^2`
- `CD = 0.080505 - 0.03026*CL + 0.86495*CL^2`
- AoA: 15 deg through 2 km/s, linear to 40 deg at 5 km/s, then 40 deg
- paper antenna: `+Y_B = [0;1;0]`
- beam half-angle: 45 deg
- TDRS: Earth-fixed lon/lat 77/0 deg, radius 42,164 km
- `q <= 20 kPa`, `n <= 2.5`, `|bank| <= 60 deg`,
  `|bank_rate| <= 10 deg/s`
- initial altitude 80 km; terminal altitude/speed 25 km/800 m/s
- BZC is active from the initial state until the first downward 60 km
  crossing and remains active during short skips above 80 km

Important paper inconsistency:

- Table 3 gives a heat reference of `160,000 W`.
- Fig. 8(c) appears to place the boundary near `1.6 MW/m^2`.
- Eq. (32) uses `k_q*sqrt(rho)*V^3.15`, but `k_q` is not provided.

The current Sutton-Graves result is therefore not used to claim paper heat
compliance.

### Saito et al. — CAPSULE

Source: `tmp/pdfs/capsule.pdf`

The current mode is now explicitly described as a 60 kg reduced surrogate:

- mass: 60 kg, intentionally replacing the paper's 150 kg
- reference area: 0.554 m²
- `Cd = 1.3`
- nominal constant `L/D = 0.25`
- beta: 83.310 kg/m²
- paper 150 kg beta with the same `Cd*A`: 208.28 kg/m²

Confirmed paper details now stored or documented:

- four Table 4 entry cases, including geodetic altitude 121.01–121.15 km,
  speed 7.97 km/s, FPA -1.16/-1.17 deg, longitude, latitude, azimuth,
  and 4,852–4,920 km range-to-go
- target latitude/longitude 22.74/158.78 deg
- first 0.20 g drag-acceleration crossing, Mach 3 guidance cutoff, and
  240 m/s termination
- density uncertainty ±10 percent
- explicit-case `Cd` uncertainty ±10 percent
- RPC-case `Cd` uncertainty ±20 percent
- recession-style `L/D` change about 10 percent
- capsule-family envelope `L/D = 0.20–0.30`, kept distinct from uncertainty

The paper uses proprietary `CA,CN(Mach,alpha)` data and an unpublished
Mach-10 trim AoA. Exact capsule aerodynamics cannot be reconstructed from
the PDF.

The paper heating law is Detra-Kemp-Riddell Eq. (28), not Sutton-Graves.
The 1 MW/m² value is an HSRC allowable heat rate; about 200 MJ/m² is an
observed shallow-entry load, not a declared allowable total-load limit.

## Implemented corrections in this continuation

### SPACEPLANE

- Removed silent clipping of negative `CL` from the exact paper polynomial.
- Added two explicit relay modes:
  - `MISSION_TARGET_DYNAMIC_ORBIT` — default requested generalization
  - `PAPER_TDRS_STATIC_EARTH_FIXED` — analytic paper-reference GEO relay
- Added paper TDRS validation and a regression at 42,164 km radius.
- Implemented stateful `PAPER_BZC`/`BLACKOUT_ONLY` gating for mission runs
  that begin above the paper's 80 km initial point: wait for the first
  downward 80 km entry, remain active through later skips above 80 km, then
  exit permanently at the first downward 60 km crossing.
- Retained `ALTITUDE_BAND` only as a non-paper literal 60–80 km option.
- Added published initial scenarios and terminal-state references.
- Renamed the interpretation of the bank scan to
  `instantaneous_bank_geometric_reachable`; it is not guidance feasibility.
- A disabled bank scan now remains explicitly unevaluated (`NaN`) instead
  of treating the prescribed bank angle as the best attainable result.
- Included active antenna-geometry failures in overall constraint status and
  the main-run warning.
- Stored and documented the paper's internally inconsistent heat references.

### CAPSULE

- Added the full Table 4 entry-state provenance and paper target coordinates.
- Added the paper capsule inertia and 150 kg ballistic-coefficient reference.
- Labeled `use_paper_entry_conditions` as `ALTITUDE_AND_FPA_ONLY`; it does not
  force speed, position, azimuth, or range-to-go.
- Made the 0.20 g activation marker latch after first crossing.
- Disabled the generic 20 km stop by default in capsule mode so the paper's
  240 m/s criterion has priority.
- Added a zero-altitude `GROUND_IMPACT` safety event to prevent propagation
  below the spherical surface if 240 m/s is never reached.
- Added `run.reentry.capsule.altitude_termination_enabled` as an explicit
  optional safety-stop control.
- Separated density, explicit `Cd`, RPC `Cd`, recession `L/D`, and family
  `L/D` references.
- Documented that entry-interface separation is a project assumption; the
  paper separates after its deorbit/standby/de-spin/release-attitude sequence
  and before entry interface.

### Shared rigor/status changes

- Heating model and paper heating model are now separately reported.
- A non-comparable heat model yields:
  - `heat_constraint_evaluated = false`
  - `heat_constraint_satisfied = NaN`
  - aggregate paper-path status `NaN` unless another evaluated constraint
    is definitely violated
- Raw Sutton-Graves values remain available for relative sensitivity work.
- Added tri-state overall constraint reporting:
  - definite violation when an evaluated path or antenna constraint fails
  - success only when every configured constraint is evaluated and passes
  - `NaN` when compliance cannot be established
- Applied the same Kleene-style ordering to partially evaluated antenna/BZC
  intervals: any known failure is `false`; otherwise a fully evaluated pass
  is `true`; only genuinely unresolved cases remain `NaN`.
- Disabled altitude termination now reports
  `terminal_altitude_error_m = NaN`; capsule plots and environment summaries
  advertise the 240 m/s nominal stop and 0 km ground-safety event instead of
  presenting the inactive 20 km stop as operational.
- Corrected the documentation error that recommended
  `co_rotate_atmosphere=false`. Both papers' rotating-Earth entry equations
  correspond to the ECI implementation with `co_rotate_atmosphere=true`.
- Updated README and assumptions/limitations to distinguish paper reference,
  reduced implementation, and surrogate diagnostic.

## Verification completed

### Static analysis and focused regression

Command:

```powershell
matlab -batch "files={'Mission_Config.m','Mission_Run_Config.m','Main_Mission_Simulator.m','Reentry_Propagator.m','validation/Validate_Reentry_Propagator.m','validation/Compare_Capsule_Separation_Sensitivity.m'}; n=0; for k=1:numel(files), issues=checkcode(files{k},'-id'); fprintf('%s: %d issues\n',files{k},numel(issues)); n=n+numel(issues); end; assert(n==0,'checkcode reported issues'); addpath('validation'); results=Validate_Reentry_Propagator(); assert(results.passed);"
```

Result:

- Code Analyzer: 0 issues in all six files.
- `Validate_Reentry_Propagator`: PASS.
- terminal-altitude refinement error: `-2.39229e-05 m`
- coarse/fine event-time difference: `0 s`
- synthetic spaceplane beta: `525.549 kg/m²`
- capsule beta: `83.310 kg/m²`
- initial synthetic aft-antenna RAAP: `123.613 deg`

New regression coverage includes:

- negative zero-AoA `CL` retained
- BZC inactive before the first downward 80 km entry
- BZC activation at 80 km and persistence during a post-entry skip above it
- permanent BZC exit after the first downward 60 km crossing
- deterministic initialization inside the BZC and below its 60 km exit
- BZC unevaluated below 60 km
- known antenna failure retained under partial geometry evaluation
- disabled instantaneous bank scan retained as unevaluated
- static paper TDRS radius
- non-comparable paper heat status
- latched 0.20 g marker
- capsule altitude termination disabled by default
- recession-style `L/D` uncertainty provenance

### Full SPACEPLANE pipeline

Result: PASS, exit code 0.

- termination: `TERMINAL_ALTITUDE`, 20 km
- entry duration: about 1819.85 s
- max dynamic pressure: about 31.406 kPa
- max aero load: about 2.568 g
- aft RAAP: about 112.94–178.91 deg
- target LOS eventually occulted
- overall evaluated constraints correctly reported as violated
- heat output explicitly labeled `SUTTON_GRAVES_SURROGATE`
- active and total-accounted final mass: 1638.67 kg

This is expected evidence that the default bank=0 trajectory is not a
constraint-satisfying optimum.

### Full CAPSULE pipeline

Result: PASS, exit code 0.

- initial stack: 2060 kg
- entry-interface stack: 1821.636 kg
- separated capsule: 60 kg
- residual carrier: 1761.636 kg, `NOT_PROPAGATED`
- termination: `PARACHUTE_SPEED`
- duration: about 977.79 s
- event altitude: about 26.06 km
- max Sutton-Graves surrogate heat: about 52.298 W/cm²
- total surrogate heat: about 242.191 MJ/m²
- max dynamic pressure: about 2.081 kPa
- max aero load: about 2.625 g
- active final mass: 60 kg
- total-accounted mass: 1821.636 kg

Assertions passed for 240 m/s priority, non-comparable heat status, capsule
mass, residual-carrier status, and total mass conservation.

### Non-model command errors encountered

- One SPACEPLANE invocation failed before execution because nested
  PowerShell/MATLAB quote escaping broke a string comparison. The corrected
  command passed.
- The first CAPSULE integration completed physically but the final assertion
  referenced a nonexistent local field name. The assertion was corrected to
  `reentry_atmo_info.carrier_post_separation_status`; the repeated full run
  passed.

## Deliberately not implemented

- GPOPS-II/two-phase pseudospectral SPACEPLANE optimization
- bank-rate control and `integral(abs(bank_rate))` objective
- terminal longitude/latitude/range targeting and published-figure
  cross-validation
- capsule HTS/RCS six-phase deorbit and attitude dynamics
- capsule explicit guidance Eqs. (13)–(22) / Table 5
- capsule RPC guidance Eqs. (23)–(25), lift/drag estimator, terminal bank
  guidance, and bank reversals
- proprietary capsule `CA,CN(Mach,alpha)` database and trim AoA
- paper-compatible thermochemical heating calculations
- RF link budget, gain pattern, polarization, coding margin, service
  availability, and plasma attenuation
- arbitrary-altitude separation with dual capsule/carrier propagation
- validated composite chaser-capsule aerodynamics
- geodetic WGS-84 entry dynamics, rarefied flow, ablation, and TPS response

## Recommended next work

1. Add a separate paper-regression runner for Zhang Tables 2–3 using the
   80 km initial states, 25 km/800 m/s terminal references, paper TDRS, and
   `PAPER_TOP`. This would validate the forward equations without claiming
   optimization.
2. Run capsule timestep convergence at 0.5, 0.1, and 0.01 s.
3. Add automatic explicit (15-case) and RPC-reference (225-case) uncertainty
   campaign generation, clearly labeled as reduced-model sensitivity rather
   than guidance reproduction.
4. Obtain or bound composite stack `Cd*A` before selecting a physical
   pre-entry separation epoch.
5. Obtain real antenna/link data before setting a finite communication range.

## Exact next starting point

```powershell
git status --short
git diff --check
git diff -- Main_Mission_Simulator.m Mission_Config.m Mission_Run_Config.m Reentry_Propagator.m
git diff -- README.md docs/ASSUMPTIONS_AND_LIMITATIONS.md docs/REENTRY_PAPER_MODELS.md
git diff -- validation/Validate_Reentry_Propagator.m CODEX_CHECKPOINT.md
```

Then decide whether to commit the nine modified files. The default run
configuration is already restored to `SPACEPLANE`.
