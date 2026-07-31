# Codex Checkpoint — Re-entry Paper Models

Date: 2026-07-31  
Workspace: `C:\Users\user\Desktop\Rendezvous_and_Reentry`

## Pause status

- Work was paused at the user's request.
- No MATLAB process is running (`MATLAB_PROCESSES=NONE`).
- All sub-agents have completed; no delegated task remains active.
- No files were staged, committed, or pushed.
- The two source PDFs were not modified. The Git-tracked workspace copies
  `tmp/pdfs/capsule.pdf` and `tmp/pdfs/spaceplane.pdf` were restored after
  temporary extraction artifacts were removed.

## Source papers

- `C:\Users\user\Desktop\IAAT\자료\Reentry\1-s2.0-S027311772501347X-main.pdf`
  - Zhang et al., spaceplane/RLV blackout-zone communication constraint.
- `C:\Users\user\Desktop\IAAT\자료\Reentry\1-s2.0-S0094576524008038-main.pdf`
  - Saito et al., controlled small-capsule Earth re-entry.

## Current working-tree files

Modified:

- `Main_Mission_Simulator.m`
- `Mission_Config.m`
- `Mission_Run_Config.m`
- `README.md`
- `Reentry_Propagator.m`
- `validation/Validate_Reentry_Propagator.m`

New/untracked:

- `docs/REENTRY_PAPER_MODELS.md`
- `validation/Compare_Capsule_Separation_Sensitivity.m`
- `CODEX_CHECKPOINT.md`

## Completed implementation

### Vehicle-mode and mission integration

- Added selectable `SPACEPLANE` and `CAPSULE` re-entry modes.
- Default remains `SPACEPLANE`.
- Added a requested 60 kg capsule mode.
- In capsule mode the default initial stack is:
  `base chaser wet mass + 60 kg capsule`.
- Added `add_to_chaser_initial_mass=false` for configs whose imported initial
  mass already includes the capsule.
- Default separation is at the 120 km Entry Interface with zero separation
  impulse and continuous position/velocity.
- Added entry-stack/capsule/residual-carrier mass validation and ledger.
- Residual carrier trajectory is explicitly marked `NOT_PROPAGATED`.
- Mission budget now distinguishes actively propagated mass from total
  accounted mass so the carrier does not appear to vanish after separation.
- Replaced the former duplicated chaser relay with the independently orbiting
  Target satellite state.

### SPACEPLANE paper model

- Implemented the paper RLV polynomial aerodynamics:
  - `CL = -0.041065 + 0.016292*alpha + 0.0002602*alpha^2`
  - `CD = 0.080505 - 0.03026*CL + 0.86495*CL^2`
- Implemented the speed-scheduled AoA profile:
  - 15 deg through 2,000 m/s
  - linear 15→40 deg from 2,000 to 5,000 m/s
  - 40 deg above 5,000 m/s
- AoA is reevaluated at every RK4 substage.
- Explicit scalar AoA override forces a constant AoA.
- Corrected schedule endpoint handling for non-monotonic user-defined AoA
  values by clamping speed/Mach to the independent-variable endpoints.
- Added paper path-constraint diagnostics:
  - dynamic pressure 20 kPa
  - aerodynamic load 2.5 g
  - heat flux 160 kW/m²
  - bank magnitude 60 deg
  - bank rate 10 deg/s
- Implemented RAAP transformation from the half-velocity frame to body frame
  using AoA and bank.
- Added two antenna mounts:
  - requested aft antenna `[-1;0;0]` (default)
  - exact paper upper-body antenna `[0;1;0]`
- Communication geometry now requires:
  - RAAP inside the beam cone
  - configured range interval
  - no Earth occultation
- Added prescribed-bank RAAP and instantaneous best-bank feasibility scan.
- Added continuous or 60–80 km blackout-only constraint scope.
- An unvisited blackout interval now reports `NOT_EVALUATED`/`NaN`, not a
  vacuous success.
- No attitude/bank control is silently applied; these remain diagnostics.

### CAPSULE paper model

- Implemented requested mass 60 kg with paper:
  - reference area `0.554 m²`
  - `Cd = 1.3`
  - nominal `L/D = 0.25` with published family range 0.20–0.30
  - resulting nominal ballistic coefficient `83.310 kg/m²`
- Stored paper reference provenance:
  - paper capsule mass 150 kg
  - paper total spacecraft mass 330 kg
  - satellite reference area 0.640 m²
  - satellite `Cd = 2.2`
- Added paper Entry Interface references:
  - altitude 120 km
  - speed 7.97 km/s (diagnostic reference)
  - FPA −1.16 deg
- Added paper phase/termination markers:
  - aeroassist activation at drag acceleration 0.20 g
  - guidance window through Mach 3
  - parachute criterion at 240 m/s
- Added density, Cd, and L/D scale factors for uncertainty studies.
- Added paper 1.0 MW/m² heat-rate limit and 200 MJ/m² shallow-entry
  heat-load reference diagnostics.
- Exact paper CA/CN tables and Mach-10 trim angle were not published;
  therefore the code uses an explicitly labeled reduced model instead of
  fabricated tables.

### Shared numerical/diagnostic changes

- Added exact terminal-altitude and terminal-speed event refinement.
- Added explicit termination reasons.
- Added atmosphere-relative FPA and ballistic-coefficient reporting.
- Added relay history names while retaining deprecated chaser aliases.
- Added selectable re-entry gravity model:
  - default `CENTRAL_SPHERICAL` for paper-style 3-DOF propagation
  - optional `J2`
- Fixed the default R-BAR final EI injection mass ledger:
  `charge_final_reentry_fuel=true`.
- Added transparent `reference_area_m2`, `cd`, and `nominal_ld` study
  overrides.

## Tests and commands already executed

### Final static and focused regression test

Executed after the latest AoA-endpoint, blackout-status, R-BAR fuel-ledger,
mass-budget, and central-gravity changes:

```powershell
matlab -batch "files={'Mission_Config.m','Mission_Run_Config.m','Main_Mission_Simulator.m','Reentry_Propagator.m','validation/Validate_Reentry_Propagator.m','validation/Compare_Capsule_Separation_Sensitivity.m'}; n=0; for k=1:numel(files), issues=checkcode(files{k},'-id'); fprintf('%s: %d issues\n',files{k},numel(issues)); n=n+numel(issues); end; assert(n==0,'checkcode reported issues'); addpath('validation'); results=Validate_Reentry_Propagator(); assert(results.passed);"
```

Result:

- MATLAB Code Analyzer: 0 issues in all six checked `.m` files.
- `Validate_Reentry_Propagator`: PASS.
- Terminal-altitude refinement error: `-2.39229e-05 m`.
- Coarse/fine event-time difference: `0 s`.
- Spaceplane test ballistic coefficient: `525.549 kg/m²`.
- 60 kg capsule ballistic coefficient: `83.310 kg/m²`.
- Synthetic initial antenna RAAP: `123.613 deg`.
- The regression suite includes:
  - terminal event accuracy and timestep comparison
  - explicit timeout classification
  - co-rotating-atmosphere option behavior
  - paper AoA schedule and scalar override
  - non-monotonic AoA endpoint clamp
  - antenna range failure
  - Earth occultation
  - blackout interval not evaluated
  - relay/chaser history compatibility
  - capsule mass, beta, L/D, EI reference
  - rejection when capsule mass exceeds entry-stack mass

### Full SPACEPLANE pipeline

Executed successfully:

```powershell
matlab -batch "set(0,'DefaultFigureVisible','off'); Main_Mission_Simulator;"
```

Result at that checkpoint:

- Exit code 0.
- Phase 4 terminated at 20 km with `TERMINAL_ALTITUDE`.
- Open-loop default violated paper constraints:
  - maximum dynamic pressure about 31.685 kPa
  - maximum aero load about 2.541 g
  - aft-antenna RAAP about 112.94–178.92 deg, outside a 45 deg half-cone
- This is expected evidence that the current bank=0 open-loop trajectory is
  not a constrained optimum.

This full run occurred before the final central-gravity/default-budget and
P1 reporting fixes. Focused regression passed afterward, but the exact final
revision has not yet received another full SPACEPLANE run.

### Full CAPSULE pipeline

Executed successfully using a temporary shadow run config; the temporary file
was deleted afterward.

Result at that checkpoint:

- Initial stack: 2,060 kg.
- EI stack after burns: 1,821.636 kg.
- Separated capsule: 60 kg.
- Residual carrier: 1,761.636 kg, marked `NOT_PROPAGATED`.
- Phase 4 termination: `PARACHUTE_SPEED`.
- Duration: about 972.17 s under the then-current J2 setting.
- Final altitude: about 26.06 km.
- Peak heat flux: about 52.231 W/cm².
- Total heat load: about 241.565 MJ/m².
- Peak dynamic pressure: about 2.083 kPa.
- Peak aerodynamic load: about 2.628 g.

This full run also occurred before the final central-gravity/default-budget
and P1 reporting fixes. The exact final revision still needs a full CAPSULE
pipeline rerun.

### Final capsule/stack ballistic sensitivity

Executed after switching the default re-entry gravity model to central
spherical gravity:

```powershell
matlab -batch "addpath('validation'); comparison=Compare_Capsule_Separation_Sensitivity([],1200); assert(comparison.separated_mass_kg==60);"
```

Result:

- 60 kg capsule: `beta = 83.31 kg/m²`.
- 2,060 kg stack with the same capsule CdA:
  `beta = 2860.32 kg/m²`.
- At the common time 1,046.5 s:
  - position difference about 1,691.935 km
  - velocity-vector difference about 7,271.999 m/s
  - altitude about 26.078 km versus 144.573 km
- Current-chaser CdA surrogate (`Cd=2.2`, `A=4.0 m²`):
  `beta = 234.09 kg/m²`.
- At the common time 792.5 s:
  - capsule versus surrogate position difference about 380.070 km
  - velocity-vector difference about 4,709.684 m/s
  - altitude about 65.868 km versus 20.000 km

These are ballistic-coefficient sensitivity cases, not validated
composite-stack aerodynamics.

## Interpretation preserved for the next session

- If the capsule remains attached until the Entry Interface and zero
  separation impulse is assumed, capsule and stack have exactly the same
  position/velocity at separation by construction.
- Before EI, mass alone cannot determine drag divergence:
  `a_D = q*Cd*A/m = q/beta`; composite `CdA` is required.
- After separation, the assumed 60 kg capsule beta is much smaller than either
  tested stack surrogate, so state divergence can become large.
- The aft antenna does not “track” the Target by itself. The model evaluates
  whether the prescribed AoA/bank trajectory satisfies geometry and whether
  some instantaneous bank within ±60 deg could satisfy RAAP. A bank-history
  optimizer or guidance law is still needed to maintain the constraint.

## Remaining work

Required next checks:

1. Rerun the full default SPACEPLANE pipeline on the exact current revision.
2. Rerun the full CAPSULE pipeline on the exact current revision and confirm:
   - base chaser + capsule stack addition
   - burn mass depletion
   - EI split conservation
   - `Remaining_Mass_kg` versus `Total_Accounted_Mass_kg`
   - `PARACHUTE_SPEED` termination
3. Run `git diff --check` and review the final diff after those integrations.
4. Decide whether to stage/commit; nothing is currently staged.

Useful additional regression work:

- RAAP exactly aligned and exactly on the 45 deg cone boundary.
- Automated Main-level stack-add → burn-depletion → split conservation test.
- Timestep convergence for the full paper capsule case, especially
  `dt = 0.5, 0.1, 0.01 s`.
- Density/Cd/L/D uncertainty matrix execution and result table.

Deliberately not implemented:

- GPOPS-II/pseudospectral SPACEPLANE trajectory optimization.
- Closed-loop bank/attitude control.
- Capsule explicit or predictor-corrector guidance.
- RF link budget, antenna gain pattern, frequency/polarization/coding margin,
  or plasma attenuation model.
- Arbitrary-altitude capsule separation and dual capsule/carrier propagation.
- Validated composite chaser-capsule aerodynamics.
- Missing capsule CA/CN/Mach trim database.
- Ablation, TPS response, rarefied-flow model, and WGS-84 geodetic dynamics.

## Exact next starting point

Start by reading this file, then:

```powershell
git status --short
matlab -batch "set(0,'DefaultFigureVisible','off'); Main_Mission_Simulator;"
```

After the default SPACEPLANE run, create a temporary/shadow
`Mission_Run_Config.m` selecting `CAPSULE` without editing the saved default,
run the full simulator, delete the temporary override, and inspect the final
mass-budget columns. Do not assume the earlier full-pipeline numbers are the
final central-gravity values.
