# Paper Reproduction Framework

## Purpose

This framework keeps two questions separate:

1. Can the mission simulator propagate a candidate atmospheric-entry state
   with one shared, tested physical kernel?
2. Can the published Zhang and Saito results be reproduced numerically from
   the information printed in the papers?

The answer to the first question can be tested inside this repository. The
answer to the second is currently only partial because both papers omit inputs
that materially determine their trajectories. The study runners therefore do
not tune hidden assumptions until a published plot happens to match. Every
quantity is instead marked as published, convention-dependent, surrogate, or
unavailable.

## Architecture

```text
Main_Mission_Simulator
        |
        v
Reentry_Propagator  -- mission policy, relay propagation, BZC latch,
        |              history, constraint summary, terminal events
        v
+reentry_core       -- shared 7-state dynamics, RK4, atmosphere/aero,
                       event refinement, LOS, RAAP, antenna geometry

+paperstudies/+zhang -- published Zhang conditions and equation audit
+paperstudies/+saito -- published Saito conditions and equation audit
        |                      |
        +----------+-----------+
                   v
             +reentry_core
          optional surrogate forward runs
```

`Reentry_Propagator` retains its original public signature. Existing mission
runs therefore continue to use the same entry interface. Its numerically pure
physics was extracted into the `+reentry_core` MATLAB package so that paper
studies can exercise the same implementation without copying a second dynamics
model. Paper-specific optimization, guidance, relay policy, and reporting do
not belong in that physical kernel.

## Reproduction labels

The study outputs use the following meanings:

- `EXACT` or `PUBLISHED_EXACT`: directly transcribed public constants,
  tables, equations, switches, or geometry.
- `CONVENTION_REQUIRED`: the paper gives an equation or state but omits a
  frame, unit, interpolation, or scaling convention needed to evaluate it.
- `SURROGATE` or `SURROGATE_ONLY`: a transparent project assumption makes a
  forward calculation possible; it is not attributed to the paper.
- `UNAVAILABLE`: a material input is not printed and no defensible exact
  numerical reconstruction can be made.

These are provenance labels, not accuracy scores. An exactly transcribed
equation can still be insufficient for an exact trajectory when its inputs are
unavailable.

## Running the studies

From the repository root in MATLAB:

```matlab
% Fast selected-regression audit: source anchors, public equation examples,
% coordinate construction, case enumeration, and provenance.
report = Run_Paper_Reproduction_Suite();

% Individual paper audits.
zhang = paperstudies.zhang.run();
saito = paperstudies.saito.run();

% Optional reduced forward propagation through the shared physical kernel.
% These results remain explicitly classified as surrogates.
options.forward = true;
report = Run_Paper_Reproduction_Suite(options);
```

The default suite is intentionally fast, runs selected deterministic
self-tests, and does not run an optimizer or claim to recreate the published
flight histories. It does not independently prove every transcribed table cell;
the full PDFs were manually audited, while regression anchors protect the most
important constants and algebra. The optional forward mode uses the
paper initial conditions where possible, but fills unavailable vehicle or
aerodynamic inputs only with the assumptions listed in the returned
`forward.assumptions` and paper-specific status/manifest structures.
Both adapters force the repository-owned ISA76 implementation so results do
not change merely because MATLAB Aerospace Toolbox is present; this density
implementation remains a project surrogate where the paper does not publish
its numerical atmosphere implementation.

Run the regression checks with:

```matlab
addpath validation
Run_All_Reentry_Validations

% Or run the checks separately:
Validate_Reentry_Core_Equivalence
Validate_Reentry_Propagator
Validate_Paper_Reproduction
```

For reproducible tests, add only the repository root and the `validation`
folder. Avoid `addpath(genpath(project_root))`, because the `legacy` directory
contains old functions with names that may shadow active code.

## Evidence supplied by this repository

### Shared physical implementation

`Validate_Reentry_Core_Equivalence` compares the extracted physical kernel
against an independent frozen transcription of the pre-extraction equations.
It covers central and J2 gravity, rotating and non-rotating atmosphere options,
spaceplane and capsule aerodynamics, RK4 integration, event refinement, LOS,
RAAP, and fixed end-to-end reference states. This establishes that the
refactor did not intentionally change the existing mission model. It is a
software regression result, not external validation of the model assumptions.

### Zhang study

The Zhang package contains:

- all three Table 2 initial cases and the published terminal box;
- the speed-scheduled angle of attack and polynomial `CL/CD` equations;
- the spherical dimensionless equations in a standalone algebraic helper;
- the paper TDRS position, upper-body antenna, 45 degree cone, RAAP
  transformation, and latched 80-to-60 km blackout-zone rule;
- the published non-thermal path limits and the two-phase OCP structure;
- an optional open-loop constant-bank forward adapter.

It does not reconstruct the paper's optimal trajectories. Vehicle mass,
reference area, density model, numerical nondimensionalization constants,
heating coefficient, initial bank, GPOPS-II transcription settings, optimized
bank history, and raw trajectory data are unavailable. The heat constraint is
also internally ambiguous between Table 3 and Fig. 8, so it is not converted
into a paper-compliance Boolean.

### Saito study

The Saito package contains:

- Tables 1, 3, 4, 5, 6, 7, 8, 9, 10, and 11 in machine-readable form;
- all four entry-interface state constructors and coordinate round-trip
  diagnostics;
- published algebraic helpers for Eqs. 10-26, with Eq. 17 distance scaling
  and unpublished explicit-guidance constants marked convention-required;
- the 15-case explicit-guidance and 225-case RPC uncertainty grids;
- published activation, cutoff, termination, integration, control, and
  guidance rates;
- an optional 150 kg, constant-`Cd`, constant-`L/D` capsule propagation.

It cannot reproduce the published landing dispersions, 6-DOF deorbit, or
thermal results exactly. The proprietary HSRC `CA/CN(Mach, AoA)` database,
Mach-10 trim angle, cross-range deadband, fading-filter coefficient, predictor
perturbation and settings, terminal-bank transition details, RW/RCS controller,
Eq. 17 distance redimensionalization, explicit-guidance `V0/K1/K2/K3`, the
individual Table-6 entry states/ranges, and complete heating-law inputs are
missing. The printed Table 4 range-to-go
values are also inconsistent with the printed coordinates under a conventional
great-circle Earth radius; the package reports that discrepancy instead of
hiding it with a fitted radius.

## What may and may not be claimed

The current tests support these claims:

- the mission and paper-study paths use the same atmospheric-entry physical
  implementation;
- the extracted kernel is numerically equivalent to the previous integrated
  implementation over the tested cases;
- the named published tables are machine-readable, while selected source
  anchors, public-equation examples, state conversions, and uncertainty-case
  enumeration are regression tested;
- explicitly named surrogate cases can be used for sensitivity studies and
  preliminary-parameter generation.

They do not yet support these claims:

- full numerical reproduction of either paper's trajectories or footprints;
- independent proof that every transcribed source cell is correct solely from
  the automated selected-regression audit;
- independent validation merely because a surrogate curve looks similar to a
  published figure;
- validation of paper heating limits with the current Sutton-Graves diagnostic;
- flight-quality guidance, navigation, control, or aerodynamic fidelity.

The defensible research workflow is therefore: reproduce what is public,
report the missing inputs, perform convergence and sensitivity studies with
the shared kernel, and reserve trajectory-level validation claims until raw
reference histories or the omitted aerodynamic/guidance data become available.

## Extending toward full reproduction

For Zhang, add the exact vehicle and density/nondimensionalization data and an
independent GPOPS-II-compatible OCP transcription, then compare time histories
against numeric reference data rather than digitized plots alone. For Saito,
add the HSRC aerodynamic database and all omitted guidance/controller
parameters before comparing landing dispersions. New data should be stored
beside its source metadata and assigned a provenance label; it should not
silently replace a surrogate in `Mission_Config`.
