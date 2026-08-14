# Zhang blackout-zone communication study

This package separates values that are exactly available in Zhang et al.
from assumptions needed by the project's shared re-entry propagator.

```matlab
audit = paperstudies.zhang.run();
status = paperstudies.zhang.status();
test = paperstudies.zhang.selftest(true);
```

The default `run()` performs only deterministic equation, configuration,
and initial-state checks. It does not propagate a trajectory or solve an
optimal-control problem. `published_reference` exposes the terminal box and
nonthermal constraint limits even when no trajectory comparison is run.

An optional open-loop study is available with:

```matlab
forward = paperstudies.zhang.run(struct('forward', true));
```

That result is always labeled `SURROGATE_OPEN_LOOP_FORWARD`. It uses the
published Table 2 cases, AoA and aerodynamic equations, paper TDRS,
`PAPER_TOP` antenna, and `PAPER_BZC` interval, but it uses project assumptions
for mass, reference area, atmosphere, numerical Earth rotation, heating, and
a constant bank angle. It does not reproduce GPOPS-II or invent an optimized
bank history.

The paper does not publish the vehicle mass, reference area, density model,
heating coefficient, complete nondimensionalization constants, initial bank,
GPOPS-II transcription settings, or raw optimized histories. Table 3 and
Fig. 8 also disagree on the heating limit and its units. See `status()` for
the complete machine-readable blocker list.
