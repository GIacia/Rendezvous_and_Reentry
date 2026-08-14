# Saito paper-study package

This package separates the public Saito et al. (2025) conditions from the
normal mission configuration. Its default run audits the printed data and
equations without propagating a trajectory:

```matlab
audit = paperstudies.saito.run();
manifest = paperstudies.saito.status();
test = paperstudies.saito.selftest();
```

`config.m` stores Tables 1, 3, 4, 5, 6, and 7-11 plus the published re-entry
switches and execution rates. `entry_state.m` constructs the four Table-4
states; `roundtrip_entry_state.m` checks their coordinate reconstruction. The
Table-6 helper generates the 27 exact ignition epochs. It deliberately leaves
the individual entry states and range-to-go values as `NaN`: the paper only
publishes the 2369-7031 km aggregate bounds and approximate spacing. The two
uncertainty helpers generate all 15 explicit-guidance and 225 RPC cases. The
`eq*.m` functions expose the public dynamics, force conversion, explicit
guidance, predictor-corrector, filtering, terminal-bank, and coefficient
algebra without inventing the omitted scheduling and controller data.

Eqs. 17-22 are available as separate helpers:

```matlab
Rref_printed = paperstudies.saito.eq17_reference_range(LD, V, g, rE);
Dref = paperstudies.saito.eq18_reference_drag(LD, V, g, rE);
hdot_ref = paperstudies.saito.eq19_reference_altitude_rate(LD, V0, V);
F1 = paperstudies.saito.eq20_gain_F1(K1, V, g, rE);
F2 = paperstudies.saito.eq21_gain_F2(K2, Rref, LD);
F3 = paperstudies.saito.eq22_gain_F3(K3, Rref, LD);
```

Their algebra is exact, but their numerical use is
`CONVENTION_REQUIRED`: Eq. 17 does not print the redimensionalization needed
to obtain the km values in Table 5, and the paper does not publish
`V0`, `K1`, `K2`, or `K3`. `eq17_reference_range` therefore returns the
expression exactly as printed and does not silently multiply by an assumed
Earth radius.

An optional reduced trajectory exercises the shared project propagator:

```matlab
result = paperstudies.saito.run(struct( ...
    'forward', true, ...
    'entry_case', 1));
```

That path uses the published 150 kg mass, 0.554 m^2 area, and baseline
`Cd=1.3`, but assumes constant `L/D=0.25` and zero-degree trim. It is always
reported as `SURROGATE_ONLY` and is not a reconstruction of the paper's
landing footprint, guidance performance, or heating.

Full numerical reproduction remains blocked by the unavailable proprietary
HSRC `CA/CN(Mach,AoA)` database, Mach-10 trim angle, cross-range deadband,
explicit-guidance `V0/K1/K2/K3` values and distance scaling, fading-filter
coefficient, predictor settings, terminal-bank transition, RW/RCS control
implementation, the 27 individual Table-6 entry states, and complete
Detra-Kemp-Riddell inputs. The
package also reports, rather than fits away, the roughly 63-65 km inconsistency
between Table-4 coordinates and its printed spherical range-to-go values.
