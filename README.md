# tmfem

A from-scratch C++17 finite-element solver for 2D transient heat conduction, linear elasticity,
von Mises elastoplasticity and one-way thermo-mechanical coupling. Every capability is checked by a
validation test against an analytical solution or a reference result. It is a learning and validation
project; it does not replace ANSYS or Abaqus.

**Status:** phases 1 to 3 of 5 validated:
- 1D transient heat conduction and diffusion
- 2D linear elasticity with T3/T6 triangles, plane stress and plane strain
- small-strain von Mises plasticity with isotropic hardening (tabulated and Ludwik): radial return,
  consistent tangent, global Newton with load stepping, plane strain and plane stress

Thermo-mechanical coupling is not built yet.

## Results

Every number below is written by the test that asserts it, into `results/validation/`, and the
table is regenerated from those files by `python3 scripts/build_results_table.py` after `ctest`.
Reference IDs and tolerances are defined in [docs/REFERENCE_VALUES.md](docs/REFERENCE_VALUES.md).
The error column is relative (%) except where the tolerance is absolute or the reference is zero.
Report-only rows are compared with a value that is not a target, such as an estimate or the source's
parameters. Their difference is shown in parentheses, and "n/a" marks rows with no comparison value.

### Phases 1 to 3

<!-- results:begin -->
| ID | Benchmark | Quantity | Reference | tmfem | Error | Tolerance | Test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1-a | Steady diffusion flux (ex. 1) | J finite element | 2.4e-09 | 2.4e-09 | +0% | 1e-12 relative | `P1-a steady flux ex1` |
| P1-a | Steady diffusion flux (ex. 1) | J formula | 2.4e-09 | 2.4e-09 | +0% | 1e-12 relative | `P1-a steady flux ex1` |
| P1-b | Carburizing, erf solution (ex. 2) | time analytic [s] | 25370 | 25370.2 | +0.000968% | 0.5% | `P1-b carburizing time to 0.8 wt% at 0.5 mm` |
| P1-b | Carburizing, erf solution (ex. 2) | time finite element [s] | 25370 | 25373.7 | +0.0144% | 0.5% | `P1-b carburizing time to 0.8 wt% at 0.5 mm` |
| P1-c | Equivalent time (ex. 3) | t at 500 C [h] | 110.4 | 110.417 | +0.0151% | 0.001 relative | `P1-c equivalent time at 500 C` |
| P1-d | Arrhenius D(T) (ex. 4) | D(1000 K) [m2/s] | 2.067e-13 | 2.0669e-13 | -0.00461% | 0.005 relative | `P1-d Arrhenius D(1000 K), D(800 K), x at 800 K` |
| P1-d | Arrhenius D(T) (ex. 4) | D(800 K) [m2/s] | 2.14e-15 | 2.1396e-15 | -0.0187% | 0.005 relative | `P1-d Arrhenius D(1000 K), D(800 K), x at 800 K` |
| P1-d | Arrhenius D(T) (ex. 4) | x at 800 K [mm] | 1.577 | 1.57702 | +0.00138% | 0.005 relative | `P1-d Arrhenius D(1000 K), D(800 K), x at 800 K` |
| P1-e | Boron predeposition dose (ex. 5) | dose finite element [at/cm2] | 4.062e+15 | 4.06214e+15 | +0.00336% | 0.005 relative | `P1-e boron predeposition dose` |
| P1-e | Boron predeposition dose (ex. 5) | dose formula [at/cm2] | 4.062e+15 | 4.06217e+15 | +0.00406% | 0.005 relative | `P1-e boron predeposition dose` |
| P1-f | Thin-film surface concentration (ex. 5) | surface C at 1200 C [at/cm3] | 6.97e+19 | 6.97382e+19 | +0.0548% | 0.005 relative | `P1-f thin-film surface concentration at 1200 C` |
| P1-g | Arrhenius fit (ex. 5) | D(1300 C) [cm2/s] | 1.74e-11 | 1.74154e-11 | +0.0884% | 0.005 relative | `P1-g Arrhenius fit Q, D0, D(1300 C)` |
| P1-g | Arrhenius fit (ex. 5) | D0 [cm2/s] | 3.11 | 3.1119 | +0.0609% | 0.005 relative | `P1-g Arrhenius fit Q, D0, D(1300 C)` |
| P1-g | Arrhenius fit (ex. 5) | Q [kJ/mol] | 338.9 | 338.867 | -0.00977% | 0.005 relative | `P1-g Arrhenius fit Q, D0, D(1300 C)` |
| P1-h | Thin film on a half-space | max relative dose drift | 0 | 3.69059e-12 | +3.69e-12 | 1e-10 | `P1-h thin film dose conservation` |
| P1-h | Thin film on a half-space | peak decay slope | -0.5 | -0.499996 | +4.1e-06 | +/- 0.01 | `P1-h thin film peak decays as t^-1/2` |
| P1-i | Steel plate parameters | Biot number | 0.00183 | 0.00183333 | +3.33e-06 | +/- 5e-6 (last quoted digit) | `P1-i steel plate Biot number and time constant` |
| P1-i | Steel plate parameters | time constant [s] | 549.5 | 549.5 | +0 | +/- 0.05 s (last quoted digit) | `P1-i steel plate Biot number and time constant` |
| P1-j | Steel plate, FE vs lumped ODE | max abs inner face FE - lumped [K] | 0 | 0.0330652 | +0.0331 | 0.05 K | `P1-j steel plate 1D FE vs lumped ODE` |
| P1-k | Steel plate, through-thickness difference | max through-thickness difference [K] | 0.0255348 | 0.0255832 | +0.19% | 5% | `P1-k steel plate through-thickness difference` |
| P1-m | 1D heat convergence orders | forced time order BE consistent mass | 1 | 0.998348 | -0.00165 | +/- 0.15 | `P1-m temporal order with time-dependent forcing, backward Euler` |
| P1-m | 1D heat convergence orders | forced time order BE lumped mass | 1 | 0.997965 | -0.00204 | +/- 0.15 | `P1-m temporal order with time-dependent forcing, backward Euler` |
| P1-m | 1D heat convergence orders | forced time order CN consistent mass | 2 | 2.0002 | +0.0002 | +/- 0.15 | `P1-m temporal order with time-dependent forcing, Crank-Nicolson` |
| P1-m | 1D heat convergence orders | forced time order CN lumped mass | 2 | 2.00019 | +0.000191 | +/- 0.15 | `P1-m temporal order with time-dependent forcing, Crank-Nicolson` |
| P1-m | 1D heat convergence orders | space order consistent mass | 2 | 1.99561 | -0.00439 | +/- 0.15 | `P1-m spatial order consistent mass` |
| P1-m | 1D heat convergence orders | space order lumped mass | 2 | 2.00849 | +0.00849 | +/- 0.15 | `P1-m spatial order lumped mass` |
| P1-m | 1D heat convergence orders | time order BE | 1 | 0.98331 | -0.0167 | +/- 0.15 | `P1-m temporal order backward Euler` |
| P1-m | 1D heat convergence orders | time order CN | 2 | 2.00048 | +0.000483 | +/- 0.15 | `P1-m temporal order Crank-Nicolson` |
| P2-a | Patch test, distorted mesh | T3 max relative error (u, sigma) | 0 | 5.81102e-15 | +5.81e-15 | 1e-10 relative | `P2-a patch test T3 distorted mesh` |
| P2-a | Patch test, distorted mesh | T6 max relative error (u, sigma) | 0 | 7.96939e-14 | +7.97e-14 | 1e-10 relative | `P2-a patch test T6 distorted mesh` |
| P2-b | Uniform bar in tension | T3 max rel error of end displacement | 0 | 1.61639e-14 | +1.62e-14 | 1e-8 relative | `P2-b uniform bar end displacement` |
| P2-b | Uniform bar in tension | T6 max rel error of end displacement | 0 | 1.22986e-13 | +1.23e-13 | 1e-8 relative | `P2-b uniform bar end displacement` |
| P2-c | Kirsch plate, stress concentration | Kt finest T3 (report only) | 3 | 2.97304 | -0.899% | not asserted | `P2-c Kirsch Kt on finest T6 mesh` |
| P2-c | Kirsch plate, stress concentration | Kt finest T6 | 3 | 2.99969 | -0.0103% | 1% | `P2-c Kirsch Kt on finest T6 mesh` |
| P2-d | Kirsch plate, r/a = 2 and 4 | max abs stress error, 20 rows x 3 components [MPa] | 0 | 0.00441135 | +0.00441 | 2 MPa | `P2-d Kirsch stresses at r/a = 2 and 4` |
| P2-e | Kirsch plate, r/a = 1.5 | max abs stress error, 10 rows x 3 components [MPa] | 0 | 0.0117644 | +0.0118 | 2 MPa | `P2-e Kirsch stresses at r/a = 1.5` |
| P2-f | Kirsch convergence orders | T3 plane strain L2 order (levels 1-4) | 2 | 1.96584 | -0.0342 | +/- 0.15 | `P2-f Kirsch convergence T3` |
| P2-f | Kirsch convergence orders | T3 plane strain energy order (levels 1-4) | 1 | 0.990579 | -0.00942 | +/- 0.15 | `P2-f Kirsch convergence T3` |
| P2-f | Kirsch convergence orders | T3 plane stress L2 order (levels 1-4) | 2 | 1.97459 | -0.0254 | +/- 0.15 | `P2-f Kirsch convergence T3` |
| P2-f | Kirsch convergence orders | T3 plane stress energy order (levels 1-4) | 1 | 0.992547 | -0.00745 | +/- 0.15 | `P2-f Kirsch convergence T3` |
| P2-f | Kirsch convergence orders | T6 plane strain L2 order (levels 1-4) | 3 | 3.12065 | +0.121 | +/- 0.15 | `P2-f Kirsch convergence T6` |
| P2-f | Kirsch convergence orders | T6 plane strain energy order (levels 1-4) | 2 | 1.96535 | -0.0346 | +/- 0.15 | `P2-f Kirsch convergence T6` |
| P2-f | Kirsch convergence orders | T6 plane stress L2 order (levels 1-4) | 3 | 3.10864 | +0.109 | +/- 0.15 | `P2-f Kirsch convergence T6` |
| P2-f | Kirsch convergence orders | T6 plane stress energy order (levels 1-4) | 2 | 1.97369 | -0.0263 | +/- 0.15 | `P2-f Kirsch convergence T6` |
| P3-a | Ludwik fit to the true tensile curve | LSQ K [MPa], 95% CI +/- 39.23, jackknife SE 49.82 (reference = source) | 559.27 | 520.875 | (-6.87%) | report only | `P3-a Ludwik least-squares fit with uncertainty` |
| P3-a | Ludwik fit to the true tensile curve | LSQ m, 95% CI +/- 0.03985, jackknife SE 0.1269 (reference = source) | 0.4077 | 0.361757 | (-11.3%) | report only | `P3-a Ludwik least-squares fit with uncertainty` |
| P3-a | Ludwik fit to the true tensile curve | LSQ residual standard deviation s [MPa] | n/a | 4.57928 | n/a | report only | `P3-a Ludwik least-squares fit with uncertainty` |
| P3-a | Ludwik fit to the true tensile curve | LSQ sigma_y0 [MPa], 95% CI +/- 11.66, jackknife SE 101.7 (reference = source) | 310.862 | 309.135 | (-0.555%) | report only | `P3-a Ludwik least-squares fit with uncertainty` |
| P3-a | Ludwik fit to the true tensile curve | log-log K [MPa] (sigma_y fixed) | 548 | 547.608 | -0.0715% | 0.5% | `P3-a Ludwik log-log fit with fixed sigma_y` |
| P3-a | Ludwik fit to the true tensile curve | log-log m (sigma_y fixed) | 0.385 | 0.384864 | -0.0353% | 0.5% | `P3-a Ludwik log-log fit with fixed sigma_y` |
| P3-a | Ludwik fit to the true tensile curve | source m vs our log-log m (report only) | 0.4077 | 0.384864 | (-5.6%) | report only | `P3-a Ludwik log-log fit with fixed sigma_y` |
| P3-b | Material point, uniaxial, tabulated hardening | eps_xx at 457.24 MPa | -0.01535 | -0.0153496 | +0.00246% | 0.5% | `P3-b uniaxial strains at 457.24 MPa, tabulated hardening` |
| P3-b | Material point, uniaxial, tabulated hardening | eps_yy at 457.24 MPa | 0.03158 | 0.031581 | +0.00317% | 0.5% | `P3-b uniaxial strains at 457.24 MPa, tabulated hardening` |
| P3-b | Material point, uniaxial, tabulated hardening | eps_zz at 457.24 MPa | -0.01535 | -0.0153496 | +0.00246% | 0.5% | `P3-b uniaxial strains at 457.24 MPa, tabulated hardening` |
| P3-c | Material point, Ludwik along the true curve | max |driver - data| / s, LSQ Ludwik, 7 points | 0 | 1.5942 | +1.59 | +/- 2 s | `P3-c Ludwik material point along the true curve` |
| P3-d | Newton convergence rate | Newton rate, last qualifying triple, plane-strain bar, Ludwik | 2 | 2.12018 | +6.01% | >= 1.8 | `P3-d Newton convergence rate on the plane-strain bar` |
| P3-d | Newton convergence rate | Newton rate, last qualifying triple, plane-stress bar, Ludwik | 2 | 2.61767 | +30.9% | >= 1.8 | `P3-d Newton convergence rate on the plane-stress bar` |
| P3-e | Plane-strain bar FE vs material point | max relative difference FE vs driver (T3, T6, tabulated, Ludwik) | 0 | 1.11596e-12 | +1.12e-12 | 1e-6 relative | `P3-e plane-strain bar FE equals material-point driver` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 3: hole-edge eps_yy | n/a | 0.0019402 | n/a | report only | `P3-f elastoplastic hole plate, plane strain (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 3: max alpha (at the edge) | n/a | 0.000407141 | n/a | report only | `P3-f elastoplastic hole plate, plane strain (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 3: max sigma_yy in the field [MPa], at r = 1.074 a | n/a | 371.569 | n/a | report only | `P3-f elastoplastic hole plate, plane strain (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 3: plastic zone extent (r - a) / a | n/a | 0.0735324 | n/a | report only | `P3-f elastoplastic hole plate, plane strain (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 3: sigma_yy at the point nearest the edge [MPa] | n/a | 358.096 | n/a | report only | `P3-f elastoplastic hole plate, plane strain (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 4: hole-edge eps_yy | n/a | 0.00193939 | n/a | report only | `P3-f elastoplastic hole plate, plane strain, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 4: max alpha (at the edge) | n/a | 0.000411336 | n/a | report only | `P3-f elastoplastic hole plate, plane strain, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 4: max sigma_yy in the field [MPa], at r = 1.075 a | n/a | 371.746 | n/a | report only | `P3-f elastoplastic hole plate, plane strain, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 4: plastic zone extent (r - a) / a | n/a | 0.0748925 | n/a | report only | `P3-f elastoplastic hole plate, plane strain, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane strain, T6 level 4: sigma_yy at the point nearest the edge [MPa] | n/a | 357.93 | n/a | report only | `P3-f elastoplastic hole plate, plane strain, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 3: hole-edge eps_yy (reference = Neuber) | 0.00264675 | 0.0023565 | (-11%) | report only | `P3-f elastoplastic hole plate, plane stress (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 3: max alpha (at the edge) (reference = Neuber) | 0.00109767 | 0.000811179 | (-26.1%) | report only | `P3-f elastoplastic hole plate, plane stress (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 3: max sigma_yy in the field [MPa], at r = 1.175 a | n/a | 331.13 | n/a | report only | `P3-f elastoplastic hole plate, plane stress (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 3: plastic zone extent (r - a) / a | n/a | 0.174643 | n/a | report only | `P3-f elastoplastic hole plate, plane stress (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 3: sigma_yy at the point nearest the edge [MPa] (reference = Neuber) | 321.314 | 318.749 | (-0.798%) | report only | `P3-f elastoplastic hole plate, plane stress (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 4: hole-edge eps_yy (reference = Neuber) | 0.00264675 | 0.00235585 | (-11%) | report only | `P3-f elastoplastic hole plate, plane stress, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 4: max alpha (at the edge) (reference = Neuber) | 0.00109767 | 0.000815221 | (-25.7%) | report only | `P3-f elastoplastic hole plate, plane stress, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 4: max sigma_yy in the field [MPa], at r = 1.174 a | n/a | 331.231 | n/a | report only | `P3-f elastoplastic hole plate, plane stress, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 4: plastic zone extent (r - a) / a | n/a | 0.174059 | n/a | report only | `P3-f elastoplastic hole plate, plane stress, level 4 (report)` |
| P3-f | Elastoplastic hole plate, S = 140 MPa (report only) | plane stress, T6 level 4: sigma_yy at the point nearest the edge [MPa] (reference = Neuber) | 321.314 | 318.675 | (-0.821%) | report only | `P3-f elastoplastic hole plate, plane stress, level 4 (report)` |
<!-- results:end -->

Notes:
- **P1-k:** the original target "below 0.01 K" is physically infeasible for the 26 to 74 °C test
  signal. It holds only for swings under 19 K. With the user's approval it was replaced by the
  quasi-steady bound `(L/k)/R_tot · max|T_sa − T_in| = 0.0255 K` at 5%. The derivation is in
  [docs/learning/phase1.md](docs/learning/phase1.md) section 9.
- **P1-j:** 0.033 K at the prescribed dt = 300 s (0.026 K at dt = 150 s) contains two parts:
  - backward-Euler time error, visible as an hourly sawtooth from the hourly-sampled input
  - the physical difference between the inner face and the plate mean, of order the Biot number
- **P1-j and P1-k use a synthetic sol-air signal** `T_sa = 50 − 24 cos(2πt/day)` °C. The comparison
  with ANSYS (P1-l) is skipped until the real input file and the ANSYS curve are available.
- **Halved-resolution reruns:** the P1-b, P1-e, P1-h and P1-j tests also run at half the mesh size
  and/or half the time step, and assert that the result moves by less than 1/5 of the tolerance.
- **P1-m measures the observed orders on four refinements.** The time order is measured both on an
  unforced eigenmode problem and on a forced problem with a source, time-dependent Dirichlet and
  Robin conditions.

| | |
| --- | --- |
| ![Carburizing: FE vs erf solution](docs/figures/phase1_erf_profile.png) | ![Thin-film decay and dose conservation](docs/figures/phase1_thin_film_decay.png) |
| ![Steel plate inner face vs lumped model](docs/figures/phase1_steel_plate_inner_face.png) | ![Convergence orders](docs/figures/phase1_convergence.png) |

Notes on phase 2 (derivations in [docs/learning/phase2.md](docs/learning/phase2.md)):
- **Kirsch set-up:** a quarter model [0, 5a]² with symmetry conditions. The exact Kirsch tractions are
  applied on the outer edges, so the finite model represents the infinite plate exactly (section 9).
  Remote stress 140 MPa along y, plane stress, unit thickness. Each point stress is computed in the
  element that contains the point, and averaged over all such elements when the point lies on a shared
  edge or node (the Kt node (a, 0), for example). Nodal averages are used only for plots and the VTU file.
- **Meshes:** generated by Gmsh 4.12.1 from `cases/kirsch/kirsch.geo` and committed. There are five
  levels (`-clscale` 1 to 1/16), with curved T6 edges on the hole. "Finest" means level 4.
- **P2-f, disclosed change:** the first run fitted levels 0 to 3, with level 4 not yet generated.
  - The T6 L2 order came out at 3.25 (plane stress) and 3.27 (plane strain), outside 3 ± 0.15.
  - The pairwise orders, 3.53, 3.17 and 3.12, were falling toward 3. A scratch level 4 gave 3.03.
  - So level 0 is pre-asymptotic. The error-norm quadrature was checked and is not the cause.
  - Adding level 4 alone is not enough: a fit over all five levels gives 3.19 (plane stress) and 3.21
    (plane strain), still outside the tolerance. The pass depends on leaving level 0 out of the fit.
  - With the user's approval, level 4 was added and the order is now fitted on levels 1 to 4, for T3
    and T6 alike. Level 0 is still computed and printed by the test.
- **P2-c..e on level 3**, the finest mesh before level 4 was added: Kt = 2.99876 (−0.041%), and the
  largest point-stress errors were 0.036 MPa (P2-d) and 0.031 MPa (P2-e).
- **P2-c also reports T3:** 2.973 on level 4 (−0.9%). It is not asserted, because the reference row is
  for T6.

| | |
| --- | --- |
| ![Kirsch plate von Mises stress](docs/figures/phase2_kirsch_von_mises.png) | ![Kirsch profiles vs exact](docs/figures/phase2_kirsch_profiles.png) |
| ![Kirsch convergence orders](docs/figures/phase2_convergence.png) | |

Notes on phase 3 (derivations in [docs/learning/phase3.md](docs/learning/phase3.md)):
- **Material data:**
  - E = 207.423 GPa and ν = 0.3.
  - The tabulated hardening curve is `references/data/tensile_true.csv` from the yield point
    (310.56 MPa) on, with α = ε − σ/E. The first α comes out at −3.2e-6 because the digitised elastic
    slope is 210 GPa, and is clipped to 0.
  - The Ludwik source parameters are 310.86 MPa, 559.27 MPa and 0.4077.
- **P3-a, the fit:**
  - The log-log fit with σ_y fixed reproduces the reference (m = 0.385, K = 548 MPa).
  - The deliverable is a 3-parameter least-squares fit on stress: σ_y0 = 309.1 ± 11.7 MPa,
    K = 520.9 ± 39.2 MPa, m = 0.362 ± 0.040 (95% intervals, covariance, 5 dof). The source values are
    in the table next to these.
  - The leave-one-out jackknife gives much larger standard errors (102 MPa, 50 MPa, 0.13).
    Without the yield point σ_y0 drops to 193 MPa: no data lie between εᵖ = 0 and 0.0086, so one point
    fixes the knee. That is a limit of the 8-point data set, not of the solver.
- **P3-b is nearly tautological.** 457.24 MPa is a data point, and the tabulated law passes through it,
  so ε_yy is a round trip. The informative check is the lateral strain: −νσ/E elastic plus −α/2
  plastic (incompressibility).
- **P3-d, disclosed change.**
  - The plan asserted the Newton rate on the Kirsch plate (T6 level 2, plane strain, one step from 100
    to 140 MPa).
  - With a consistent tangent, Newton went from 1.9e-6 to 1.4e-10 to the 7e-14 round-off floor once
    the plastic set settled. That left too few iterates for the pre-fixed rate criterion: no qualifying
    triple.
  - With the user's approval, the rate is asserted on the FE bar, where the active set is fixed by
    construction. The criterion is unchanged. The hole-plate sequence is still printed and plotted.
  - The bar is uniform, so this rate checks the tangent along the loading path: it amounts to a scalar
    Newton on σ_yy(ε_yy). It still fails with a continuum tangent.
  - The full global tangent is checked separately. The test "Global tangent matches finite differences of
    the internal force" compares K_T·v with a central difference of f_int on a non-uniform plastic
    Kirsch state, to 1e-6.
- **P3-c, the ±2s band is a weak test for the least-squares fit:** with 8 points and 3 parameters, no
  residual can exceed √5·s = 2.24s. The source's Ludwik parameters miss the band at 3 of the 7 points
  (−2.95s, −2.99s, −2.26s). This is printed by the test, not asserted.
- **P3-f is a report, not a validation.**
  - Set-up: the quarter Kirsch model with the elastic Kirsch tractions on the 5a boundary (an
    approximation once a plastic zone exists), tabulated hardening, and 28 load steps.
  - Plane stress is the thin plate of the source's Abaqus model. On T6 level 3, at the notch root, FE
    against Neuber's estimate gives:

    | Quantity | FE | Neuber | Difference |
    | --- | --- | --- | --- |
    | edge strain ε_yy | 0.002357 | 0.002647 | −11% |
    | equivalent plastic strain | 0.000811 | 0.001098 | −26% |
    | stress at the point nearest the edge | 318.7 MPa | 321.3 MPa | −0.8% |

    The level-4 values are in the table. Neuber overestimating the strain is the usual direction.
  - The largest σ_yy in the field (331 MPa) is not at the root: it sits at the plastic-zone boundary,
    r ≈ 1.17a, where σ_xx > 0.
  - The plastic-zone extent is 0.159a, 0.167a, 0.175a and 0.174a on levels 1 to 4: +4.6% from level 2 to 3,
    then −0.3% from level 3 to 4.
  - In plane strain the notch root is biaxial (σ_zz ≈ 124 MPa), so uniaxial Neuber does not apply and
    is listed only for comparison.
  - The level-4 rows come from tests tagged `slow`.

| | |
| --- | --- |
| ![Ludwik fits and the uniaxial material point](docs/figures/phase3_tensile_fit.png) | ![Newton residual histories](docs/figures/phase3_newton_history.png) |
| ![Plastic zone at the hole edge](docs/figures/phase3_plastic_zone.png) | |

## Build and test

Requires a C++17 compiler, CMake >= 3.16 and Eigen 3 (`sudo apt install libeigen3-dev`).
Catch2 and nlohmann/json are downloaded at pinned versions during configuration.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run one phase with `ctest --test-dir build -L phase1`, or one reference row with
`ctest --test-dir build -L P2-c`.

Run a case from the repository root (paths inside case files are relative to it):

```bash
./build/src/tmfem cases/steel_plate/synthetic.json   # writes results/steel_plate/synthetic_faces.csv
python3 scripts/plot_phase1.py                       # phase 1 figures, after ctest -L phase1
./build/src/tmfem cases/kirsch/kirsch_t6.json        # writes results/kirsch/kirsch_t6.vtu (ParaView)
python3 scripts/plot_phase2.py                       # phase 2 figures, after the case and ctest -L phase2
./build/src/tmfem cases/kirsch/kirsch_plastic_plane_stress.json  # elastoplastic plate: VTU + quadrature-point CSV
python3 scripts/plot_phase3.py                       # phase 3 figures, after ctest -L phase3
```

The two level-4 elastoplastic plate tests take a few minutes each. Skip them with `ctest -LE slow`.

The `.msh` files are committed. `bash scripts/make_meshes.sh` regenerates them and needs Gmsh 4.12.1.

## Limitations

Out of scope by design: 3D, large deformation, dynamics, temperature-dependent properties,
contact, radiation, real particle-shower energy maps.

Phase 1 (1D heat and diffusion):
- **Linear elements, one dimension:** lateral faces are adiabatic.
- **Constant properties:**
  - k, ρc and D are constant in time and temperature.
  - The Arrhenius `D(T)` is evaluated, not coupled into the solver.
- **Linear convection only, no radiation.**
- **Fixed time step:** changing dt refactorises the matrix, and there is no adaptive error control.
- **The steel-plate validation is against a lumped model with a synthetic input.** It is not yet
  compared with the source's ANSYS run.

Phase 2 (2D linear elasticity):
- **Small strain, linear isotropic material, one material per run, no body force.**
- **Unit thickness:** reactions are in N per metre of depth.
- **Triangles only (T3, T6):** no quadrilaterals. On curved T6 the 3-point stiffness rule is not exact,
  but it keeps the optimal convergence order.
- **Point stresses come straight from the elements, averaged where elements meet:** there is no stress recovery scheme beyond nodal
  averaging for plots.
- **The Kirsch benchmark uses the infinite-plate tractions:** it is not a finite-width plate like the
  source's Abaqus model.
- **Mesh dependence:** another Gmsh version would produce slightly different meshes, and so slightly
  different numbers.

Phase 3 (elastoplasticity):
- **Small strain, isotropic hardening only:** there is no kinematic hardening, so no Bauschinger effect.
  The model is rate-independent and isothermal.
- **The hole plate uses elastic Kirsch tractions** on the 5a boundary. Hole-edge stresses come from the
  nearest quadrature point, and the edge strain comes from the displacement field.
- **Tabulated hardening has slope jumps** at the data points, so quadratic Newton convergence is
  asserted with Ludwik only.
- **The fit uncertainty describes 8 digitised points and model misfit,** not material scatter.

## Credits

Benchmark cases adapted from M. R. Kabiri, ENSAM Meknès lab materials (with permission);
analytical solutions from the textbooks cited next to each test.
