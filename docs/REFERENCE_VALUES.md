# Reference values for validation

Authoritative test targets. Every value was recomputed by
`references/verify_reference_values.py`; where the source document disagrees, the
"Source says" column shows it and the recomputed value is the target. Never adjust a
test to match a "Source says" value that is marked wrong.

## Phase 1: heat and diffusion

| ID | Case | Quantity | Target | Source says | Tolerance |
| --- | --- | --- | --- | --- | --- |
| P1-a | Diffusion ex. 1 (unit test) | Steady flux J | 2.4e-9 kg/(m2 s) | 2.4e-12 (wrong, mm not converted) | 1e-12 relative |
| P1-b | Carburizing ex. 2, C0 = 0.25, Cs = 1.2 wt%, D = 1.6e-11 m2/s | Time for 0.8 wt% at 0.5 mm | 25 370 s (7.047 h) | 7.025 h (erf table interpolation) | 0.5% |
| P1-c | Ex. 3 (unit test) | Equivalent time at 500 C | 110.4 h | 110.4 h | 0.1% |
| P1-d | Ex. 4, D0 = 1.8e-5 m2/s, Q = 152 kJ/mol (unit test) | D(1000 K), D(800 K), x at 800 K | 2.067e-13 m2/s, 2.140e-15 m2/s, 1.577 mm | 2.04e-13, 2.11e-15, 1.575 mm | 0.5% |
| P1-e | Ex. 5, boron in Si, t = 360 s | Dose M | 4.062e15 at/cm2 | 4.06e15 | 0.5% |
| P1-f | Ex. 5, thin film at 1200 C | Surface concentration | 6.97e19 at/cm3 | 2.8e19 (matches t = 36 min) | 0.5% |
| P1-g | Ex. 5, Arrhenius fit | Q, D0, D(1300 C) | 338.9 kJ/mol, 3.11 cm2/s, 1.74e-11 cm2/s | 339, 3.2 m2/s (unit slip), 5.86e-11 m2/s | 0.5% |
| P1-h | Thin film, numerical | Dose conservation; peak ~ t^-1/2 | exact | n/a | 1e-10 relative; slope -0.5 +/- 0.01 |
| P1-i | Steel plate (steel_plate_case.json) | Biot number; time constant | 0.00183; 549.5 s | not given | exact formulas |
| P1-j | Steel plate, synthetic T_sa | 1D FE vs lumped ODE (RK4, small step) | agreement | n/a | max 0.05 K |
| P1-k | Steel plate, synthetic T_sa | Through-thickness temperature difference | below 0.01 K | n/a | bound |
| P1-l | Steel plate, real ttab (when provided) | Inner-face peaks; RMS vs ANSYS | report | about 51, 53, 60 C | report only, no tuning |
| P1-m | Any analytic case | Observed order: BE time, CN time, space | 1, 2, 2 | n/a | +/- 0.15 |

Physics note for P1-l: with a 9-minute time constant the plate cannot carry heat from one
day to the next, so rising daily peaks must come from the input signal. The source report's
explanation (accumulated heat, steel's high conductivity) does not hold; conduction is
0.05% of the total thermal resistance.

## Phase 2: linear elasticity

| ID | Case | Quantity | Target | Source says | Tolerance |
| --- | --- | --- | --- | --- | --- |
| P2-a | Patch test, distorted mesh, T3 and T6 | Linear field reproduced | exact | n/a | 1e-10 relative |
| P2-b | Uniform bar | End displacement | sigma L / E | n/a | 1e-8 relative |
| P2-c | Kirsch, S = 140 MPa, exact outer tractions | sigma_tt at r = a, theta = 90 deg (from load axis) | 420 MPa (Kt = 3) | 420 MPa | 1% on finest T6 mesh |
| P2-d | Kirsch | All rows at r/a = 2 and 4 in kirsch_sigma140MPa.csv | CSV values | match | 2 MPa |
| P2-e | Kirsch | Rows at r/a = 1.5 | CSV values (e.g. theta = 0: 25.93, -10.37, 0) | 2100, -3080, 0 (wrong) | 2 MPa |
| P2-f | Kirsch | Convergence: displacement L2 order; energy-norm order | T3: 2, 1; T6: 3, 2 | n/a | +/- 0.15 |

Angle convention: theta is measured from the loading axis y, as in the source. At r = a,
sigma_tt = -S at theta = 0 and 3S at theta = 90 deg.

## Phase 3: elastoplasticity

Data: `references/data/tensile_true.csv`, `material_tensile.json`. E = 207.423 GPa, nu = 0.3.

| ID | Case | Quantity | Target | Source says | Tolerance |
| --- | --- | --- | --- | --- | --- |
| P3-a | Ludwik fit, sigma_y = 310.56 MPa fixed, log-log | m, K | 0.385, 548 MPa | 0.4077, 559.27 MPa (different point choice) | report both; least-squares fit with uncertainty is the deliverable |
| P3-b | Material point, uniaxial stress, tabulated hardening | Strains at sigma = 457.24 MPa | eps_yy 0.03158, eps_xx = eps_zz -0.01535 | -0.01535, 0.03158, -0.01535 (Abaqus) | 0.5% |
| P3-c | Material point, Ludwik hardening | Stress along the true curve, eps_p > 0 | data points | data | within the fit's residual band |
| P3-d | Newton iterations | Observed convergence rate | quadratic | n/a | rate >= 1.8 |
| P3-e | Plane-strain uniform bar | FE vs material-point driver with eps_zz = 0 | identical | n/a | 1e-6 relative |
| P3-f | Kirsch plate at S = 140 MPa, elastoplastic | Peak stress and strain at hole edge; plastic zone size; Neuber estimate | report | Abaqus von Mises map only | report only |

## Phase 4: thermo-mechanical

| ID | Case | Quantity | Target | Tolerance |
| --- | --- | --- | --- | --- |
| P4-a | Bar in plane stress, held in x, free laterally, uniform dT | sigma_xx, sigma_yy | -E alpha dT, 0 | 1e-8 relative |
| P4-b | Free body, uniform dT | Stress; displacement | 0; alpha dT x | 1e-10 relative |
| P4-c | Thin disk, T = dT (1 - r^2/b^2), plane stress | sigma_r / (E alpha dT) at r/b = 0, 0.5, 1 | -0.25, -0.1875, 0 | 1% |
| P4-d | Same | sigma_theta / (E alpha dT) at r/b = 0, 0.5, 1 | -0.25, -0.0625, 0.5 | 1% (absolute 0.0025 near zero crossings) |
| P4-e | 2D Gaussian heat spreading | Variance s0^2 + 2 alpha t; energy | closed form | 0.5%; 1e-10 relative |

P4-c and P4-d come from sigma_r = (E alpha / b^2) int_0^b T r dr - (E alpha / r^2) int_0^r T r dr
and sigma_theta = (E alpha / b^2) int_0^b T r dr + (E alpha / r^2) int_0^r T r dr - E alpha T
(Timoshenko and Goodier), checked numerically by the verification script.

## Stretch S1: cracked strip

See `references/data/crack_case.md`. Target: K_I = sigma sqrt(pi a) / sqrt(cos(pi a / (2b)))
for the chosen a/b, within 2% with the J-integral. The source's 840 vs 922 MPa sqrt(m) cannot
be rechecked until its geometry is known.
