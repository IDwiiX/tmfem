# tmfem

A from-scratch C++17 finite-element solver for 2D transient heat conduction, linear elasticity,
von Mises elastoplasticity and one-way thermo-mechanical coupling. Every capability is checked by a
validation test against an analytical solution or a reference result. It is a learning and validation
project; it does not replace ANSYS or Abaqus.

**Status:** phase 1 of 5 validated: 1D transient heat conduction and diffusion. 2D elasticity,
plasticity and thermo-mechanical coupling are not built yet.

## Results

Every number below is written by the test that asserts it, into `results/validation/`, and the
table is regenerated from those files by `python3 scripts/build_results_table.py` after `ctest`.
Reference IDs and tolerances are defined in [docs/REFERENCE_VALUES.md](docs/REFERENCE_VALUES.md).
The error column is relative (%) except where the tolerance is absolute or the reference is zero.

### Phase 1: 1D heat and diffusion

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
```

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

## Credits

Benchmark cases adapted from M. R. Kabiri, ENSAM Meknès lab materials (with permission);
analytical solutions from the textbooks cited next to each test.
