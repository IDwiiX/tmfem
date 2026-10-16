# tmfem specification

A from-scratch C++17 finite-element solver for 2D thermo-mechanics, built in five phases.
Each phase ends with passing validation tests. Reference numbers live in
`docs/REFERENCE_VALUES.md`; source data lives in `references/data/`.

## 1. Context

CERN designs beam-intercepting devices (targets, collimators, beam dumps). Their workflow:
a Monte Carlo code (FLUKA) computes where the beam deposits energy, then a finite-element
code (ANSYS) computes the temperature field and the resulting thermal stress. tmfem
rebuilds the second step from scratch, validated against analytical solutions and against
Abaqus/ANSYS results from M. R. Kabiri's ENSAM Meknès lab materials. It does not replace
ANSYS and must never claim to.

## 2. Scope

In scope: 1D and 2D transient heat conduction and diffusion; 2D linear elasticity (plane
stress, plane strain); small-strain von Mises plasticity with isotropic hardening;
one-way thermal-to-mechanical coupling; VTU output for ParaView.

Out of scope: 3D, large deformation, dynamics and stress waves from fast heating,
temperature-dependent properties, contact, radiation, real particle-shower energy maps.
The report states these limits explicitly.

## 3. Architecture

- `mesh`: read Gmsh `.msh` 4.1 ASCII (nodes, elements, physical groups); element geometry.
- `fem`: shape functions (1D linear, T3, T6), Gauss quadrature, Jacobians, sparse assembly
  via Eigen triplets, Dirichlet elimination, Neumann/Robin boundary integrals.
- `thermal`: theta-method time integrator shared by 1D and 2D.
- `mechanics`: constitutive models behind one interface
  (`stress(strain_increment, history) -> stress, tangent, new_history`), global
  Newton-Raphson with load stepping.
- `coupling`: thermal strain from a temperature field on the same mesh.
- `io`: VTU writer, CSV writer, time-series reader (CSV `t_s,value` with linear interpolation).
- Cases are described by small JSON files in `cases/`, never hard-coded in `main`.

Solvers: `Eigen::SimplicialLDLT` for symmetric systems, `Eigen::SparseLU` otherwise.

## 4. Phases

### Phase 1: 1D transient heat and diffusion

Build:
- 1D linear finite elements on a graded or uniform grid; consistent and lumped capacity.
- Theta-method: theta = 1 (backward Euler) and theta = 0.5 (Crank-Nicolson).
- Boundary conditions: Dirichlet, prescribed flux, convective `-k dT/dn = h (T - T_inf(t))`
  with `T_inf(t)` from a time-series file.
- Diffusion reuses the same code with `rho*c = 1`, `k = D`.
- Arrhenius helper `D(T) = D0 exp(-Q / (R T))`.

Validation (all in REFERENCE_VALUES.md, section Phase 1):
- V1.1 Carburizing, semi-infinite: `C = Cs - (Cs - C0) erf(x / (2 sqrt(D t)))`. Domain length
  at least 10 sqrt(D t_end). Time to reach 0.8 wt% at 0.5 mm.
- V1.2 Thin-film source on a half-space (zero flux at x = 0):
  `C = M / sqrt(pi D t) * exp(-x^2 / (4 D t))`. Start from the analytic profile at a small
  t0 > 0 (not a delta). Check total dose conserved to 1e-10 relative and the peak decaying
  as t^(-1/2).
- V1.3 Arrhenius and equivalent-time exercises as pure unit tests.
- V1.4 Steel plate, 5 mm, convective both faces (parameters in `references/data/steel_plate_case.json`).
  Until the real ttab file arrives, drive it with a synthetic daily sol-air signal and compare
  with the lumped ODE `rho c L dT/dt = h_e (T_sa - T) + h_i (T_in - T)` solved with RK4 at a
  small step. Use backward Euler here: dt = 300 s is about 200 times the through-thickness
  diffusion time, where Crank-Nicolson rings.
  When `references/data/ttab.csv` and `references/data/ansys_inner_face.csv` exist, a second
  test compares against the ANSYS curve and reports the RMS difference; skip it otherwise.
- V1.5 Convergence: observed temporal order about 1 (BE) and 2 (CN), spatial order about 2,
  each within +/- 0.15, measured on at least three refinements.

Done when: all phase1 tests pass; figures in `docs/figures/`: erf profile vs numeric,
thin-film decay, steel plate inner face vs lumped model, convergence plots.

### Phase 2: 2D linear elasticity

Build:
- Gmsh reader with physical groups for boundaries and regions; `.geo` files in `cases/`.
- T3 and T6 elements; plane stress and plane strain elasticity matrices. T6 must be
  isoparametric with mid-side nodes on the hole boundary (Gmsh second-order meshes do
  this); straight-sided T6 on a curved boundary caps the convergence order.
- Traction loading by edge integration; Dirichlet by elimination.
- Stress at Gauss points plus nodal averaging; von Mises; VTU output.

Validation (section Phase 2):
- V2.1 Patch test: a linear displacement field is reproduced to 1e-10 relative on a
  distorted mesh, for T3 and T6.
- V2.2 Uniform bar in tension: end displacement `sigma L / E`.
- V2.3 Kirsch plate with a hole. Quarter model using symmetry, hole radius a, outer square of
  side 10a. Apply the EXACT Kirsch tractions on the outer boundary so the finite model
  represents the infinite plate; then the FE solution converges to the analytic field.
  Remote stress 140 MPa along y. Compare against `references/data/kirsch_sigma140MPa.csv`
  (theta measured from the load axis y, as in the source). Report the stress
  concentration factor (exact value 3).
- V2.4 Convergence: displacement L2 error order about 2 (T3) and 3 (T6); energy-norm order
  about 1 (T3) and 2 (T6).

Done when: phase2 tests pass; T6 finest mesh gives Kt within 1% of 3 and every table point at
r = 1.5a, 2a and 4a within 2 MPa; VTU of von Mises around the hole; convergence plot.

### Phase 3: Elastoplasticity

Build:
- Small-strain J2 plasticity with isotropic hardening, radial return mapping, consistent
  (algorithmic) tangent. Two hardening laws behind one interface:
  - Tabulated: piecewise-linear yield stress vs plastic strain built from
    `tensile_true.csv` (eps_p = eps_true - sigma / E). This is what the Abaqus run in the
    source uses: its strain at 457 MPa equals a data point exactly.
  - Ludwik: `sigma_eq = sigma_y + K eps_p^m`, fitted below.
  At 457 MPa a Ludwik fit misses the data by 1.3% (our fit) to 3.0% (source fit) in stress,
  which becomes 11% to 27% in plastic strain because the curve is flat there. So P3-b uses the tabulated law; Ludwik is
  checked against the data within its fit residual (P3-c).
- Global Newton-Raphson with load stepping; history variables at Gauss points.
- Plane strain first. Plane stress (local iteration enforcing sigma_zz = 0) only after plane
  strain passes.
- A material-point driver for uniaxial stress (iterate lateral strains until lateral
  stresses vanish): this is the tensile test.
- Parameter identification: least-squares fit of (sigma_y, K, m) to
  `references/data/tensile_true.csv` with an uncertainty estimate (covariance or bootstrap).
  Report next to the source values (sigma_y 310.86 MPa, K 559.27 MPa, m 0.4077).

Validation (section Phase 3):
- V3.1 Material point with tabulated hardening reproduces the true stress-strain curve; at
  sigma = 457.24 MPa with E = 207.4 GPa and nu = 0.3, strains eps_yy, eps_xx = eps_zz match
  the reference to 0.5%.
- V3.2 Newton converges quadratically (observed rate at least 1.8 in the last iterations of
  a step), which proves the tangent is consistent.
- V3.3 Uniform plane-strain bar in FE equals the material-point driver with eps_zz = 0
  to 1e-6 relative at every load step.
- V3.4 Hole plate at remote 140 MPa: the elastic peak (420 MPa) exceeds yield, so a plastic
  zone forms. Report peak strain and stress at the hole edge and the plastic-zone extent;
  compare with Neuber's rule `sigma * eps = (Kt S)^2 / E` solved on the same hardening curve
  as the FE run. Use plane stress here once it exists (a thin plate, like the source's
  Abaqus model); until then run plane strain and say so in the output.
  Neuber is an estimate: report the difference, never force agreement.
- Optional: unload to zero and report residual stress at the hole edge.

Done when: phase3 tests pass; figures: fitted curve with data and uncertainty band, Newton
residual history, plastic zone map.

### Phase 4: Thermo-mechanical coupling and the beam-heated block

Build:
- 2D transient heat on T3/T6 with volumetric source q(x, y, t) and convective boundaries,
  reusing the phase 1 integrator and phase 2 assembly.
- Thermal strain `alpha (T - T_ref)` in the mechanical solve, for plane stress and plane
  strain. Derive both thermal load vectors in `docs/learning/phase4.md` before coding.

Validation (section Phase 4):
- V4.1 Bar in plane stress, held in x at both ends and free laterally, uniform rise dT:
  `sigma_xx = -E alpha dT`, `sigma_yy = 0`.
- V4.2 Unconstrained body, uniform dT: zero stress to 1e-10 relative, displacement
  `alpha dT x`.
- V4.3 Thin solid disk (plane stress), radius b, `T(r) = dT (1 - r^2 / b^2)`:
  `sigma_r = -(E alpha dT / 4)(1 - r^2/b^2)`, `sigma_theta = -(E alpha dT / 4)(1 - 3 r^2/b^2)`
  (Timoshenko and Goodier). Within 1% at r = 0, b/2, b on the finest mesh, plus convergence order.
- V4.4 2D heat: an initial Gaussian of width s0 in a large domain spreads as a Gaussian of
  variance s0^2 + 2 alpha_th t with conserved energy. Check against the closed form.

Demo D4 (not a validation, label it as such): a rectangular copper block, then isotropic
graphite, heated by a Gaussian "beam spot" of rms size s_beam depositing energy E_pulse
over a short pulse, adiabatic during the pulse, convective cooling after. Output
temperature and von Mises maps and a yield check; sweep s_beam. Take every material
property from a cited source and record it in `cases/beam_block/README.md`. State that real
energy maps come from FLUKA or Geant4 and that the Gaussian is a stand-in.

Done when: phase4 tests pass; demo figures and a one-paragraph interpretation.

### Phase 5: Report and polish

- README: one-paragraph summary, results table (benchmark, quantity, reference, tmfem,
  error, test name), key figures, build instructions, limitations, credits.
- GitHub Actions: build and run ctest on ubuntu-latest for every push.
- Technical report in `report/` (LaTeX, 8 to 10 pages): context, methods, validation,
  demo, limitations, references.
- Credits: benchmark cases adapted from M. R. Kabiri, ENSAM Meknès lab materials (with
  permission); analytical solutions from cited textbooks.

### Stretch goals (only after phase 5 is presentable)

- S1 Center-cracked strip, K_I by displacement correlation and by domain J-integral, against
  `K = sigma sqrt(pi a) / sqrt(cos(pi a / (2 b)))` (Feddersen). Show the J-integral beats
  single-point extraction.
- S2 OpenMP assembly with a timing table.
- S3 Small ML surrogate (Python) predicting peak temperature and von Mises from beam
  parameters, trained on solver sweeps, with held-out error.

## 5. Conventions

- Voigt order [xx, yy, xy] in 2D with engineering shear strain; zz tracked separately for
  plane strain plasticity.
- Quadrature: exact for the element's stiffness (T3: 1 point, T6: 3 points; 6 points for
  nonlinear residuals).
- Sign convention: tension positive.
- Output files store SI units except stresses, which are written in MPa with the unit in
  the column header.

## 6. Definition of done (whole project)

A stranger can clone the repo, run three commands, and see every validation test pass;
every number in the README traces to a test; the report states all limitations.
