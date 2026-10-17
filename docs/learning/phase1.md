# Phase 1: 1D transient heat conduction and diffusion

This note derives everything the phase 1 code does. Each section ends with the place in the
code that implements it, so the derivation and the implementation can be checked against
each other.

Contents:
- 0: what is computed
- 1-10: derivation
- 11: algorithm
- 12: validation cases
- 13: what went wrong
- 14: limitations

## 0. What this phase computes

tmfem phase 1 computes how temperature (or a diffusing concentration) evolves in time through
the thickness of a slab, given its material, its initial state, and what happens at its two faces
(fixed temperature, fixed heat flux, or convection to a fluid whose temperature can follow a
measured time series).

It does this with linear finite elements in space and the theta method in time (backward Euler
or Crank-Nicolson).

Every capability is checked against a closed-form answer:
- the erf and thin-film diffusion solutions
- Arrhenius exercises
- a lumped model of a thin steel plate
- the measured orders of convergence

## 1. Strong form

Heat conduction in a slab of thickness L (only x varies):

$$\rho c\,\frac{\partial T}{\partial t} = \frac{\partial}{\partial x}\left(k\,\frac{\partial T}{\partial x}\right) + q(x,t), \qquad 0 < x < L,\qquad T(x,0) = T_0(x).$$

- `rho c` is the volumetric heat capacity (J/(m3 K)), `k` the conductivity (W/(m K)), `q` a volumetric source (W/m3).
- Each end has an outward normal n: n = -1 at x = 0 and n = +1 at x = L, so `dT/dn = n dT/dx`.

Three kinds of boundary condition:

| Type | Condition | Meaning |
|---|---|---|
| Dirichlet | `T = g(t)` | temperature imposed |
| Flux | `k dT/dn = q_in(t)` | heat flux `q_in` (W/m2) *entering* the body |
| Convective (Robin) | `-k dT/dn = h (T - T_inf(t))` | Newton cooling to a fluid at `T_inf` with film coefficient h |

**Diffusion is the same equation.** Fick's second law `dC/dt = d/dx(D dC/dx)` is the heat equation
with `rho c = 1`, `k = D`, `T = C`, `q = 0`. The code has a single solver; the diffusion tests just
pass `rho_c = 1`, `k = D`.

## 2. Weak form

Take a test function w(x) that vanishes wherever T is imposed (Dirichlet ends). Multiply the PDE by
w and integrate over the slab:

$$\int_0^L w\,\rho c\,\dot T\,dx = \int_0^L w\,(k T')'\,dx + \int_0^L w\,q\,dx.$$

Integrate the conduction term by parts:

$$\int_0^L w\,(kT')'\,dx = \big[w\,k\,T'\big]_0^L - \int_0^L k\,w'\,T'\,dx.$$

The boundary bracket is `w(L) k T'(L) - w(0) k T'(0)`. At either end that equals `w k dT/dn`
(at x = 0, `dT/dn = -T'`). Insert the boundary conditions:

- Dirichlet end: w = 0, so the term disappears (the unknown reaction flux never appears).
- Flux end: `w k dT/dn = w q_in`.
- Convective end: `w k dT/dn = -w h (T - T_inf)`.

The term `-w h T` contains the unknown, so it moves to the left-hand side. The weak form is: find T
(satisfying the Dirichlet conditions) such that for all admissible w

$$\int_0^L \rho c\,w\,\dot T\,dx + \int_0^L k\,w'\,T'\,dx + \sum_{\text{Robin}} h\,w\,T
= \int_0^L w\,q\,dx + \sum_{\text{flux}} w\,q_{in} + \sum_{\text{Robin}} h\,w\,T_\infty(t).$$

Only first derivatives appear, so continuous piecewise-linear functions are admissible.

## 3. Galerkin discretisation

Write `T_h(x,t) = sum_j N_j(x) T_j(t)` with nodal "hat" functions N_j, and take w = N_i for every
free node i:

$$M\,\dot{\mathbf T} + (K + H)\,\mathbf T = \mathbf F(t)$$

$$M_{ij} = \int \rho c\,N_i N_j\,dx,\quad K_{ij} = \int k\,N_i' N_j'\,dx,\quad
H_{ii} = h \text{ at a Robin node},\quad F_i = \int N_i q\,dx + q_{in}\big|_{\text{flux node}} + h T_\infty(t)\big|_{\text{Robin node}}.$$

In 1D the boundary is a point and the shape function there equals 1, so the Robin "boundary
integral" is the single number h added to one diagonal entry. In 2D (phase 4) the same term becomes
`int h N_i N_j ds` along boundary edges.

Code: `src/thermal/heat1d.cpp` (`build_operators`).

## 4. Linear element matrices

Element [x1, x2], length `h_e = x2 - x1`, reference coordinate xi in [-1, 1]:

$$N_1 = \tfrac{1-\xi}{2},\quad N_2 = \tfrac{1+\xi}{2},\quad x = N_1 x_1 + N_2 x_2,\quad
dx = \tfrac{h_e}{2}\,d\xi,\quad \frac{dN}{dx} = \frac{dN}{d\xi}\frac{d\xi}{dx} = \frac{1}{h_e}(-1,\ 1).$$

**Stiffness** (k constant per element):

$$K^e_{ij} = \int_{-1}^{1} k\,\frac{dN_i}{dx}\frac{dN_j}{dx}\,\frac{h_e}{2}\,d\xi
= k\,\frac{(\pm1)(\pm1)}{h_e^2}\,h_e \;\Rightarrow\; K^e = \frac{k}{h_e}\begin{bmatrix}1 & -1\\ -1 & 1\end{bmatrix}.$$

**Consistent capacity:**

$$\int_{-1}^1 N_1^2\,\tfrac{h_e}{2}d\xi = \tfrac{h_e}{2}\int_{-1}^1 \tfrac{(1-\xi)^2}{4}d\xi = \tfrac{h_e}{2}\cdot\tfrac{2}{3} = \tfrac{h_e}{3},\qquad
\int_{-1}^1 N_1N_2\,\tfrac{h_e}{2}d\xi = \tfrac{h_e}{2}\int_{-1}^1 \tfrac{1-\xi^2}{4}d\xi = \tfrac{h_e}{6}$$

$$M^e_{\text{cons}} = \frac{\rho c\,h_e}{6}\begin{bmatrix}2 & 1\\ 1 & 2\end{bmatrix}.$$

**Lumped capacity** (row-sum lumping): each row of `M_cons` sums to `rho c h_e / 2`, so

$$M^e_{\text{lump}} = \frac{\rho c\,h_e}{2}\begin{bmatrix}1 & 0\\ 0 & 1\end{bmatrix}.$$

Both versions have total `sum_ij M_ij = rho c h_e`, the element's heat capacity per unit area, so
both conserve energy. K has zero row sums, because a uniform temperature carries no heat flux.

**Uniform source:** `F^e_i = int q N_i dx = q h_e / 2` for each node.

The code computes K and M by 2-point Gauss quadrature (exact for the quadratic integrands
here) and the source vector by 3-point Gauss, and the unit test `element matrices match closed form` checks them against the formulas above.

Code: `src/fem/line2.cpp`, `src/fem/quadrature.cpp`.

## 5. Theta method in time

Let `A = K + H`. Evaluate the semi-discrete equation as a weighted average of the two ends of the
step [t_n, t_{n+1}], with `dt = t_{n+1} - t_n`:

$$M\,\frac{\mathbf T^{n+1}-\mathbf T^n}{\Delta t} + A\,[\theta\mathbf T^{n+1} + (1-\theta)\mathbf T^n] = \theta\mathbf F^{n+1} + (1-\theta)\mathbf F^n$$

$$\Rightarrow\quad (M + \theta\Delta t\,A)\,\mathbf T^{n+1} = (M - (1-\theta)\Delta t\,A)\,\mathbf T^n + \Delta t\,[\theta\mathbf F^{n+1} + (1-\theta)\mathbf F^n].$$

- theta = 1: backward Euler (BE). First order: the local error is `dt^2/2 * T''`.
- theta = 1/2: Crank-Nicolson (CN). Second order: it is the trapezoidal rule, symmetric about the
  midpoint, so the `dt^2` term of the Taylor expansion cancels.

**Stability.** Decompose into generalised eigenmodes `A v = lambda M v` (lambda >= 0). Each mode
evolves by the amplification factor

$$G(\lambda\Delta t) = \frac{1 - (1-\theta)\lambda\Delta t}{1 + \theta\lambda\Delta t},$$

while the exact factor is `exp(-lambda dt)`. `|G| <= 1` for every `lambda dt >= 0` exactly when
theta >= 1/2 (unconditional stability).

**Why CN rings and BE does not.** As `lambda dt -> infinity`:
- `G_CN -> -1`: a stiff mode flips sign every step and hardly decays. That is an oscillation, not a
  diffusion.
- `G_BE = 1/(1 + lambda dt) -> 0`: stiff modes are damped at once (L-stability).

In the steel plate the lowest through-thickness conduction mode has
`lambda ~ pi^2 k/(rho c L^2) = pi^2 * 60/(7850*434*0.005^2) ~ 7 s^-1`. With dt = 300 s,
`lambda dt ~ 2000` and `G_CN ~ -0.999`. Any sharp feature in the forcing would ring for hundreds of
steps, which is why P1-j uses BE.

**Monotonicity.** With a step initial condition, as in carburizing (surface jumps to Cs at t = 0):
- With consistent mass, BE can undershoot below C0 near the front unless `dt >= rho c h^2/(6k)`.
  The off-diagonal of `M + dt K` must be non-positive for the discrete maximum principle:
  `rho c h/6 - dt k/h <= 0`.
- Lumped mass has a zero off-diagonal in M, so it is monotone for any dt.

P1-b and P1-e therefore use lumped mass. That choice is made for this reason, before seeing any
result.

Code: `src/thermal/theta_integrator.cpp`.

## 6. Dirichlet conditions by elimination

Partition the nodes into free (f) and constrained (c). With `S = M + theta dt A` and the
right-hand side b from section 5:

$$\begin{bmatrix}S_{ff} & S_{fc}\\ S_{cf} & S_{cc}\end{bmatrix}\begin{bmatrix}\mathbf T_f\\ \mathbf g^{n+1}\end{bmatrix} = \begin{bmatrix}\mathbf b_f\\ \cdot\end{bmatrix}
\quad\Rightarrow\quad S_{ff}\,\mathbf T_f^{n+1} = \mathbf b_f - S_{fc}\,\mathbf g^{n+1}.$$

The constrained rows are dropped, because they hold the unknown reaction fluxes. `S_ff` is
symmetric positive definite (M SPD, K positive semi-definite, H >= 0), so it is factorised once
per dt with `Eigen::SimplicialLDLT` and reused every step.

Code: `src/fem/dirichlet.cpp`.

## 7. Discrete conservation (thin film, P1-h)

With zero-flux ends there is no Dirichlet, no Robin and no source. Then:
- H = 0 and F = 0.
- K has zero row sums, so `1^T K = 0` (K is symmetric).

Multiply the theta step by `1^T`:

$$\mathbf 1^T M\,\mathbf T^{n+1} + \theta\Delta t\,\underbrace{\mathbf 1^T K}_{0}\,\mathbf T^{n+1} = \mathbf 1^T M\,\mathbf T^n - (1-\theta)\Delta t\,\underbrace{\mathbf 1^T K}_{0}\,\mathbf T^n.$$

So the discrete dose `1^T M T` is constant for any theta, any dt and any mesh. It changes only by
the linear solver's residual, which is why it can be checked to 1e-10.

## 8. Analytic solutions used by the tests

**Carburizing (semi-infinite, constant surface concentration).**
- The equation is invariant under `x -> lambda x, t -> lambda^2 t`, so C depends only on
  `eta = x / (2 sqrt(D t))`.
- The ODE in eta gives `C = Cs - (Cs - C0) erf(eta)`.
- The time to reach concentration C at depth x is `t = x^2/(4 D z^2)` with `erf(z) = (Cs - C)/(Cs - C0)`.
- A finite domain behaves like the half-space if its length is at least `10 sqrt(D t)`:
  `erfc(5) ~ 1.5e-12`, so the far end never feels the front.

**Dose from a constant source.** `M = int_0^inf (C - C0) dx`, here with C0 = 0. Since
`int_0^inf erfc(eta) d eta = 1/sqrt(pi)`, this gives `M = 2 Cs sqrt(D t / pi)`.

**Thin film on a zero-flux half-space.**
- For an infinite line, the Gaussian `M/(2 sqrt(pi D t)) exp(-x^2/(4Dt))` spreads a dose M.
- A zero-flux wall at x = 0 acts as a mirror: an image source adds the same Gaussian. This
  doubles the concentration:

$$C(x,t) = \frac{M}{\sqrt{\pi D t}}\,e^{-x^2/(4Dt)},\qquad \int_0^\infty C\,dx = M,\qquad C(0,t) = \frac{M}{\sqrt{\pi D t}} \propto t^{-1/2}.$$

- The test starts from this profile at t0 > 0 (a delta function cannot be represented on a mesh).

**Arrhenius.** `D(T) = D0 exp(-Q/(R T))` with T in kelvin.
- From two measurements, `ln(D2/D1) = (Q/R)(1/T1 - 1/T2)`, so `Q = R ln(D2/D1)/(1/T1 - 1/T2)` and
  `D0 = D1 exp(Q/(R T1))`.
- Profiles depend only on `D t`, so equal `D t` means an equal profile: `t2 = t1 D1/D2`
  (equivalent time) and `x2 = x1 sqrt(D2/D1)` at equal time.

## 9. Steel plate: lumped model and through-thickness difference

Integrate the heat equation over the thickness, using the two convective conditions:

$$\rho c\,L\,\frac{d\bar T}{dt} = h_e\,(T_{sa} - T(0)) + h_i\,(T_{in} - T(L)).$$

- The Biot number `Bi = h_e L/k = 0.00183` compares the conduction resistance `L/k` with the film
  resistance `1/h_e`.
- When `Bi << 1`, `T(0) ~ T(L) ~ T_mean`, giving the lumped ODE `rho c L dT/dt = h_e (T_sa - T) + h_i (T_in - T)`.
- Its time constant is `tau = rho c L/(h_e + h_i) = 549.5 s`.

Quasi-steady (forcing slow compared with the 1.4 s diffusion time), the same flux crosses the
three resistances in series, `R_tot = 1/h_e + L/k + 1/h_i`. So:

$$T(0) - T(L) = \frac{L/k}{R_{tot}}\,(T_{sa} - T_{in}) = 5.32\times10^{-4}\,(T_{sa} - T_{in}).$$

For the synthetic signal `T_sa = 50 - 24 cos(2 pi t / day)` °C (26 to 74 °C), the maximum
difference is `5.32e-4 * 48 K = 0.0255 K`. That is the P1-k target. The older "below 0.01 K" holds
only while `|T_sa - T_in| < 19 K`.

The formula is quasi-steady: it ignores heat stored in the plate. Storage makes the outer flux
differ from the inner one by `rho c L dT/dt`, and `T(0) - T(L) ~ (L/k) (q_e + q_i)/2`. The relevant
parameter is the plate's time constant relative to the forcing period,
`omega tau = 2 pi * 549.5 s / 86400 s = 0.04`. At the daily peak, `dT/dt` is near zero, so the
storage correction is small there. A first-order estimate gives a few tenths of a percent, inside the
5% tolerance of P1-k.

## 10. Measuring the order of convergence

If `e(h) ~ C h^p`, then on two levels with ratio r = h_k/h_{k+1}:

$$p = \frac{\log(e_k/e_{k+1})}{\log r}.$$

The tests use four or more levels, assert the least-squares slope of log e against log h (or
log dt), and print every pairwise p.

Expected with linear elements:
- `L2` error `O(h^2)`.
- Energy (H1-seminorm) error `O(h)`.
- BE `O(dt)`, CN `O(dt^2)`.

**Isolating the time error.** On a uniform mesh with zero Dirichlet ends:
- The nodal vector `v_i = sin(pi x_i/L)` is an exact eigenvector of both K and M (consistent or
  lumped), with discrete eigenvalue `lambda_h`. This holds because the stencils are symmetric
  three-point stencils and `sin(a - b) + sin(a + b) = 2 sin(a) cos(b)`.
- So the semi-discrete solution is exactly `exp(-lambda_h t) v`, and after n steps the theta method
  gives `G(lambda_h dt)^n v`.
- Comparing the two measures the time-integration error alone, with zero spatial error.
- The spatial order is measured separately against the exact PDE solution
  `exp(-pi^2 alpha t/L^2) sin(pi x/L)`, with dt small enough that the time error is under 1% of
  the spatial error. The test checks this ratio rather than assuming it.

**The eigenvector step, written out.** Uniform mesh, spacing h, interior node i, `v_i = sin(i phi)`
with `phi = pi h / L`. Using `v_{i-1} + v_{i+1} = 2 cos(phi) v_i`:

$$(K v)_i = \frac{k}{h}(-v_{i-1} + 2v_i - v_{i+1}) = \frac{2k}{h}(1-\cos\phi)\,v_i,\qquad
(M_{\text{cons}} v)_i = \frac{\rho c\,h}{6}(v_{i-1} + 4v_i + v_{i+1}) = \frac{\rho c\,h}{3}(2+\cos\phi)\,v_i,$$

$$(M_{\text{lump}} v)_i = \rho c\,h\,v_i.$$

The boundary nodes are Dirichlet zeros, and `sin(0) = sin(pi) = 0` closes the pattern. Hence:

$$\lambda_h^{\text{cons}} = \frac{6k}{\rho c\,h^2}\,\frac{1-\cos\phi}{2+\cos\phi} = \lambda\left(1 + \frac{\phi^2}{12} + O(\phi^4)\right),\qquad
\lambda_h^{\text{lump}} = \frac{2k}{\rho c\,h^2}(1-\cos\phi) = \lambda\left(1 - \frac{\phi^2}{12} + O(\phi^4)\right),$$

with `lambda = pi^2 k / (rho c L^2)`.

Both errors are `O(h^2)`, which is the spatial order P1-m measures. They have opposite signs:
- Consistent mass makes the mode decay slightly too fast.
- Lumped mass makes it decay slightly too slowly.

## 11. The algorithm, step by step

1. **Read the case.** `src/io/case_heat1d.cpp:load_heat1d_case` (line 46).
   - Parses the JSON.
   - Converts °C to K at the input boundary.
   - Wraps the outer air temperature series in a `TimeSeries` (`src/io/timeseries.cpp:TimeSeries::operator()`,
     linear interpolation, throws outside the data).
2. **Build the mesh.** `src/mesh/interval_mesh.cpp:graded_interval` (line 13).
   - Node `i` sits at `x_i = L (r^i - 1)/(r^n - 1)`; `r = 1` gives the uniform mesh.
3. **Element matrices.** By 2-point Gauss (`src/fem/quadrature.cpp:gauss_legendre`):
   - `src/fem/line2.cpp:element_stiffness` (line 16)
   - `src/fem/line2.cpp:element_capacity` (line 28; row-sum lumping when `MassType::Lumped`)
   - `src/fem/line2.cpp:element_source` (line 43, 3-point Gauss)
4. **Assemble M and A = K + H, and define F(t) and the Dirichlet values g(t).**
   `src/thermal/heat1d.cpp:build_operators` (line 40):
   - A Robin end adds h to `A(node, node)` and `h T_inf(t)` to `F(node)`.
   - A flux end adds `q_in(t)` to `F(node)`.
   - A Dirichlet end is recorded as a constrained dof.
5. **Partition the dofs** into free and constrained: `src/fem/dirichlet.cpp:DirichletPartition`
   (`block_ff`, `block_fc` at lines 44 and 48).
6. **Factorise once per dt.** `src/thermal/theta_integrator.cpp:ThetaIntegrator::factorize` (line 23):
   - Forms `S = M + theta dt A` and `R = M - (1 - theta) dt A`.
   - LDLT-factorises `S_ff`.
   - `set_dt` (line 18) refactorises when dt changes (the thin film uses dt = 0.01 t).
7. **One step.** `src/thermal/theta_integrator.cpp:ThetaIntegrator::step` (line 34):
   - `b = R T^n + dt [theta F^{n+1} + (1 - theta) F^n]`
   - `rhs_f = b_f - S_fc g^{n+1}`
   - Solve `S_ff T_f = rhs_f`, then put free and constrained values back together.
8. **March.** `src/thermal/theta_integrator.cpp:march` (line 45). Time is `t_n = t0 + n dt`,
   computed fresh each step rather than accumulated, so no rounding drift builds up over thousands of steps.
9. **Post-process:**
   - probes by `src/mesh/interval_mesh.cpp:interpolate` (line 41)
   - fluxes by `src/thermal/heat1d.cpp:element_flux` (line 95)
   - CSV by `src/io/csv_writer.cpp:write_csv`
   - CLI run by `src/app/main.cpp:run_heat1d` and `src/io/case_heat1d.cpp:run_heat1d_case` (line 73)

The steady solve `src/thermal/theta_integrator.cpp:solve_steady` (line 57) is step 7 with `M = 0`
and theta = 1. P1-a and the unit tests use it.

## 12. Validation cases: what each proves and what bug it would catch

| ID | What it proves | A bug it would catch |
|---|---|---|
| P1-a | Steady assembly + Dirichlet elimination reproduce a linear profile exactly; flux sign convention `J = -D dC/dx` | wrong sign of the flux; Dirichlet values applied to the wrong node; mm not converted to m (the source's own 1000x error) |
| P1-b | Transient diffusion with a Dirichlet surface converges to the erf solution; time interpolation of a threshold crossing | wrong theta weighting; a far boundary that leaks (domain too short); consistent-mass undershoot near the front |
| P1-c, d, g | Arrhenius algebra, `D t` scaling, °C to K | using °C in the exponent (a 20x error in D); `ln` vs `log10`; cm²/s vs m²/s (the source's D0 slip) |
| P1-e | Predeposition dose formula, and the FE dose `1^T M C` converging to it | missing factor 2 in `2 Cs sqrt(D t/pi)`; dose computed with the wrong mass matrix |
| P1-f | Thin-film surface formula with the correct time | the source's error (t = 36 min instead of 6 min) |
| P1-h | Exact discrete conservation with zero-flux ends; correct self-similar decay | K without zero row sums; a non-symmetric assembly; a natural BC accidentally turned into Dirichlet; dt mishandled when it changes |
| P1-i | Biot number and time constant from the data file | typo in parameters; the case file drifting from the reference data (the test checks every value) |
| P1-j | Robin BCs with a time-varying `T_inf`, end to end through the JSON case and the time series | `h T_inf` added with the wrong sign or at the wrong node; `T_inf` evaluated at `t_n` instead of `t_{n+1}` in BE; °C/K mix-up (would give 273 K offsets) |
| P1-k | The FE resolves the tiny through-thickness gradient correctly | conduction term scaled wrongly (e.g. `k` vs `k/L`), which would change the 0.0255 K directly |
| P1-l | (pending data) comparison with ANSYS | report only |
| P1-m | BE order 1, CN order 2, space order 2 for both mass matrices | CN implemented as BE (order 1); F evaluated at the wrong time level (order drops); a wrong element matrix (order drops or no convergence) |

The P1-b, P1-e, P1-h and P1-j tests also rerun at half the step size and half the mesh size, and assert
that the answer barely moves (less than 1/5 of the tolerance). This shows that the resolution was not
picked to land on the reference value.

## 13. What went wrong while building this phase, and how it was caught

Honest record: all P1 reference tests passed on the first run. No solver bug surfaced, so I cannot
claim the tests "caught" one; the convergence orders and the halved-resolution reruns are what give
confidence beyond a single passing number. What did go wrong:

1. **An infeasible reference value (P1-k).**
   - "Below 0.01 K" was found before any code, by a resistance estimate:
     `(L/k)/R_tot * 48 K = 0.0255 K`.
   - It was raised as a question and replaced only after explicit review. The test now
     checks the analytic bound. The FE gives 0.02558 K vs 0.02553 K.
2. **Two unit-test tolerances below what floating point allows.**
   - `erfinv` at z = 3.5: the inverse amplifies rounding by `exp(z^2) ~ 2e5`.
   - The steady Robin system: cond(A) about 1.7e5.
   - Both failed at 1e-13. The fix was in the tests, which now check what is attainable. The code was
     right. Lesson: a tolerance must be justified by an error estimate, not chosen to look strict.
3. **The plot script read a comment line as the CSV header.** numpy `genfromtxt(names=True)` takes
   the first line even if it is a comment. It crashed immediately and was fixed by stripping comments first.
4. **Environment.** A user-installed numpy 2 broke the system matplotlib. The scripts re-exec with
   `python3 -s`.
5. **The git repository was rooted one directory too high** (`~/projects`). As a result,
   `.github/workflows/ci.yml` was not at the repository root and CI never ran during phase 0.
   It was re-rooted with `git subtree split`.
6. **A visible artefact that is not a bug.** In the steel-plate figure, FE − lumped wiggles with a
   one-hour period.
   - BE lags the forcing by about dt/2, and the hourly-sampled T_sa changes slope every hour.
   - Halving dt changes the inner face by only 0.0068 K, which is consistent with this explanation.

## 14. Limitations and the assumptions behind them

- **1D only.** Lateral faces are adiabatic. For the steel plate this is the source's own model
  (a 10 cm plate, 5 mm thick).
- **Constant properties.**
  - k, rho c and D are constant within an element. The solver does not let them depend on temperature.
  - The Arrhenius `D(T)` is used only in the exercises, at one temperature at a time.
  - Non-isothermal diffusion would need nonlinear iteration.
- **Linear elements only:** second order in space, at best.
- **Boundary conditions.**
  - Convection is linear in T; there is no radiation (that would need a nonlinear Robin term).
  - Only the two ends of the slab carry conditions.
- **Time stepping.**
  - Fixed dt per factorisation; a dt change costs a refactorisation.
  - There is no adaptive error control.
  - CN is reserved for smooth data, because it rings on stiff modes.
- **Steel plate input.**
  - The input signal is synthetic until the real ttab arrives (P1-l is skipped).
  - The agreement shown is with a lumped model, not with ANSYS.
  - The source's claim that rising daily peaks come from "accumulated heat" cannot hold with a
    9-minute time constant.
- **The P1-m eigenmode trick needs a uniform mesh** (the sine is an exact eigenvector only there).
  The spatial order on graded meshes is not measured separately; P1-h and P1-b show that graded
  meshes converge.
