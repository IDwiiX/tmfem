# Phase 4: 2D transient heat and one-way thermo-mechanical coupling

This note derives what the phase 4 code does. Each section ends with the place in the code that
implements it.

Contents:
- 0: what is computed
- 1-5: 2D transient heat
- 6-9: thermal strain and the thermal load vector
- 10-11: analytic solutions and the demo model
- 12-16: algorithm, validation cases, what went wrong, limitations and interview questions

Notation: alpha is the coefficient of thermal expansion (1/K). The thermal diffusivity is written
`a = k/(rho c)`, which SPEC calls alpha_th. Phase 3 used alpha for the hardening variable; nothing here
refers to it. dT = T - T_ref. m = [1, 1, 0]^T in Voigt order [xx, yy, xy].

## 0. What this phase computes

A temperature field T(x, y, t) in a 2D body heated by a volumetric source and cooled by convection,
followed by the stress that the thermal expansion causes. The coupling is one-way: temperature
changes stress, but stress does not change temperature (the thermoelastic coupling term
T0 * beta * d(tr eps)/dt is negligible for the slow processes considered here). The mechanics is
quasi-static and linear elastic.

## 1. Strong form of 2D heat conduction

In Omega with thickness t_z (the plate) or unit depth (a long body):

| Equation | Meaning |
|---|---|
| `rho c dT/dt = div(k grad T) + q(x, t) - s(T)` | energy balance per unit volume |
| `T = g(t)` on Gamma_D | imposed temperature |
| `-k dT/dn = h (T - T_inf(t))` on Gamma_R | convection on an edge (n outward) |
| `-k dT/dn = -q_in(t)` on Gamma_N | imposed flux into the body |
| `T(x, 0) = T0(x)` | initial condition |

The face sink `s(T)` exists only for a thin plate cooled on both faces. Integrating the 3D heat
equation through the thickness, with `-k dT/dz = h_f (T - T_inf)` on each face and T nearly uniform in z
(through-thickness Biot `h_f t_z/(2k)` much less than 1: for example 1.3e-4 for copper and 4.8e-4 for graphite
with h_f = 10 W/(m^2 K) and t_z = 10 mm), gives

$$s(T) = \frac{2h_f}{t_z}\,(T - T_\infty).$$

Each face takes `h_f (T - T_inf)` per unit area; two faces, divided by the thickness, give a loss per unit
volume. With no faces cooled, s = 0. The sink is part of the solver and is checked by P4-f (ii) and (iii).
The demo does not use it: its block is long along z (section 11), so it has no cooled faces in the model.

## 2. Weak form

Multiply by a test function w (w = 0 on Gamma_D), integrate over Omega and use
`div(w k grad T) = w div(k grad T) + k grad w . grad T` with the divergence theorem:

$$\int_\Omega w\,\rho c\,\dot T\,dA + \int_\Omega k\,\nabla w\cdot\nabla T\,dA + \int_{\Gamma_R} h\,w\,T\,ds + \int_\Omega \tfrac{2h_f}{t_z}\,w\,T\,dA
= \int_\Omega w\,q\,dA + \int_{\Gamma_R} h\,w\,T_\infty\,ds + \int_{\Gamma_N} w\,q_{in}\,ds + \int_\Omega \tfrac{2h_f}{t_z}\,w\,T_\infty\,dA.$$

- The boundary term from integration by parts is `int w k dT/dn ds`. On Gamma_R it equals `-int w h (T - T_inf) ds`.
  The part with T moves to the left (it becomes a matrix), and the part with T_inf stays on the right
  (a load).
- An edge with no condition has `k dT/dn = 0`: insulated, or a symmetry line. This is natural, so
  nothing is assembled. The quarter models rely on it.
- The face sink has the same structure as Robin, but over the area instead of the edge.

## 3. Galerkin discretisation and element matrices

`T_h = sum_j N_j T_j` with the T3 or T6 shape functions of phase 2 (isoparametric, the same map).
With w = N_i:

$$M\dot{\mathbf T} + (K + H)\,\mathbf T = \mathbf F(t)$$

| Matrix | Element integral | Integrand degree (straight T3 / T6) | Rule |
|---|---|---|---|
| M (capacity) | `int rho c N_i N_j dA` | 2 / 4 | 3-point / 6-point (exact) |
| K (conductivity) | `int k grad N_i . grad N_j dA` | 0 / 2 | 1-point / 3-point (exact) |
| H (Robin) | `int_edge h N_i N_j ds` | 2 / 4 | 5-point Gauss (exact on straight edges) |
| H_f (face sink) | `int (2h_f/t_z) N_i N_j dA` | as M | as M |
| F (source) | `int q N_i dA` | depends on q | 6-point |
| F (Robin, flux, face) | `int h T_inf N_i ds`, `int q_in N_i ds`, `int (2h_f/t_z) T_inf N_i dA` | 1 / 2 | 5-point / 6-point |

A = K + H + H_f. This gives the same `ThermalOperators {M, A, F(t), Dirichlet}` as in 1D, so the theta
method (phase1.md section 5) is reused unchanged.

**Why no lumped mass for T6.** For a straight T6 the corner functions integrate to zero
(`int N_corner dA = 0`, while each mid-side function gets area/3). So row-sum lumping gives the
corners zero capacity, and M becomes singular. Phase 4 therefore uses the consistent mass for both
elements. The price is that the discrete maximum principle is not guaranteed for sharp initial data
(phase1.md section 5). The cases here start from smooth fields or from T = T0 with a smooth source,
so this does not matter.

Exactness checks that the unit tests use:
- `sum_i N_i = 1`, so `1^T M 1 = int rho c dA`, `1^T H 1 = int_edge h ds` and `1^T F_q = int q dA`.
- K annihilates constants: `K 1 = 0`.

Code: `src/thermal/heat2d.cpp`.

## 4. Discrete energy balance

Stored energy (J per metre of depth, or per plate with t_z): `E_h(t) = int rho c (T_h - T_ref) dA = 1^T M (T - T_ref 1)`.
This holds because `sum N_i = 1`. Multiply one theta step by `1^T`:

$$\mathbf 1^T M(\mathbf T^{n+1} - \mathbf T^n) = \Delta t\,[\theta\,P^{n+1} + (1-\theta)\,P^n],\qquad P = \mathbf 1^T\mathbf F - \mathbf 1^T A\,\mathbf T.$$

Since `1^T K = 0`, P contains only the source power, minus the Robin and face losses
`int h (T_h - T_inf) ds + int (2h_f/t_z)(T_h - T_inf) dA`, plus the imposed flux.
So the discrete energy changes by exactly the time-integrated discrete power. The only error is the
linear solver's residual. This is the 2D version of phase1.md section 7, and it is why P4-e (energy
to 1e-10) and P4-f (ii) can be asserted at round-off.

*What this balance does not prove.* Multiplying the theta step by `1^T` gives the balance for any M, A and F,
including wrong ones. So it checks the integrator and the bookkeeping, not the integrals. The numerics review
showed this by passing the check with a load tripled and an extra sink added to A only. Two checks can fail:
- **Load/matrix consistency (in the P4-f (ii) test).** With q = 0, a body at the ambient temperature must be
  in equilibrium: `A (T_inf 1) = F(t)` for spatially uniform T_inf(t). Here `K 1 = 0`, and H and H_f applied
  to `T_inf 1` must give exactly the Robin and face loads. The test asserts this on the
  curved Robin rim and the face sink to 1e-13 relative to `sum |A_ij| T_inf`. A load with the wrong h or the
  wrong 2/t_z fails it at order 1.
- **The integrals themselves:** the unit checks of section 3 and the exact steady slab (section 10).

With Dirichlet nodes the balance also contains the unknown reaction flux, so these checks use no
Dirichlet conditions.

**Uniform plate cooled through its faces only.** If T0 is uniform and the edges are insulated, the uniform
field stays uniform, because K acting on a constant is zero. Also `M 1 = rho c/(2h_f/t_z) * H_f 1`,
since both use the same N_i N_j integral. Each step then multiplies `T - T_inf` by

$$G = \frac{1 - (1-\theta)\lambda\Delta t}{1 + \theta\lambda\Delta t},\qquad \lambda = \frac{2h_f}{\rho c\,t_z},$$

which P4-f (iii) checks exactly. That is the amplification factor of phase1.md section 5, now for a real mode.

*Floating point limits this check, not the method.* The first version used T_inf = 300 K, T0 = 400 K and
k = 400 W/(m K), and missed 1e-10 by up to 6e-6. Two causes, both measured:
- After 20 steps the excess is `G^20 x 100 K = 0.03 K`, carried by T of about 300 K. Forming `T - T_inf`
  loses about 4 digits.
- The uniform mode does not involve conduction (`K 1 = 0`), but with `k dt/(rho c h^2) = 5.6e4` the matrix
  `M + theta dt A` is dominated by K. Round-off in `K 1`, of size eps times the entries of K, then enters the
  uniform mode. The error grew in proportion to k (k = 0.4, 4, 40, 400).

The test therefore runs in excess temperature (T_inf = 0) with k = 0.4 W/(m K). The discrete operator for
the mode is the same, and G^n does not depend on either choice.

## 5. Time integration

Unchanged from phase 1: `(M + theta dt A) T^{n+1} = (M - (1-theta) dt A) T^n + dt [theta F^{n+1} + (1-theta) F^n]`,
factorised once per (A, dt). A change of operators, such as the demo going from adiabatic during the pulse to
convective afterwards, builds a new integrator. Stage ends fall on step boundaries, and a stage evaluates
its own source at both ends of each of its steps. So the pulse stage sees the source "on" at every t_n
including t_n = tau, and the cooling stage sees it "off" from its first step on. The deposited energy is then
`q tau` exactly for any theta. Without stages, a switch-off inside a step would give that step a theta-weighted
fraction of the power: the energy would be wrong by up to `q dt`. The demo still uses BE: the jump in the
source and the change of operators excite the high-frequency modes, which CN does not damp (its amplification
tends to -1 for large lambda dt, phase1.md section 5) and which then ring in the stress.

Code: `src/thermal/theta_integrator.cpp` (unchanged), `src/io/case_thermomech2d.cpp` (stages).

## 6. Thermal strain in 3D

A free isotropic element heated by dT expands equally in all directions:
`eps_th = alpha dT I` (no shear). Hooke's law acts on the elastic part of the strain:

$$\boldsymbol\sigma = \mathbb C : (\boldsymbol\varepsilon - \alpha\,\Delta T\,\mathbf I) = \lambda\,\mathrm{tr}\,\boldsymbol\varepsilon\,\mathbf I + 2\mu\boldsymbol\varepsilon - (3\lambda + 2\mu)\,\alpha\Delta T\,\mathbf I,\qquad 3\lambda + 2\mu = \frac{E}{1-2\nu} = 3K_{bulk}.$$

Two consequences used below:
- A uniform dT in an unconstrained body gives `eps = alpha dT I` and zero stress.
- If expansion is prevented completely, `sigma = -E alpha dT/(1-2nu) I`. That is the largest possible
  thermal stress, and it blows up as nu -> 1/2.

## 7. Plane stress (thin plate)

`sigma_zz = sigma_xz = sigma_yz = 0`. Solve `sigma_zz = 0` for eps_zz:

$$\varepsilon_{zz} = -\frac{\nu}{1-\nu}(\varepsilon_{xx}+\varepsilon_{yy}) + \frac{1+\nu}{1-\nu}\,\alpha\Delta T.$$

Substitute back. The in-plane stress is the plane-stress D of phase 2 applied to the in-plane elastic
strain `eps - alpha dT m`:

$$\boldsymbol\sigma = D_{ps}(\boldsymbol\varepsilon - \alpha\Delta T\,\mathbf m) = D_{ps}\boldsymbol\varepsilon - \beta_{ps}\,\Delta T\,\mathbf m,\qquad D_{ps}\mathbf m = \frac{E}{1-\nu^2}(1+\nu)\,\mathbf m = \frac{E}{1-\nu}\mathbf m,\quad \beta_{ps} = \frac{E\alpha}{1-\nu}.$$

eps_zz is output only. It does not enter the 2D solve.

## 8. Plane strain (long body)

`eps_zz = 0`, so the out-of-plane expansion is prevented and that direction carries stress:

$$\boldsymbol\sigma = D_{pe}\boldsymbol\varepsilon - \beta_{pe}\,\Delta T\,\mathbf m,\qquad \beta_{pe} = 3K_{bulk}\alpha = \frac{E\alpha}{1-2\nu},$$

$$\sigma_{zz} = \lambda(\varepsilon_{xx}+\varepsilon_{yy}) - \frac{E\alpha\Delta T}{1-2\nu} = \nu(\sigma_{xx}+\sigma_{yy}) - E\alpha\Delta T.$$

The last equality uses `lambda = 2 nu (lambda + mu)` and `3K(1 - 2nu) = E`. It is the phase 2 formula plus
`-E alpha dT`, which von Mises must include.

Equivalent form: `sigma = D_pe (eps - (1+nu) alpha dT m)`. In plane strain the effective in-plane free
expansion is `(1+nu) alpha dT`, larger than alpha dT, because the prevented z expansion squeezes
the material sideways through Poisson's effect.

Both cases fit one form:

$$\boldsymbol\sigma = D\boldsymbol\varepsilon - \beta\,\Delta T\,\mathbf m,\qquad \beta = \begin{cases}E\alpha/(1-\nu) & \text{plane stress}\\ E\alpha/(1-2\nu) & \text{plane strain.}\end{cases}$$

Code: `src/coupling/thermal_strain.cpp` (`thermal_stress_coefficient`, `sigma_zz`).

## 9. The thermal load vector

Put section 8 into the virtual work of phase2.md section 2:

$$\int_\Omega \boldsymbol\varepsilon(w)^T D\,\boldsymbol\varepsilon(u)\,dA = \int_\Omega \boldsymbol\varepsilon(w)^T\,\beta\Delta T\,\mathbf m\,dA + \int_{\Gamma_t} w\cdot\bar t\,ds.$$

So, with `eps(w) = B w_e`:

$$K\mathbf u = \mathbf f_t + \mathbf f_{th},\qquad \mathbf f_{th}^e = \int_{\Omega_e} B^T\mathbf m\,\beta\,\Delta T_h\,dA,\qquad \Delta T_h = \textstyle\sum_a N_a T_a - T_{ref}.$$

Plane stress and plane strain differ only by beta (and by D on the left).

**What B^T m is.** Row 2a of B^T m is `dN_a/dx` and row 2a+1 is `dN_a/dy`. So the nodal thermal force is
`f_a = int beta dT grad N_a dA`. Integrating by parts:

$$\int_\Omega \beta\Delta T\,\nabla N_a\,dA = \int_{\partial\Omega} N_a\,\beta\Delta T\,\mathbf n\,ds - \int_\Omega N_a\,\nabla(\beta\Delta T)\,dA.$$

This is Duhamel's analogy: a temperature field acts like a body force `-grad(beta dT)` plus an outward
normal traction `beta dT` on the boundary. A uniform dT gives only the boundary push.

**Quadrature.** For T3, B is constant and dT_h is linear, so the integrand has degree 1 and the 1-point rule is
exact. For straight T6, B is linear and dT_h quadratic, so the integrand has degree 3. The 3-point stiffness rule
(degree 2) is not enough, so the 6-point degree-4 rule is used.

With a uniform dT, P4-b holds to round-off on any mesh, curved T6 included. The exact solution
`u* = alpha_eff dT x` is linear, so it lies in the element space, and `B u* = alpha_eff dT m` is constant
with `D m alpha_eff = beta m`. Hence `K u* = int B^T (beta dT m) dA`, integrated with the stiffness rule, and
`f_th` is the same integral with the thermal-load rule. Both integrands are `beta dT adj(J)^T dN/dxi`. That
is a polynomial in xi even on a curved cell: degree 0 for T3, degree 1 + 1 = 2 for T6, because adj(J) is
linear in xi and dN/dxi is too. The `1/det J` in `dN/dx` cancels against the `det J` of `dA`. Both rules are
exact for it, so `f_th = K u*` to round-off and the stress vanishes.

**Stress recovery.** `sigma_h(xi) = D B(xi) u_e - beta dT_h(xi) m`, at Gauss points, at any point and in
the error norms. Forgetting the second term is the classic bug: in plane stress, a free, uniformly heated
body would then show a stress of `E alpha dT/(1-nu)`.

**Out-of-plane stress and von Mises.** In plane strain, von Mises needs
`sigma_zz = nu (sigma_xx + sigma_yy) - E alpha dT_h` (section 8), with dT_h evaluated at the same point as the
in-plane stress. The phase 2 formula `nu (sigma_xx + sigma_yy)` silently drops the thermal term. For a hot
spot that term dominates: at the centre it is twice the in-plane stress (section 11). For nodal output,
sigma_zz is computed from the nodal-averaged in-plane stress and the nodal dT. Because dT_h is continuous at
a node and the average is linear, this equals the average of the per-cell sigma_zz values. Test:
`tests/phase4/thermal_output.cpp`, which fails with the phase 2 formula (relative error 0.7 for the held bar,
1.0 for the free body).

**Dependence on nu.** With in-plane loads, the Airy stress function of a plane-stress body satisfies
`lap^2 phi = -E alpha lap T`. If the boundary is traction-free (symmetry lines are equivalent) and the body is
simply connected, nothing else depends on nu. So the plane-stress stresses in the disk (P4-c/d) depend on
`E alpha` only. In plane strain the same equation holds with `E -> E/(1-nu^2)` and `alpha -> (1+nu) alpha`.
So `E alpha -> E alpha/(1-nu)`: the in-plane stresses are the plane-stress ones times `1/(1-nu)`, and
`sigma_zz` follows from section 8. Displacements depend on nu in both cases. So in the plane-strain demo the
assumed nu of graphite does matter, through 1/(1-nu). The demo reports the sensitivity (section 11).

Code: `src/coupling/thermal_strain.cpp` (`element_thermal_load`, `assemble_thermal_load`),
`src/mechanics/stress_recovery.cpp`.

## 10. Analytic solutions used by the tests

**Held bar (P4-a).** u_x = 0 at both ends and free in y. Then u_x is 0 everywhere, so eps_xx = 0, and
sigma_yy = 0 fixes eps_yy:
- Plane stress: `eps_yy = (1+nu) alpha dT`, `sigma_xx = -E alpha dT`.
- Plane strain: `eps_yy = (1+nu)/(1-nu) alpha dT`, `sigma_xx = sigma_zz = -E alpha dT/(1-nu)`.

Both are linear displacement fields, which T3 and T6 reproduce exactly.

**Free body (P4-b).** `u = alpha_eff dT (x - x_0)`, with `alpha_eff = alpha` (plane stress) or `(1+nu) alpha`
(plane strain). In-plane stress is 0. In plane strain, `sigma_zz = -E alpha dT`. The supports
(rollers on the symmetry lines x = 0 and y = 0 of the quarter plate in `thermal_free_body.cpp`; a roller
on the whole left edge plus the corner held in y in `thermal_output.cpp`) remove the rigid-body motion
without restraining an expansion about the origin.

**Thin disk (P4-c, P4-d).** Axisymmetric plane stress, a solid disk of radius b with a free edge
(Timoshenko and Goodier, Theory of Elasticity, section on thin circular disks):

$$\sigma_r = \frac{E\alpha}{b^2}\int_0^b T r\,dr - \frac{E\alpha}{r^2}\int_0^r T r\,dr,\qquad
\sigma_\theta = \frac{E\alpha}{b^2}\int_0^b T r\,dr + \frac{E\alpha}{r^2}\int_0^r T r\,dr - E\alpha T.$$

With `T = dT (1 - r^2/b^2)`: `int_0^r T r dr = dT (r^2/2 - r^4/(4 b^2))`, `int_0^b = dT b^2/4`. So

$$\frac{\sigma_r}{E\alpha\Delta T} = -\frac14\Big(1 - \frac{r^2}{b^2}\Big),\qquad \frac{\sigma_\theta}{E\alpha\Delta T} = -\frac14\Big(1 - \frac{3r^2}{b^2}\Big).$$

Displacement, from `u_r = r eps_theta = r [(sigma_theta - nu sigma_r)/E + alpha T]`:

$$u_r = \frac{\alpha\Delta T\,r}{4}\Big[(3+\nu) - (1+\nu)\frac{r^2}{b^2}\Big],\qquad \mathbf u = u_r\,(x, y)/r.$$

Check at r = 0: `eps = alpha dT (3+nu)/4 = (-(1-nu)/4 + 1) alpha dT`. Both forms agree.
Cartesian stress: `sxx = sr cos^2 + st sin^2`, `syy = sr sin^2 + st cos^2`, `sxy = (sr - st) sin cos`.
The nodal interpolant of T is used, just as a solver-computed T would be. Its error is O(h^{p+1}) in L2,
so it does not reduce the expected orders (P4-g: L2 order p+1, energy order p).

**Gaussian spreading (P4-e).** For `dT/dt = a lap T` in the plane, starting from `T0 exp(-r^2/(2 s0^2))`:

$$T(r, t) = T_0\,\frac{s_0^2}{s^2}\,\exp\!\Big(-\frac{r^2}{2s^2}\Big),\qquad s^2 = s_0^2 + 2at.$$

The energy `rho c T0 2 pi s0^2` is constant, and the variance per coordinate is
`int x^2 T / int T = s^2`. In fact `d/dt int x^2 T dA = a int x^2 lap T dA = 2a int T dA`
(integrate by parts twice; the boundary terms vanish), so the variance grows like `2at` for any
initial shape.

*A caution about what P4-e proves.* On a straight-sided T6 mesh, x^2 is in the finite-element space.
Taking w = x^2 in the semi-discrete equation then reproduces the identity above exactly. The only
leftover is the boundary term at x = L, which is about exp(-18) with L = 6 s_end. So for T6 the
variance check tests little beyond the initial interpolation and the energy bookkeeping. For T3 it is
a genuine test. P4-h (L2 error order of the whole profile) is the real accuracy check of the 2D heat
solver. That is why it is proposed.

**Separating time and space error in P4-h.** The error at t_end has two parts:
`e = e_space + e_time` with `||e_space|| ~ C_s h^{p+1}` (L2, element degree p) and, for Crank-Nicolson,
`||e_time|| ~ C_t dt^2`. To measure the spatial order, the time part must shrink at least as fast as the
spatial one. Otherwise the fit sees `dt^2` and reports a mixture. Choosing `dt_l = dt_0 (h_l/h_0)^{(p+1)/2}`
gives `dt_l^2 ~ h_l^{p+1}`, so the ratio `e_time/e_space` stays roughly constant across levels. The fitted
slope is then the combined space-time order along this path. It equals `p + 1` only if both parts have their
theoretical orders, so a wrong spatial order still shows (a T6 spatial order of 2.5 would give about 2.6).
For T3 this is `dt ~ h`, and for T6 `dt ~ h^{3/2}`. The test chooses `dt_0` once (10 steps on level 0) from this
argument, not from the result. It estimates the time part on level 2 as `(4/3) ||T(dt) - T(dt/2)||`. This is
the norm of the difference; a difference of norms would underestimate it, and the first version did. The
time part is 2% of the error for T3 but about half for T6. The numerics review measured 0.45, 0.47, 0.50,
0.52, 0.53 on levels 0-4, so the ratio is not exactly constant. Its Richardson-corrected, space-only orders
are 3.03 (T6) and 1.95 (T3). The meshes halve h nominally, so the step counts are `N_l = N_0 2^{l (p+1)/2}`,
rounded.

**Steady slab with source and convection (P4-f i).** On [0, L], with uniform q, `k T'' + q = 0`, and Robin h at
both faces:

$$T(x) = T_\infty + \frac{qL}{2h} + \frac{q\,x(L-x)}{2k}.$$

Check at x = 0: `-k dT/dn = k T'(0) = qL/2 = h (T(0) - T_inf)`. The profile is quadratic, so it lies in the T6 space. The source and the
Robin terms are integrated exactly on straight edges, so the T6 solution equals it to round-off.
Any error in H, the Robin load or the source quadrature would show up.

## 11. The demo model (D4) and its assumptions

This is a demonstration, not a validation. A Gaussian beam of rms size s travels along z through a long
rectangular block, and the energy deposit is taken as uniform along z: `E_l` joules per metre of beam
path during a pulse of length tau. The 2D model is the block's cross-section:

$$q(x, y, t) = \frac{E_l}{2\pi s^2\,\tau}\exp\!\Big(-\frac{r^2}{2s^2}\Big),\quad 0 \le t \le \tau;\qquad q = 0 \text{ after.}$$

Real deposits vary along z (shower build-up) and come from FLUKA or Geant4. The Gaussian is a stand-in.

**Why plane strain, not plane stress.** Plane stress needs the thickness along the beam to be much
smaller than the in-plane length scale of T, which is s. For a block several centimetres long and
s = 0.5 to 4 mm, the opposite holds. Away from the block's ends, a slice is in *generalised* plane
strain: `eps_zz = eps_0`, a constant set by zero net axial force. Equilibrium of a traction-free cross-section
gives `int (sigma_xx + sigma_yy) dA = 0`. So `int sigma_zz dA = 0` requires `eps_0 = alpha mean(dT)`, and
`sigma_zz` shifts by the uniform amount `E alpha mean(dT)`. For the deposit below in a 40 mm x 40 mm
copper section, `mean(dT) = E_l/(rho c A) = 0.18 K` against a peak of 46 K at s = 1 mm. Plane strain
(`eps_0 = 0`) therefore differs from generalised plane strain by a uniform `sigma_zz` shift of
`E alpha mean(dT)`. Relative to the peak `|sigma_zz| = E alpha dT0/(1-nu)` this is
`(1-nu) mean(dT)/dT0 = (1-nu) 2 pi s^2/A`, using `dT0/mean(dT) = A/(2 pi s^2)` (A = full section area).
It does not depend on E_l or E. For copper (nu = 0.35, A = 40 mm x 40 mm): 0.06% at s = 0.5 mm, 0.26% at
1 mm (0.65 x 0.18 K/46.2 K), 1.0% at 2 mm and 4.1% at 4 mm. Graphite (nu = 0.2) is about 20% higher.
Plane strain is also exact for a block whose ends are held. **The real state is generalised plane strain.**
The demo uses plane strain and prints this correction for each s. The demo uses plane strain and states this. A plane-stress run would
model a thin foil (t_z much smaller than s), a different device.

**Setup (approved: plane strain, annealed copper).**

| Item | Value | Reason |
|---|---|---|
| Cross-section | 40 mm x 40 mm, quarter model 20 x 20 mm, symmetry on x = 0 and y = 0 | (s/W)^2 <= 0.04 for s <= 4 mm |
| Deposit | E_l = 1 kJ/m (for example 10 J over 10 mm), tau = 10 us, same for both materials | stand-in; puts annealed Cu just above yield (74 vs 69 MPa) at s = 1 mm |
| Beam size sweep | s = 0.5, 1, 2, 4 mm | mesh size at the spot <= s/10 |
| Thermal BCs | pulse: adiabatic (insulated, no sink); cooling: Robin h = 10 W/(m^2 K), T_inf = 293.15 K on the outer edges | free convection in air, 2 to 25 W/(m^2 K) (Incropera et al., Table 1.1) |
| Initial / reference | T0 = T_ref = 293.15 K | stress-free at room temperature |
| Mechanics | plane strain, symmetry rollers on x = 0 and y = 0, outer edges free; linear elastic | one-way coupling, quasi-static |
| Time | pulse: BE, 10 steps of 1 us; cooling: BE, geometric steps to 1 s | stresses relax over s^2/a, about 10 ms |

**Material data.** Room-temperature values. Every row is cited, and gaps are marked.

| Property | Copper, OFE (UNS C10100, annealed) | Graphite, isotropic (SGL SIGRAFINE R4550) |
|---|---|---|
| rho (kg/m^3) | 8940 (0.323 lb/in^3) [CDA] | 1830 [SGL] |
| c (J/(kg K)) | 385 (0.092 Btu/(lb F)) [CDA] | 709: **assumption**, ideal graphite Cp = 8.517 J/(mol K) at 298.15 K, M = 12.011 g/mol [JANAF]; not measured on R4550 |
| k (W/(m K)) | 391 (226 Btu/(ft h F) at 68 F) [CDA] | 105 at 20 C [SGL] |
| alpha (1/K) | 16.9e-6 (9.4e-6/F, 68-212 F) [CDA] | 4.2e-6 (20-200 C) [SGL] |
| E (GPa) | 123.5 +/- 0.7, annealed polycrystal [L&N, Table 10]. CDA lists 117 (17 000 ksi), a rounded engineering value | 11.5 (dynamic) [SGL] |
| nu | 0.350 +/- 0.009, annealed polycrystal [L&N, Table 10] | 0.2: **assumption**, not in the TDS; sensitivity 0.15 to 0.25 reported |
| Strength (MPa) | yield 69, annealed OS050 [CDA, via MakeItFrom]. Temper-dependent: 69 to 400 | flexural 60, compressive 130 [SGL]. No tensile value given |

- [CDA] Copper Development Association, alloy C10100 data sheet (copper.org). The values above are the CDA
  table as reproduced by distributors; copper.org itself could not be fetched when this was written.
- [SGL] SGL Carbon, SIGRAFINE R4550 technical data sheet, TDS R4550.00, 04 2019 ("typical average values").
- [L&N] H. M. Ledbetter and E. R. Naimon, Elastic properties of metals and alloys. II. Copper,
  J. Phys. Chem. Ref. Data 3 (1974) 897-935, doi:10.1063/1.3253150. Table 10 lists recommended
  room-temperature polycrystal values, averaged over annealed specimens.
- [JANAF] M. W. Chase, NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9 (1998),
  table C-002, C1(ref) (graphite), janaf.nist.gov/tables/C-002.txt.
- [Incropera] F. P. Incropera et al., Fundamentals of Heat and Mass Transfer, Table 1.1.

The final values go in `cases/beam_block/README.md` and in the case JSON, with these citations.

**Estimates the demo prints next to the FE result (sanity checks, not validation rows):**
- *Adiabatic pulse.* The diffusion length `sqrt(a tau)` for tau = 10 us is 34 um (Cu) and 28 um (graphite),
  much less than s. So the end-of-pulse field is the deposit itself: `dT0 = E_l/(2 pi s^2 rho c)` at the
  centre. The domain stores exactly `E_l/4` per metre (quarter model, section 4).
- *Quasi-static stress.* The longitudinal sound speed in copper is about 4.7 km/s, so a 1 mm spot is
  crossed in about 0.2 us, much less than tau. Inertia and stress waves are neglected (out of scope). This
  approximation gets worse for shorter pulses and smaller spots, and the README says so.
- *Centre stress, infinite plane strain.* From section 10 with b -> infinity and the factor 1/(1-nu) of
  section 9: `sigma_r = sigma_theta = -E alpha dT0/(2(1-nu))`, `sigma_zz = -E alpha dT0/(1-nu)`, so
  von Mises `= E alpha dT0/(2(1-nu))`. The finite section differs by terms of order (s/W)^2.

| s (mm) | Cu dT0 (K) | Cu von Mises at centre (MPa) | Graphite dT0 (K) | Graphite sigma_zz at centre (MPa) |
|---|---|---|---|---|
| 0.5 | 185 | 297 | 491 | -29.6 |
| 1 | 46.2 | 74.2 | 123 | -7.4 |
| 2 | 11.6 | 18.6 | 30.7 | -1.9 |
| 4 | 2.9 | 4.6 | 7.7 | -0.5 |

- *Material comparison.* At equal deposit the centre stress scales as `E alpha/((1-nu) rho c)`, the
  thermal-shock figure of merit (lower is better), set against strength. Copper: 0.93 Pa m^3/J.
  Graphite: 0.047. Graphite gets 2.65 times hotter for the same deposit but carries about 20 times less
  stress. That, together with its high sublimation temperature, is why graphite is used for dumps and absorbers.
- *Strength check.* Copper is ductile, so use von Mises against yield. Above yield (s <= 1 mm here) the
  elastic result is only an indicator, because thermoplasticity is not built. Graphite is brittle, so use
  the maximum principal stress against a tensile strength and the minimum against the compressive strength (130 MPa).
  The TDS gives only flexural strength (60 MPa), which overestimates the tensile strength of graphite. So a
  ratio sigma_1/60 MPa is reported as optimistic, and no conversion factor is invented.
- *Validity of constant properties (the 780 K limit).* With constant properties, graphite at s = 0.5 mm reaches
  about 780 K. JANAF gives Cp = 19.83 J/(mol K) at 800 K against 8.517 at 298.15 K, so c is 2.3 times higher at
  780 K. Integrating the JANAF c(T) for the same deposit (6.4e8 J/m^3 at the centre) gives dT0 = 319 K
  (T = 612 K) instead of 491 K: constant c overestimates dT0 by 54% in that row. The rows s >= 1 mm
  (dT0 <= 123 K) are less affected but not exact. alpha and E of graphite also vary with T. All of this is out
  of scope (constant properties), and the demo flags the s = 0.5 mm row.
  Temperature-dependent properties are out of scope, and the demo says so for that row.
- **Not modelled:** temperature-dependent properties, radiation, real FLUKA/Geant4 energy maps, variation along
  the beam, plasticity under thermal load, fatigue from repeated pulses, dynamic stress waves.

## 12. The algorithm, step by step

**2D heat** (`src/thermal/heat2d.cpp`, `theta_integrator.cpp`):
1. Read the Gmsh mesh (T3 or T6, isoparametric).
2. Assemble M (capacity) and A = K + H + H_f. Each element integral uses a rule that is exact on straight
   cells (section 3).
3. Precompute one unit load vector per time function: Robin `∫h N ds`, flux `∫N ds`, face `∫(2h_f/t_z) N dA`.
   `F(t)` is then `Σ value(t) × unit vector`, plus the source, which is re-integrated at every t.
4. Factorise `M + θ dt A` once per (operators, dt), with LDLT. Step:
   `(M + θ dt A) Tⁿ⁺¹ = (M − (1−θ) dt A) Tⁿ + dt [θ Fⁿ⁺¹ + (1−θ) Fⁿ]`, Dirichlet by elimination.

**Thermal stress** (`src/coupling/thermal_strain.cpp`, `src/mechanics/elastic2d.cpp`):
1. Interpolate the nodal T onto the mechanical mesh. It is the same mesh, so this is the identity.
2. `f = f_traction + Σ_e ∫ Bᵀ m β ΔT_h dA`: 1 point for T3, 6 points for T6 (section 9).
3. Solve `K u = f` with the displacement constraints eliminated.
4. Recover `σ = D B u − β ΔT_h m`. In plane strain also `σ_zz = ν(σ_xx + σ_yy) − EαΔT_h`, then von Mises
   and principal stresses.

**Demo D4** (`src/io/case_thermomech2d.cpp`), for each beam size s:
1. Pulse stage: adiabatic, source on, 10 BE steps.
2. Check: stored energy `4·1ᵀM(T−T0)` against deposited `4 τ 1ᵀF_q`.
3. Elastic solve at the end of the pulse, plus the ν re-solves for graphite.
4. Cooling stages: Robin on the outer edges, source off, BE with a fixed step per decade (10 µs to 1 s), and
   an elastic solve at each stage end.
5. Print the comparison with the section 11 estimates.

## 13. Validation cases: what each proves and which bug it would catch

| Test | What it proves | A bug it would catch |
|---|---|---|
| P4-a held bar (both planes) | β and the recovery term together give the uniaxial-strain state | wrong β (for example using Eα in plane strain); a missing −βΔT m in recovery |
| P4-b free body (curved T6) | f_th = K u* exactly, so free expansion is stress-free | thermal load integrated with a different geometry or rule than K; a sign error in f_th |
| thermal_output (P4-a/b through the case CSV) | σ_zz and von Mises in the written output carry −EαΔT | the Phase 2 σ_zz formula in output code (it did catch this) |
| P4-c/d thin disk | a non-uniform ΔT is mapped into the right stress field | f_th built from nodal T instead of ΔT_h at quadrature points; wrong polar transform |
| P4-g disk orders | thermal load and recovery are consistent to the element's order | a thermal load that is exact only for uniform ΔT, which would cap the order |
| P4-f (i) steady slab | H, the Robin load and the source vector are each right | a missing factor 2 in H or a wrong edge Jacobian (curved edges are tested in (ii)) |
| P4-f (ii) energy balance | the θ integrator's bookkeeping (an identity for any M, A, F) | a step that drops or double-counts a load; loss through the solver |
| P4-f (ii) consistency `A T_inf 1 = F` | Robin and face loads match their matrices, on a curved rim | a T_inf load with the wrong h or the wrong 2/t_z |
| P4-f (iii) face-cooled plate | the face sink has the M-structure and the integrator gives exactly G per step | H_f assembled with another rule than M; θ applied to the wrong side |
| P4-e Gaussian variance and energy | the transient conserves energy and spreads at rate 2a | a wrong k/ρc scaling, or leaking energy through the edges |
| P4-h Gaussian orders | the transient 2D heat solution converges at p + 1 | under-integrated K or M, wrong T6 node order in the heat assembly |
| D4 smoke | the demo runs and stores exactly the deposited energy | stage bookkeeping that drops or double-counts the pulse |

## 14. What went wrong while building this phase, and how it was caught

1. **The demo was first set up as plane stress.** The draft hit a 10 mm plate face-on. Plane stress needs
   t_z ≪ s, and here s = 0.5 to 4 mm. This was caught on review of the derivation, before any demo code. The
   fix is a long block in plane strain, which the user approved, with the generalised plane strain shift
   quantified (section 11).
2. **σ_zz in the case output ignored the thermal term.** The draft had `sigma_zz_thermal` but the output used
   the Phase 2 `ν(σ_xx + σ_yy)`. Found by reading the code. Test written first: it failed with relative errors
   0.70 (held bar) and 1.00 (free body), then passed at 1e-13.
3. **P4-e at level 2 missed 0.5% by 0.05%.** The variance offset was constant in time and fell 4× per level,
   so it is the O(h²) interpolation of the initial Gaussian, not the dynamics. Moved to the finest level
   with the user's approval, and disclosed in the test.
4. **P4-f (iii) missed 1e-10 by up to 6e-6.** The error grew with k at fixed λ dt, and vanished in excess
   temperature. So the cause is floating point (section 4). Re-posed with the user's approval, and disclosed.
5. **A unit-test tolerance was scaled to the answer (hL) rather than to the magnitude of the cancelling
   terms (Σ|A_ij|).** Round-off in `1ᵀA1` is set by entries of size k, so the tolerance now scales with
   them.
6. **A hand calculation in this note was wrong.** The generalised plane strain shift was first written as
   0.4% at s = 1 mm. It is 0.26%; the (1−ν) factor was missing. It was caught by recomputing with a script
   before it reached the README. Lesson: numbers in prose get the same check as numbers in tests.
7. **The numerics review found claims stronger than the tests.**
   - P4-f (ii) is an identity of the scheme. It now also asserts `A T_inf 1 = F`, which can fail.
   - The P4-h time share was estimated with a difference of norms. It is 50% for T6, not 25%, and the asserted
     slope is stated as a combined space-time order.
   - Demo cooling values had a few percent of time-step error at 10 steps per decade. It is now 40 steps, and
     the values are within 1.4% of a run with 4× smaller steps.
   - The demo's σ_zz used T0 where it should use T_ref. This had no effect, since they are equal in both cases.
   - The temperature and flux boundary conditions had no test. Now exact 1D profiles check them.
8. **Two copper sources disagree on E:** CDA gives 117 GPa, Ledbetter and Naimon 123.5 ± 0.7 GPa. The
   demo uses the primary source for both E and ν, so the pair is consistent, and records both values.

## 15. Limitations and the assumptions behind them

- One-way coupling. The thermoelastic heating term `T0 3Kα tr(ε̇)` is dropped. Its relative size is the
  coupling parameter `δ = (3Kα)² T0/(ρc (λ + 2μ))`, about 0.02 for copper at 293 K with the demo's data
  (3K = E/(1−2ν) = 412 GPa, λ + 2μ = 198 GPa). So the error is about 2%.
- Quasi-static mechanics: no inertia and no stress waves (valid when τ is much longer than s/c_L).
- Constant, temperature-independent properties. The graphite row at s = 0.5 mm overestimates ΔT by 54%.
- Linear elasticity under thermal load: no thermoplasticity, so copper above yield is only flagged.
- Plane strain for the demo. Generalised plane strain, which is the real state of a long free block, shifts
  σ_zz uniformly by 0.06% to 5%.
- Consistent capacity only (T6 corners make row-sum lumping singular). There is no discrete maximum principle
  for sharp data, though none was violated in the demo.
- A Gaussian stand-in for the energy deposition. Real maps come from FLUKA or Geant4.

## 16. Five interview questions

1. *Why is a uniformly heated, unconstrained body stress-free in your code to 1e-14, even on curved T6?*
   The exact displacement `α_eff ΔT x` is linear, so it lies in the element space. K u* and f_th are then the
   same integrand, `β ΔT adj(J)ᵀ ∂N/∂ξ`. That is a polynomial in ξ even on a curved cell, because adj(J) is
   linear, and both quadrature rules integrate it exactly.
2. *Why is β different in plane stress and plane strain, and what does that do to σ_zz?* Plane stress lets
   the body expand freely in z, so only the in-plane constraint matters: β = Eα/(1−ν). Plane strain prevents
   the z expansion. That adds the Poisson squeeze, giving β = Eα/(1−2ν), and a stress
   `σ_zz = ν(σ_xx + σ_yy) − EαΔT` that must enter von Mises.
3. *How did you know the P4-f (iii) failure was not a bug?* Its exact answer does not depend on k or T_inf.
   The error grew in proportion to k at fixed λ dt, and it vanished when the problem was posed in excess
   temperature. A discretisation bug would not scale like that. Cancellation and conditioning do.
4. *Why does graphite get hotter than copper yet survive?* ΔT ∝ 1/(ρc), and graphite's ρc is 2.65 times
   lower. Stress ∝ EαΔT/(1−ν): graphite's Eα is about 40 times lower than copper's. The figure of merit
   Eα/((1−ν)ρc) is 20 times lower. Its strengths are of the same order as copper's yield.
5. *What would change with a real FLUKA map?* The source term becomes interpolated data, not a Gaussian. The
   deposit varies along z, so a 2D section is valid only where the gradient along z is weak; otherwise the
   model needs 3D or generalised plane strain per slice. The peak ΔT would be higher, so the properties would
   need to depend on temperature.
