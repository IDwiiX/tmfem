# Phase 3: elastoplasticity (J2 plasticity with isotropic hardening)

This note derives everything the phase 3 code does. As in phase 2, each section ends with the place
in the code that implements it, so the derivation and the code can be checked against each other.

Contents:
- 0: what is computed
- 1-12: derivation
- 13: algorithm
- 14: validation cases
- 15: what went wrong
- 16: limitations
- 17: interview questions
- Model answers: at the very end

## 0. What this phase computes

Phase 3 adds a nonlinear material to the phase 2 solver: small-strain von Mises (J2) plasticity with
isotropic hardening. The hardening curve comes from the tensile test in the source. There are two
versions:
- the measured true stress-strain points, interpolated linearly ("tabulated")
- the Ludwik law `sigma_y = sigma_y0 + K alpha^m`

It provides four capabilities:
- **A material-point driver.** It runs one material point under any mix of prescribed strains and
  stresses. A uniaxial tensile test is one such mix.
- **A parameter fit.** It identifies the Ludwik parameters from the data, with an uncertainty estimate.
- **A global Newton-Raphson solver with load stepping.** It runs on the phase 2 meshes, in plane strain
  first and then in plane stress.
- **An elastoplastic hole plate.** This is the Kirsch geometry at a remote 140 MPa. The elastic peak
  of 420 MPa exceeds yield, so a plastic zone forms. It is compared with Neuber's rule.

## 1. Setting and notation

Small strain, with the strain split additively into elastic and plastic parts:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p,\qquad
\boldsymbol\sigma = \mathbb C : \boldsymbol\varepsilon^e = K_b\,\mathrm{tr}(\boldsymbol\varepsilon^e)\,\mathbf 1 + 2G\,\mathrm{dev}(\boldsymbol\varepsilon^e),$$

with `G = E / (2(1+ν))` and `K_b = E / (3(1−2ν))`.

Notation:
- The bulk modulus is called K_b so that K stays free for the Ludwik coefficient.
- `dev a = a − ⅓ tr(a) 1` and `s = dev σ` (the stress deviator).
- `‖a‖ = sqrt(a : a)`.
- The internal variables are the plastic strain εᵖ and the equivalent plastic strain α.

Data:
- E = 207.423 GPa and ν = 0.3 (`references/data/material_tensile.json`).
- Everything is in SI internally (Pa); MPa appear only in files and printouts.

## 2. The J2 yield function

$$f(\boldsymbol\sigma, \alpha) = q(\boldsymbol\sigma) - \sigma_y(\alpha) \le 0,\qquad q = \sqrt{\tfrac32\,\mathbf s:\mathbf s} = \sqrt{\tfrac32}\,\|\mathbf s\|.$$

- **Why the 3/2.** In uniaxial stress σ along y, `s = σ diag(−⅓, ⅔, −⅓)` and `s:s = ⅔ σ²`, so q = |σ|.
  The yield stress in f is then the yield stress of the tensile test.
- **Pressure-insensitive and isotropic.** f depends on σ only through s. In principal stress space the
  yield surface is a cylinder of radius `sqrt(2/3) σ_y` around the hydrostatic axis σ1 = σ2 = σ3.
  Metals yield by shear (dislocation glide), and hydrostatic pressure does not drive it. This is why
  J2 suits steel.
- **In components** (no out-of-plane shear):
  `q² = ½[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²] + 3σxy²`. This is the phase 2 `von_mises()`. One
  difference matters: in plasticity σ_zz comes from the material model. The elastic plane-strain
  formula `ν(σxx+σyy)` no longer holds once εᵖ_zz ≠ 0 (section 8).

Code: `src/mechanics/j2_plasticity.cpp` (`equivalent_stress`).

## 3. Flow rule, hardening variable, KKT conditions

**Associative flow.** The plastic strain rate is normal to the yield surface:

$$\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\frac{\partial f}{\partial\boldsymbol\sigma} = \dot\gamma\,\frac32\frac{\mathbf s}{q} = \dot\gamma\sqrt{\tfrac32}\,\mathbf n,\qquad \mathbf n = \frac{\mathbf s}{\|\mathbf s\|}.$$

The derivative: `∂(s:s)/∂σ = 2s`, because dev is a projector, so `∂q/∂σ = (3/2)(2s)/(2q) = (3/2) s/q`.

**Consequences:**
- **Plastic incompressibility.** n is deviatoric, so `tr ε̇ᵖ = 0`.
- **Equivalent plastic strain.** Define `α̇ = sqrt(2/3) ‖ε̇ᵖ‖`. Then `α̇ = sqrt(2/3) sqrt(3/2) γ̇ = γ̇`.
  The plastic multiplier *is* the rate of equivalent plastic strain.
- **Work conjugacy.** `σ : ε̇ᵖ = γ̇ sqrt(3/2) s : n = γ̇ sqrt(3/2) ‖s‖ = q α̇ = σ_y α̇` while yielding.
  So the plastic work per unit volume is `∫ σ_y dα`, and σ_y(α) is the tensile-test curve plotted
  against plastic strain (next point).

**Loading/unloading (Karush-Kuhn-Tucker) conditions** and consistency:

$$\dot\gamma \ge 0,\quad f \le 0,\quad \dot\gamma f = 0,\qquad \dot\gamma\,\dot f = 0.$$

- Either the stress is inside the surface and nothing flows (f < 0, γ̇ = 0), or it flows and stays on
  the surface (γ̇ > 0, f = 0, ḟ = 0).

**The uniaxial test.**
- In tension along y, `n ∝ diag(−½, 1, −½)`, so
  `εᵖ_yy = α` and `εᵖ_xx = εᵖ_zz = −α/2` (incompressible).
- Adding the elastic part:

$$\varepsilon_{yy} = \frac{\sigma}{E} + \alpha,\qquad \varepsilon_{xx} = \varepsilon_{zz} = -\nu\frac{\sigma}{E} - \frac{\alpha}{2}.$$

- **P3-b in numbers**, at σ = 457.2405 MPa:
  - `σ/E = 0.0022044`
  - from the table, `α = 0.031581 − 0.0022044 = 0.0293766`
  - so `ε_yy = 0.031581` and `ε_xx = ε_zz = −0.3·0.0022044 − 0.0146883 = −0.0153496`
  - these match Abaqus's −0.01535 and 0.03158
- The lateral strain is mostly plastic: half the axial plastic strain. The elastic Poisson part
  is only 0.00066.

## 4. Hardening laws

`σ_y(α)` with slope `H(α) = dσ_y/dα ≥ 0`. In uniaxial loading q = σ and α = εᵖ_axial, so a tensile
curve converts directly into a hardening curve:

$$\alpha_i = \varepsilon_{\text{true},i} - \sigma_i / E.$$

**Tabulated** (used for P3-b and P3-f).
- The points of `tensile_true.csv` from the yield point on:

  | α | σ (MPa) |
  |---|---|
  | 0 (computed −3.2e-6, clipped) | 310.5611 |
  | 0.0086197 | 394.9973 |
  | 0.0293766 | 457.2405 |
  | 0.0510362 | 490.9348 |
  | 0.0845296 | 525.4933 |
  | 0.1287043 | 557.4827 |
  | 0.1511749 | 570.6438 |
  | 0.1715201 | 581.3758 |

- σ_y is piecewise linear in α, so H is piecewise constant: 9796 MPa on the first segment down to
  527 MPa on the last.
- At a knot, H takes the right-hand segment's slope, which is the one that applies when α grows.
- Beyond the last point, the last slope is extended (never reached here).
- **Why α at the yield point is −3.2e-6:** the digitised elastic line has slope
  120.7504 / 0.000575 = 210.0 GPa, not 207.423 GPa. With E = 207.423 GPa the elastic strain of the
  yield point is slightly larger than the data's total strain. The 3e-6 offset is the data's own
  round-off. It is clipped to 0, not "fixed".

**Ludwik.**
$$\sigma_y = \sigma_{y0} + K\alpha^m,\qquad H = mK\alpha^{m-1}\ \ (\to\infty \text{ as } \alpha\to 0 \text{ when } m<1).$$
- The source writes the Hollomon law σ = K εⁿ, but it fits `ln(σ − σ_y) = (1/n) ln εᵖ + ln K`. That is
  Ludwik with `m = 1/n = 0.4077`, K = 559.27 MPa and σ_y0 = 310.86 MPa.
- The infinite initial slope is what lets Ludwik follow the sharp knee at yield. It also needs care
  in the return mapping (section 5).

**Linear** (unit tests only): `σ_y = σ_y0 + Hα`, which has closed-form answers.

**Why a Ludwik fit is 11% to 27% off in plastic strain at 457 MPa, though only 1.3% to 3% off in stress.**
- Near α = 0.03 the curve is flat: H is between 1.6 and 3 GPa, against E = 207 GPa.
- A stress misfit δσ of 6 to 14 MPa therefore becomes `δα = δσ / H ≈ 0.003 to 0.009`, which is 11% to
  27% of α = 0.029.
- Inverting a flat curve amplifies errors. That is why P3-b uses the tabulated law.

Code: `src/mechanics/hardening.cpp`.

## 5. Return mapping (backward Euler)

The model's job: given the total strain ε_{n+1} at a Gauss point (from the current global iterate)
and the converged state at the end of the last load step (εᵖ_n, α_n), return σ_{n+1}, the tangent and
the new state. The flow rule is integrated with backward Euler: the flow direction is evaluated at
the *end* of the step.

**1. Elastic trial (freeze plastic flow):**
$$\boldsymbol\sigma^{tr} = \mathbb C:(\boldsymbol\varepsilon_{n+1} - \boldsymbol\varepsilon^p_n),\quad
\mathbf s^{tr} = 2G\,\mathrm{dev}(\boldsymbol\varepsilon_{n+1} - \boldsymbol\varepsilon^p_n),\quad
q^{tr} = \sqrt{\tfrac32}\|\mathbf s^{tr}\|,\quad f^{tr} = q^{tr} - \sigma_y(\alpha_n).$$

The interface takes the total strain, not the increment. The two are identical:
`σ^tr = C:(ε_{n+1} − εᵖ_n) = σ_n + C:Δε`. The total form avoids accumulating round-off over steps.

**2. Elastic case.** If `f^tr ≤ 1e-12 σ_y(α_n)`, the trial is the answer: σ = σ^tr, the tangent is
C, and the state is unchanged. The small tolerance keeps a converged point that lies exactly on the
surface from counting as plastic.

**3. Plastic case.** Backward Euler gives:
$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\sqrt{\tfrac32}\,\mathbf n_{n+1},\qquad \alpha_{n+1} = \alpha_n + \Delta\gamma.$$

- Δεᵖ is deviatoric, so `s_{n+1} = s^tr − 2G Δεᵖ = s^tr − 2GΔγ sqrt(3/2) n_{n+1}`.
- Write both deviators as norm × direction:
  `(‖s_{n+1}‖ + 2G sqrt(3/2) Δγ) n_{n+1} = ‖s^tr‖ n^tr`.
- Both coefficients are positive, so **n_{n+1} = n^tr**. The stress returns along the radius of the
  cylinder: the "radial return".
- It also gives `‖s_{n+1}‖ = ‖s^tr‖ − 2G sqrt(3/2) Δγ`. Multiplying by sqrt(3/2):

$$q_{n+1} = q^{tr} - 3G\,\Delta\gamma.$$

**4. Consistency fixes Δγ.** The scalar equation is
$$r(\Delta\gamma) = q^{tr} - 3G\Delta\gamma - \sigma_y(\alpha_n + \Delta\gamma) = 0,\qquad r'(\Delta\gamma) = -(3G + H) < 0.$$

- **One root in a known bracket.**
  - `r(0) = f^tr > 0`.
  - `r(f^tr/3G) = σ_y(α_n) − σ_y(α_n + f^tr/3G) ≤ 0`, since σ_y does not decrease.
  - r is strictly decreasing, so there is exactly one root in `[0, f^tr/(3G)]`.
- **Linear hardening:** `Δγ = f^tr / (3G + H)` in one step.
- **General hardening:** Newton, `Δγ ← Δγ + r / (3G + H(α_n + Δγ))`, safeguarded by the bracket:
  - It starts at the upper end `f^tr/3G`, where H is finite even for Ludwik at α_n = 0. Starting at
    Δγ = 0 there would give H = ∞ and a zero step.
  - A Newton step that leaves the current bracket is replaced by bisection.
  - After each evaluation the bracket shrinks: if r > 0 the root lies to the right, otherwise to the left.
  - It stops when `|r| ≤ 1e-12 σ_y(0)`.
- **Convergence.**
  - For a concave σ_y (Ludwik with m < 1, and this table), r is convex. Newton from the right then
    lands left of or on the root and climbs monotonically.
  - On a single tabulated segment r is linear, so Newton is exact in one step.
  - If a step crosses a knot, the next one uses the new segment's slope.

**5. Update:**
$$\boldsymbol\sigma_{n+1} = \boldsymbol\sigma^{tr} - 2G\sqrt{\tfrac32}\,\Delta\gamma\,\mathbf n,\qquad
\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \sqrt{\tfrac32}\,\Delta\gamma\,\mathbf n,\qquad \alpha_{n+1} = \alpha_n + \Delta\gamma.$$

The pressure is unchanged: `p_{n+1} = p^tr = K_b tr ε_{n+1}`.

**Why it is exact for proportional loading.** If the trial direction does not change during loading
(uniaxial tension), backward Euler's n^tr is the true direction. Then one step or a hundred give the
same answer. That is why P3-b does not depend on the number of steps. Non-proportional paths (the
plane-strain bar, the hole plate) do depend on the step size, to first order.

Code: `src/mechanics/j2_plasticity.cpp` (`J2Plasticity::update`, `solve_plastic_multiplier`).

## 6. Voigt form, and where the factor 2 lives

- The code stores 4-vectors: strain `[εxx, εyy, εzz, γxy]` with **engineering** shear γ = 2εxy, and
  stress `[σxx, σyy, σzz, σxy]`.
- The zz component is explicit, because εᵖ_zz and σ_zz are nonzero in plane strain. xz and yz are
  zero in 2D.
- With m = [1, 1, 1, 0]ᵀ:

$$\mathbb I_{dev} = \mathrm{diag}(1,1,1,\tfrac12) - \tfrac13\,\mathbf m\mathbf m^T,\qquad
\mathbb C = K_b\,\mathbf m\mathbf m^T + 2G\,\mathbb I_{dev}.$$

- **The ½ in 𝕀_dev.** 𝕀_dev maps an engineering strain vector to its deviator in tensor components
  (stress-like). The xy entry must turn γ into εxy = γ/2. The ½ therefore puts G (not 2G) on the
  shear diagonal of C.
- **Check:** the {xx, yy, xy} rows and columns of C equal the phase 2 plane-strain D (unit test).
- **Norm of a stress-like vector:** `‖s‖² = sxx² + syy² + szz² + 2 sxy²` (sxy appears twice in the tensor).
- **n is stored stress-like**, so the contraction with a strain is a plain dot product:
  `n:dε = nxx dεxx + nyy dεyy + nzz dεzz + 2 nxy dεxy = nᵀ dε_eng`.
- **Plastic strain increment, engineering form:** `√(3/2) Δγ [nxx, nyy, nzz, 2nxy]`. Forgetting the 2
  makes the plastic shear strain half what it should be. Uniaxial tests cannot see this, because they
  have no shear. A unit test checks `εᵖ = ε − C⁻¹σ` in a state with shear, and that catches it.

Code: `src/mechanics/j2_plasticity.cpp` (`elastic_tensor4`, `deviator_projector`).

## 7. The consistent (algorithmic) tangent, derived in full

Global Newton needs `dσ_{n+1}/dε_{n+1}` of the *algorithm* in section 5, not of the continuum
equations. Differentiate the update with the converged state (εᵖ_n, α_n) fixed:

$$d\boldsymbol\sigma = \mathbb C:d\boldsymbol\varepsilon - 2G\sqrt{\tfrac32}\,\big(\mathbf n\,d\Delta\gamma + \Delta\gamma\,d\mathbf n\big).$$

**(i) dn.** `n = s^tr / ‖s^tr‖` and `ds^tr = 2G 𝕀_dev : dε`. The derivative of a unit vector removes
the component along itself:
$$d\mathbf n = \frac{1}{\|\mathbf s^{tr}\|}(\mathbb I - \mathbf n\otimes\mathbf n):d\mathbf s^{tr} = \frac{2G}{\|\mathbf s^{tr}\|}(\mathbb I_{dev} - \mathbf n\otimes\mathbf n):d\boldsymbol\varepsilon,$$
using `(n⊗n):𝕀_dev = n⊗n` because n is deviatoric.

**(ii) dΔγ.** Differentiate the consistency equation `r = 0`:
`dq^tr − 3G dΔγ − H dΔγ = 0`, with `H = H(α_{n+1})` and
`dq^tr = sqrt(3/2) n : ds^tr = 2G sqrt(3/2) n : dε`. So

$$d\Delta\gamma = \frac{2G\sqrt{3/2}}{3G + H}\,\mathbf n:d\boldsymbol\varepsilon.$$

**(iii) Combine.**
- `(2G sqrt(3/2))² = 6G²`.
- `2G sqrt(3/2) Δγ · 2G/‖s^tr‖ = 6G²Δγ/q^tr`, because `‖s^tr‖ = sqrt(2/3) q^tr`.

Hence

$$\boxed{\ \mathbb C^{alg} = \mathbb C - \frac{6G^2}{3G+H}\,\mathbf n\otimes\mathbf n - \frac{6G^2\Delta\gamma}{q^{tr}}\,(\mathbb I_{dev} - \mathbf n\otimes\mathbf n)
= K_b\,\mathbf 1\otimes\mathbf 1 + 2G\theta\,\mathbb I_{dev} - 2G\bar\theta\,\mathbf n\otimes\mathbf n\ }$$

$$\theta = 1 - \frac{3G\Delta\gamma}{q^{tr}},\qquad \bar\theta = \frac{3G}{3G+H} - \frac{3G\Delta\gamma}{q^{tr}}.$$

In Voigt form: `C_alg = K_b m mᵀ + 2Gθ 𝕀_dev − 2Gθ̄ n nᵀ`, with n stress-like.

**Properties.**
- **Symmetric**, because the flow is associative (∂f/∂σ is the flow direction). K_T is therefore
  symmetric, and SimplicialLDLT applies.
- **Continuum limit.** As Δγ → 0 it becomes the continuum elastoplastic tangent
  `C − 6G²/(3G+H) n⊗n`.
- **The extra term is the point.** At finite Δγ, the term `(6G²Δγ/q^tr)(𝕀_dev − n⊗n)` accounts for
  n^tr rotating as ε changes. With the continuum tangent instead, global Newton loses its quadratic
  rate and converges only linearly. The P3-d rate test exists to catch that.
- **Positive definite for H > 0.** In the three invariant subspaces the moduli are:
  - volumetric: 3K_b
  - deviatoric directions orthogonal to n: `2Gθ = 2G q_{n+1}/q^tr > 0`
  - along n: `2G(θ − θ̄) = 2G H / (3G + H) > 0`
- **Tangent at a Ludwik first yield.** At α_n = 0, `H(α_{n+1} = Δγ) = mKΔγ^(m−1)` is finite because
  Δγ > 0. Only Δγ = 0 exactly would be singular, and that case is elastic.
- **Checked numerically while planning.** A central finite difference of the return map (Ludwik, a
  non-proportional state with shear) agreed with this formula to 3e-11 relative. That check is now a
  unit test.

Code: `src/mechanics/j2_plasticity.cpp` (end of `J2Plasticity::update`), test `tests/phase3/return_mapping.cpp`.

## 8. Plane strain and plane stress

**Plane strain** (`ε_zz = γxz = γyz = 0`).
- The 3D state has only the four components [xx, yy, zz, xy], so the 3D algorithm runs exactly on
  4-vectors. Nothing is approximated.
- σ_zz and εᵖ_zz are nonzero and stay in the state.
- The 2D stress is rows {xx, yy, xy} of the 4-vector, and the 2D tangent is the {xx, yy, xy} block
  of C_alg. No condensation is needed, since dε_zz = 0.
- From `ε_zz = 0 = (σ_zz − ν(σxx+σyy))/E + εᵖ_zz`:

$$\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy}) - E\,\varepsilon^p_{zz},$$

  which tends to (σxx+σyy)/2 as incompressible plastic flow dominates. At the hole edge, σ_zz/σ_yy
  goes from ν = 0.3 (elastic) toward 0.5.

**Plane stress** (`σ_zz = 0`, ε_zz unknown).
- The 3D model wants all strains, but ε_zz is not known from the 2D displacement field. Each Gauss
  point solves for it.
- **Local Newton:**

$$\varepsilon_{zz} \leftarrow \varepsilon_{zz} - \sigma_{zz}(\varepsilon_{zz}) / C^{alg}_{zz,zz}.$$

  - Each evaluation is a full 3D return from the converged state (εᵖ_n, α_n).
  - It starts from the previous converged ε_zz (stored in the state) and stops when
    `|σ_zz| ≤ 1e-12 σ_y(0)`.
- **Condensed tangent.** Setting `dσ_zz = C_zβ dε_β + C_zz dε_zz = 0` gives
  `dε_zz = −C_zβ dε_β / C_zz`, so

$$D_{\alpha\beta} = C_{\alpha\beta} - \frac{C_{\alpha z}\,C_{z\beta}}{C_{zz}}\qquad(\alpha,\beta \in \{xx, yy, xy\}).$$

  This is exact because C_alg is the exact derivative at the converged local state.
- **Elastic check.** In the elastic limit this condensation reproduces the phase 2 plane-stress D
  (unit test).
- **Why this design.**
  - It wraps the validated 3D return for any model. Phase 4's thermal strain needs only an offset to
    ε, not new plasticity algebra.
  - The alternative, Simo and Taylor's projected plane-stress return, is faster but needs a separate
    derivation and separate tests.

Code: `src/mechanics/plane_reduction.cpp`.

## 9. Global Newton-Raphson with load stepping

**Residual** (unit thickness). The internal force minus the external force:

$$\mathbf R(\mathbf u) = \mathbf f_{int}(\mathbf u) - \lambda\hat{\mathbf f},\qquad
\mathbf f_{int} = \sum_e\sum_q B_q^T\,\boldsymbol\sigma\big(\boldsymbol\varepsilon(\mathbf u)_q;\ \text{state}_{q,n}\big)\,w_q\det J_q.$$

- λ is the load factor.
- f̂ is the traction load at λ = 1.
- Prescribed displacements are `u_c = λ ĝ`.
- In the elastic case `f_int = K u`, which recovers phase 2.

**Linearisation.**
$$K_T = \frac{\partial\mathbf f_{int}}{\partial\mathbf u} = \sum_e\sum_q B_q^T D^{alg}_q B_q\,w_q\det J_q.$$
This is the exact Jacobian of R only because D_alg is the consistent tangent of section 7.

**One load step** (λ_{k−1} → λ_k), with every Gauss point's converged state frozen.
- Start from `u⁽⁰⁾ = u_{k−1}`.
- At iteration i:
  1. Evaluate every Gauss point from `(ε(u⁽ⁱ⁾), state_n)`. Always use the *converged* state, never the
     previous iterate's. Otherwise the final answer would depend on the iteration path.
  2. Assemble R and K_T.
  3. Set `Δu_c = λ_k ĝ − u_c⁽ⁱ⁾`. It is nonzero only at i = 0.
  4. Solve `K_T,ff Δu_f = −R_f − K_T,fc Δu_c`, which is phase 2 elimination for an increment.
  5. Update `u⁽ⁱ⁺¹⁾ = u⁽ⁱ⁾ + Δu`. Evaluate again and measure `r = ‖R_f‖`.
- **Converged** when `r ≤ tol · F_ref`, with `F_ref = max(‖λ f̂_f‖, ‖f_int,c‖)`. f_int,c are the
  reactions, so the test still works for pure displacement control, where f̂ = 0. Default tol = 1e-10.
- Then commit: `state_n ← state_{n+1}` at every Gauss point.
- If 25 iterations are not enough, the solver throws. There is no hidden step cutting; a failure is
  reported, not papered over.
- **Iteration 0 is an elastic predictor.** At u_{k−1} with the converged states, every point has
  f ≤ 0, so its tangent is C.

**Why Newton converges quadratically, and when it does not.**
- If R is twice differentiable near the solution and K_T is its exact Jacobian, then
  `‖R⁽ⁱ⁺¹⁾‖ ≤ c ‖R⁽ⁱ⁾‖²`. The observed rate
  `p_i = ln(r_{i+1}/r_i) / ln(r_i/r_{i−1})` tends to 2.
- R is smooth while the active set is fixed (the same Gauss points plastic) and σ_y is smooth for
  α > 0 (Ludwik).
- When a Gauss point switches between elastic and plastic between iterates, R has a kink and the rate
  drops for that iteration.
- With the tabulated law H jumps at the knots, so the rate is printed but not asserted.

**How P3-d measures the rate.** This criterion was fixed before any run.
- **Set-up.** Kirsch T6 level 2, plane strain, Ludwik with the source parameters. Two steps:
  S = 0 → 100 MPa (elastic) and 100 → 140 MPa (one large plastic step). The large step gives many
  iterations in the asymptotic range.
- **Qualifying triples.** A triple of consecutive iterates counts only if:
  - the plastic set is identical at all three iterates, and
  - `r_{i+1} ≥ 1e-11 F_ref`. Round-off sits near 1e-15 to 1e-14 F_ref, so it moves p by less than 0.05.
- **Assertion.** At least one triple must qualify, and the last qualifying triple must have p ≥ 1.8.

**Linear solver.**
- K_T is symmetric positive definite for H > 0, so `SimplicialLDLT` is used.
- The sparsity pattern is analysed once; each iteration only refactorises.
- The phase 2 relative-pivot check guards against a missed rigid-body mode.

**Quadrature.**
- T3: one point. The strain is constant in a T3, so one point is exact even for a nonlinear material.
- T6: the 6-point degree-4 rule (SPEC §5) for both residual and tangent. They must use the same
  points, or K_T is not the derivative of R. The history variables live at those 6 points.
- The elastic stiffness from 6 points equals the 3-point one only on straight-sided T6 (a degree-2
  integrand). That is why the elastic-limit test compares on the straight bar.

Code: `src/mechanics/elastoplastic2d.cpp`, `src/fem/assembly.cpp` (`assemble_system`).

## 10. The material-point driver

The driver runs a single material point, without a mesh, under mixed control. Each component of
[xx, yy, zz, xy] is either strain-controlled (its value is prescribed) or stress-controlled (its
stress is prescribed, and the strain is unknown).

**Each step:**
- Set the strain-controlled components.
- Solve for the stress-controlled strains s by Newton, starting from the previous step's strain:
  `ε_s ← ε_s − C_ss⁻¹ (σ_s − σ_target,s)`.
- Every evaluation starts from the converged state of the previous step.
- C_ss is the stress-controlled block of C_alg, so this Newton is also quadratic.
- Stop at `max |σ_s − σ_target,s| ≤ 1e-4 Pa` (about 2e-13 relative).

**Paths used:**
- **Uniaxial stress** (the tensile test): σ_xx = σ_zz = σ_xy = 0. Either σ_yy is prescribed (P3-b)
  or ε_yy is (P3-c).
- **Plane-strain bar** (P3-e): ε_xx prescribed, ε_zz = 0, σ_yy = σ_xy = 0. This path is not
  proportional, because σ_zz/σ_xx changes from ν toward ½. The FE bar must therefore use the same
  increments to match.

**Why it matters.**
- It separates the constitutive model from the finite elements.
- If the FE bar (a uniform state) differs from the driver, the bug is in assembly or boundary
  conditions, not in the material.

Code: `src/mechanics/material_point.cpp`.

## 11. Parameter identification and its uncertainty

**Data.**
- The plastic points of the true curve: the 8 points with σ ≥ 310.5611 MPa, with α computed as in
  section 4. The first point is clipped to α = 0.
- The 0.000575 strain point is elastic and is not used.

**(a) Log-log fit with σ_y fixed** (reproduces REFERENCE P3-a).
- Linear regression of `ln(σ_i − σ_y) = ln K + m ln α_i` with σ_y = 310.5611 MPa.
- It uses the 7 points with α > 1e-6; the clipped point has ln 0.
- Result expected: m = 0.385, K = 548 MPa. The source's 0.4077 / 559.27 used a different choice of
  points.

**(b) Nonlinear least squares, all three parameters free** (the deliverable).

$$\min_\theta\ S(\theta) = \sum_i r_i^2,\qquad r_i = \sigma_i - \sigma_{y0} - K\alpha_i^m,\qquad \theta = (\sigma_{y0}, K, m).$$

- Residuals are in MPa: the misfit is measured in stress.
- Jacobian of the model, row i: `[1, α_i^m, K α_i^m ln α_i]`, with `α^m ln α → 0` at α = 0.
- Solved by Levenberg-Marquardt, `(JᵀJ + μ diag JᵀJ) δ = Jᵀr`, starting from the log-log estimate.
  μ shrinks after an accepted step and grows after a rejected one.
- It stops when no step lowers the sum of squares or the step is at round-off. The test then checks
  stationarity, `‖Jᵀr‖ ≤ 1e-10 ‖J‖‖r‖`; the measured value is 8e-12.

**Uncertainty**, all deterministic:
1. **Asymptotic covariance.**
   - `Cov = s² (JᵀJ)⁻¹`, with `s² = S/(n − p)`, n = 8 and p = 3.
   - Standard errors are the square roots of the diagonal.
   - 95% intervals: `±t(0.975, 5) SE = ±2.5706 SE`, from Student's t because n − p is small.
2. **Jackknife.**
   - Refit 8 times, each time leaving one point out.
   - `SE_jk = sqrt((n−1)/n Σ(θ_(i) − θ̄)²)`.
   - The jackknife makes no linearity assumption, and it shows which points the fit depends on.
3. **Confidence band on the curve (delta method).**
   - `Var σ̂(α) = gᵀ Cov g`, with `g = ∂σ̂/∂θ = [1, α^m, Kα^m ln α]`.
   - The band is `±t·sqrt(Var)`.

**What to expect** (from the planning check):
- σ_y0 = 309.1 ± 4.5 MPa, K = 520.9 ± 15.3 MPa, m = 0.362 ± 0.016, s = 4.6 MPa.
- The jackknife SEs are much larger: about 100 MPa, 50 MPa and 0.13.
- **Leaving out the yield point drops σ_y0 to 193 MPa.** There are no data between α = 0 and 0.0086,
  so the knee, and with it σ_y0, is fixed by a single point.
- This is a property of the data set, not of the code. A test cannot fix it; only a denser tensile
  record near yield can.
- A 2-parameter fit with σ_y0 fixed at 310.5611 MPa is reported beside it, as the better-posed
  problem.

**Honest reading.** n = 8 digitised points.
- The band reflects model misfit and the influence of individual points. It says nothing about
  specimen-to-specimen scatter, which one test cannot measure.
- E is taken as exact. An error in E shifts every α_i.

Code: `src/mechanics/tensile_fit.cpp`.

## 12. Neuber's rule at the hole

Neuber's rule estimates the notch-root stress and strain from the elastic solution. The product of
stress and strain at the root equals the elastic one:

$$\sigma\,\varepsilon = \frac{(K_t S)^2}{E},\qquad \varepsilon = \frac{\sigma}{E} + \alpha(\sigma)\quad(\text{same hardening curve}).$$

- **Solving it.** The left side increases with σ, so one root exists. It is found by bracketed
  bisection.
- **Where it applies.** It assumes a **uniaxial** stress state at the root. That holds at the edge of
  a hole in a thin plate (plane stress): σ_rr = 0 on the free surface, and σ_zz = 0.
- **Numbers** with Kt S = 420 MPa and the tabulated curve:
  - `(Kt S)²/E = 0.85044 MPa`
  - **σ_N = 321.3 MPa, ε_N = 0.002647** (α_N = 0.00110)
  - the elastic values would be 420 MPa and 0.002025
- **Plane strain.** The edge is biaxial there (σ_zz/σ_yy between 0.3 and 0.5). The uniaxial rule is
  then reported only with the label "not applicable: biaxial notch root".
- **It is an estimate.** It ignores the redistribution of load around the plastic zone. The
  difference from FE is reported, never forced.

**First yield, from the elastic Kirsch field at the edge** (σ_yy = 3S, σ_rr = 0):
- **Plane stress:** von Mises = 3S, so yield starts at S = 310.56/3 = 103.5 MPa.
- **Plane strain:** σ_zz = 0.3·3S, so von Mises = `3S sqrt(1 − ν + ν²) = 2.667 S`, and yield starts at
  S = 116.5 MPa.

So at S = 140 MPa the plate is plastic in both conditions, more so in plane stress.

Code: `src/mechanics/neuber.cpp`, test `tests/phase3/kirsch_plastic.cpp`.

## 13. The algorithm, step by step

1. **Read the case.** Mesh (.msh), E and ν, the plane condition, the hardening law, and the boundary
   conditions at load factor 1.
   - A tabulated law is built from the true curve: `α = ε − σ/E`, with the first point clipped to 0.
2. **Precompute per cell.** The dofs, and at each quadrature point (1 on T3, 6 on T6): x, `w det J`
   and B. Every point starts with the zero state (εᵖ = 0, α = 0, ε_zz = 0).
3. **For each load factor λ_k = k/N:**
   1. Linearise at `u_{k−1}` with the converged states, which gives the elastic predictor.
   2. Correct: solve `K_T,ff Δu_f = −R_f − K_T,fc Δu_c`. Δu_c applies the displacement increment, and
      only in the first correction.
   3. At every quadrature point compute `ε = B u_e`, then `plane_update`:
      - plane strain: the 3D return with ε_zz = 0
      - plane stress: the local ε_zz iteration around the 3D return
      - the 3D return itself: trial → yield check → scalar Newton for Δγ → update → C_alg
   4. Assemble f_int and K_T, and measure `‖R_f‖` against `F_ref = max(‖f_ext,f‖, ‖reactions‖)`.
   5. Repeat 2–4 until `‖R_f‖ ≤ 1e-10 F_ref`, then commit the states.
4. **Post-process.**
   - Quadrature-point stresses, σ_zz and α (CSV).
   - Nodal means of the cell-mean α and von Mises stress (VTU, plots only).
   - Total strain at a point, from the displacement field.

## 14. Validation cases: what each proves and what bug it would catch

| ID | Case | Proves | Would catch |
|---|---|---|---|
| P3-a | Log-log fit (σ_y fixed) reproduces 0.385 / 548 MPa; 3-parameter least squares with covariance and jackknife | The data conversion and the fit | E not applied in α = ε − σ/E, MPa/Pa slips, a wrong ln α derivative (synthetic recovery to 1e-8) |
| P3-b | Uniaxial stress to 457.24 MPa, tabulated | Return mapping in uniaxial stress: ε_yy is a round trip, and the lateral strains test ν and incompressibility | wrong plastic flow direction, a missing −α/2, a wrong elastic part (−νσ/E) |
| P3-c | Ludwik at the data strains; closed form and ±2s band | The driver with infinite initial slope. The band part is weak: with 8 points and 3 parameters no least-squares residual can exceed √5 s = 2.24 s, so it nearly always passes. The source fit misses it at 3 points (−2.95, −2.99, −2.26 s), printed. | Newton stuck at α = 0 (H = ∞), the closed form disagreeing with the 3D model |
| P3-d | Newton rate ≥ 1.8 on the FE bar (plane strain and plane stress) | The tangent along the loading path is the exact derivative. The bar is uniform, so global Newton reduces to a scalar Newton on σ_yy(ε_yy). That still contains the dn term, so a continuum tangent fails it. The full K_T is checked separately against central finite differences on a non-uniform plastic Kirsch state ("Global tangent matches finite differences of the internal force", 1e-6). | the continuum tangent instead of the algorithmic one (rate → 1), the wrong condensation in plane stress |
| P3-e | FE bar = material point, every step and point, 1e-6 | Assembly, history handling and BCs agree with the constitutive model | history updated from the iterate instead of the converged state, wrong quadrature weights, wrong reactions |
| P3-f | Hole plate at 140 MPa, report with Neuber | The whole solver on a curved, non-uniform problem | (report) equilibrium and yield consistency are asserted |
| units | FD tangents (local, condensed, global), elastic limits, the 6-point rule | Each piece in isolation | a missing factor 2 on the shear plastic strain (εᵖ = ε − C⁻¹σ), a wrong D in either plane condition |

## 15. What went wrong while building this phase, and how it was caught

1. **The P3-d set-up fixed in the plan produced no data.**
   - The plan asserted the rate on the Kirsch plate: T6 level 2, plane strain, one step from
     100 to 140 MPa.
   - Newton converged so fast that the criterion found no qualifying triple. The criterion needs three
     iterates with the same plastic set and the last residual ≥ 1e-11 F_ref.
     - Residuals: 0.29, 1.35e-2, 3.5e-4, 1.9e-6, 1.4e-10, then the floor at 7e-14.
     - The plastic set grew from 93 to 99 to 100 points and was fixed only from iterate 3.
     - So only two residuals lay above round-off with a fixed set.
   - Nothing pointed to a wrong tangent:
     - The global finite-difference tangent check passed.
     - The ratios r₊/r fell from 0.047 to 0.026, 0.0054 and 7.7e-5. That is superlinear; a continuum
       tangent would make them level off.
   - Changing the measurement after seeing a result is the user's call, so I asked. The user chose to
     assert on the uniform FE bar: every point is plastic from the first correction, so the active set is
     fixed by construction.
   - Result: p = 2.12 (plane strain) and 2.62 (plane stress). The Kirsch sequences are still printed and
     plotted as a diagnostic.
   - **Lesson:** a rate needs three iterates in the asymptotic range. With a consistent tangent that
     range is short, so the test problem must keep the active set fixed.
2. **A unit-test tolerance was tighter than its inputs.**
   - The plane-stress wrapper test fed the driver's strains to the wrapper and asked for |σ_xx| ≤ 1e-6 Pa.
   - The driver converges its lateral stresses only to 1e-4 Pa, and the wrapper converges σ_zz to
     1e-12 σ_y0 (3e-4 Pa). The observed σ_xx was 3e-5 to 2e-4 Pa, which is 5e-13 of σ_yy.
   - The bound was a specification error, so it became 1e-12 σ_yy, consistent with both tolerances.
     The solver did not change.
3. **A misleading printout.**
   - The rate is computed from raw residuals, but the first table printed r divided by each iterate's own
     F_ref. F_ref contains the reactions, which jump during a displacement step, so the printed ratios did
     not match the rates.
   - The table and the CSV now normalise by the final F_ref. The qualification test still uses each
     iterate's own F_ref, as fixed in the plan.
4. **Not a bug: the jackknife.**
   - The leave-one-out fits showed that the yield point alone pins σ_y0 (section 11).
   - This was reported rather than "fixed", because it is a property of an 8-point data set.
5. **Findings from the numerics review.**
   - **P3-f compared the wrong stress with Neuber.** The first table set the field maximum of σ_yy
     against Neuber's notch-root stress. That maximum lies at the plastic-zone boundary (r ≈ 1.17a in
     plane stress), where σ_xx > 0, not at the root.
     - As a result, FE stress looked 3% above Neuber while FE strain was below it, which is inconsistent
       at a uniaxial root.
     - The comparison now uses the quadrature point nearest the edge: 318.7 MPa against 321.3 MPa
       (−0.8%). The field maximum is reported separately, with its radius and no reference.
   - **The P3-d bar is uniform** and the P3-c band is nearly guaranteed. Both are now stated in section
     14 and the README, together with the global finite-difference tangent test that carries the full
     consistency claim.
   - **Clipping guard.** The tabulated law used to clip the first α to 0 whatever its value. It now
     throws if α < −1e-5, so a wrong E or σ_start can no longer hide (unit test with E = 150 GPa).
   - **The README "Error" column** showed differences from non-targets. Report-only rows now show the
     difference in parentheses, and rows without a reference show n/a.

## 16. Limitations and the assumptions behind them

- **Small strain** (additive split, no rotations). The tensile data reach a true strain of 0.17, but
  the hole plate stays below 0.003.
- **Isotropic hardening only.** There is no kinematic hardening, so no Bauschinger effect: the residual
  stress after unloading assumes reverse yielding at +σ_y(α). At −100 MPa residual stress (plane stress)
  no reverse yielding occurs, so this does not matter here.
- **Rate-independent, isothermal.** One material, homogeneous, unit thickness.
- **The hole plate uses the elastic Kirsch tractions** on the 5a boundary. Once a plastic zone exists,
  these approximate an infinite plate. The zone is small next to the 5a distance.
- **The zone size converges more slowly than the edge values.** In plane stress the extent r − a is
  0.159a, 0.167a, 0.175a and 0.174a on T6 levels 1 to 4: +4.6% from level 2 to 3, then −0.3% from 3 to 4.
  The edge strain changes by only 0.01% from level 2 to 3. It is measured at quadrature points, so its
  resolution is the local point spacing.
- **Hole-edge stresses come from the nearest quadrature point** (4e-5 m from the edge on level 3). The
  edge strain comes from the displacement field.
- **Load steps:** 28 equal steps. Halving them changes the reported values by less than 0.01%.
- **Tabulated hardening:** H jumps at the knots, so quadratic convergence is not guaranteed in general
  (P3-d asserts with Ludwik; the tabulated bar happened to give 2.13).
- **The Ludwik uncertainty** is a statement about 8 digitised points and the model misfit, not about
  material scatter.
- **Neuber** is applied only where its uniaxial assumption holds (plane stress). In plane strain it is
  printed with a "not applicable" label.

## 17. Five interview questions

1. Why is radial return exact for J2 plasticity with isotropic hardening, and what does "exact" mean
   here? When would you need more than one scalar equation?
2. Derive the extra term in the consistent tangent compared with the continuum tangent. What happens to
   Newton if you use the continuum tangent?
3. Your plane-stress implementation wraps a 3D model in a local iteration. Why is the condensed tangent
   `C_ab − C_az C_zb / C_zz` exact, and what would break if the local iteration were stopped early?
4. The Ludwik fit gives σ_y0 = 309 ± 5 MPa from the covariance but ±100 MPa from the jackknife. Which
   do you believe, and what experiment would settle it?
5. At the hole edge in plane stress, FE gives ε_yy = 0.00236 and Neuber gives 0.00265. Why does Neuber
   overestimate, and why does it not apply in plane strain?

---

## Model answers

Try the questions first. Then expand.

<details>
<summary>1. Radial return</summary>

- The backward-Euler equations make s_{n+1} parallel to s^tr: both are a positive multiple of n_{n+1}
  (section 5).
- So the direction is known before Δγ is, and the whole return reduces to one scalar equation
  `q^tr − 3GΔγ = σ_y(α_n + Δγ)`.
- "Exact" means the discrete backward-Euler equations are solved exactly (to round-off). It does not
  mean the continuum path is reproduced exactly; that is only true for proportional loading.
- More equations are needed for:
  - kinematic hardening with a nonlinear (Armstrong-Frederick) back stress, where the direction
    depends on the back stress
  - pressure-dependent yield (Drucker-Prager apex)
  - plane stress done directly, without the 3D wrapper
</details>

<details>
<summary>2. Consistent vs continuum tangent</summary>

- Differentiating `σ = σ^tr − 2G√(3/2)Δγ n` gives two terms: dΔγ (present in the continuum tangent)
  and Δγ dn.
- dn exists because the trial direction rotates when ε changes:
  `dn = (2G/‖s^tr‖)(𝕀_dev − n⊗n):dε`.
- The extra term is `−(6G²Δγ/q^tr)(𝕀_dev − n⊗n)`. It softens the deviatoric directions orthogonal
  to n by the factor θ = q_{n+1}/q^tr.
- With the continuum tangent, K_T is not the derivative of R. Newton still converges for small steps,
  but only linearly: the ratio r₊/r tends to a constant instead of to 0.
</details>

<details>
<summary>3. Plane-stress condensation</summary>

- At the converged local state σ_zz(ε_in, ε_zz) = 0. The implicit function theorem gives
  `dε_zz/dε_in = −C_zβ/C_zz`, so dσ_in/dε_in is the condensed matrix.
- It is exact only if:
  - C is the exact derivative of the return (consistent), and
  - the local equation is solved: σ_zz is zero to round-off.
- If the local iteration stops early, σ_zz ≠ 0, the state is not plane stress, and the condensed
  matrix is no longer the derivative of what the element sees. Global convergence degrades.
</details>

<details>
<summary>4. Covariance vs jackknife</summary>

- The covariance assumes the model is right, the errors are independent, and the problem is close to
  linear near the optimum.
- The jackknife shows that one point, the yield point, carries σ_y0: removing it lets σ_y0 slide to
  193 MPa while K and m compensate.
- So σ_y0 is not identified by the rest of the data. The honest answer: σ_y0 ≈ 309 to 311 MPa *given*
  the measured yield point, with the knee shape poorly constrained.
- A tensile test with dense sampling between ε = 0.15% and 1% would settle it.
</details>

<details>
<summary>5. Neuber</summary>

- Neuber assumes the local "strain energy" product σε equals the elastic one. That ignores the load
  redistribution through the plastic zone: the surrounding elastic material constrains the root.
- It is usually found to overestimate the notch strain. Here the FE edge strain is 11% lower.
- The FE value lies between the linear rule (elastic strain, 0.00203) and Neuber (0.00265). That
  ordering is common for notches, but this project shows it for one case, not in general.
- In plane strain σ_zz ≈ 0.3 to 0.5 σ_yy at the edge, so the root is biaxial. The uniaxial curve
  σ(ε) is then the wrong constitutive relation. Generalised Neuber (Hoffmann-Seeger, Glinka energy) is
  needed.
</details>
