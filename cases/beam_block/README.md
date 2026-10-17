# Demo D4: beam-heated block (not a validation)

A Gaussian beam of rms size s travels along z through a long block and deposits E_l = 1 kJ per metre of
path during a pulse of 10 µs. That is, for example, 10 J over 10 mm. The model is the 40 mm x 40 mm
cross-section, simulated as a quarter with symmetry. The pulse is adiabatic. Afterwards the outer edges
cool by convection (h = 10 W/(m² K), T_inf = 293.15 K). Each output time solves linear elasticity in
**plane strain** with the thermal strain (one-way coupling). Derivation and assumptions:
`docs/learning/phase4.md` section 11.

Run from the repository root:

```
./build/src/tmfem cases/beam_block/copper.json
./build/src/tmfem cases/beam_block/graphite.json
python3 scripts/plot_phase4.py
```

Meshes: `block.geo`, one T6 mesh per s with element size s/10 near the axis (`scripts/make_meshes.sh`).

## What this demo is not

- **Not a real energy map.** Real deposits come from FLUKA or Geant4. They vary along z (shower build-up) and
  are not Gaussian. The Gaussian here is a stand-in.
- **Not the exact mechanical state.** The real state of a long block is *generalised* plane strain: eps_zz is a
  constant set by zero net axial force. The difference from plane strain is a uniform sigma_zz shift of
  E alpha mean(dT), that is `(1 - nu) 2 pi s^2 / A` of the peak |sigma_zz|. For copper this is 0.06% at
  s = 0.5 mm, 0.26% at 1 mm, 1.0% at 2 mm and 4.1% at 4 mm (graphite: 0.08%, 0.31%, 1.3%, 5.0%). The demo
  prints the shift for each s. Plane strain is exact for a block whose ends are held.
- **No temperature-dependent properties, plasticity, dynamics or fatigue.** Above yield the elastic copper
  stresses are only an indicator. A 10 µs pulse is long compared with the 0.2 µs acoustic transit of a
  1 mm spot, so quasi-static stress is a fair approximation, but stress waves are not modelled.

## Material properties (room temperature)

| Property | Copper, OFE (UNS C10100), annealed | Graphite, isotropic (SGL SIGRAFINE R4550) |
|---|---|---|
| rho (kg/m³) | 8940 (0.323 lb/in³) [CDA] | 1830 [SGL] |
| c (J/(kg K)) | 385 (0.092 Btu/(lb °F)) [CDA] | 709 **assumption**: ideal graphite, Cp = 8.517 J/(mol K) at 298.15 K, M = 12.011 g/mol [JANAF]; not measured on R4550 |
| k (W/(m K)) | 391 (226 Btu/(ft h °F) at 68 °F) [CDA] | 105 at 20 °C [SGL] |
| alpha (1/K) | 16.9e-6 (9.4e-6/°F, 68-212 °F) [CDA] | 4.2e-6 (20-200 °C) [SGL] |
| E (GPa) | 123.5 ± 0.7, annealed polycrystal [L&N, Table 10]. CDA lists 117 (17 000 ksi), a rounded engineering value | 11.5, dynamic modulus [SGL] |
| nu | 0.350 ± 0.009, annealed polycrystal [L&N, Table 10] | 0.2 **assumption**, not in the TDS. The demo reruns nu = 0.15 and 0.25 |
| Strength (MPa) | yield 69, annealed OS050 [CDA, via MakeItFrom]. Temper-dependent: 69 to 400 | flexural 60, compressive 130 [SGL]. No tensile value is given; flexural strength overestimates tensile strength, so the tensile check is optimistic |
| Convection h (W/(m² K)) | 10, free convection in air (range 2-25) [Incropera] | same |

Sources:
- [CDA] Copper Development Association, alloy C10100 data sheet, copper.org. The values above are the CDA
  table as reproduced by distributors (for example Aviva Metals, avivametals.com/products/c10100-oxygen-free-electronic-copper);
  copper.org itself refused automated access when this was written. Yield strength of the annealed temper
  OS050 as listed on makeitfrom.com (Annealed OS050 C10100 Copper), a secondary source.
- [L&N] H. M. Ledbetter and E. R. Naimon, Elastic properties of metals and alloys. II. Copper,
  J. Phys. Chem. Ref. Data 3 (1974) 897-935, doi:10.1063/1.3253150, Table 10 (recommended room-temperature
  polycrystal values, averaged over annealed specimens).
- [SGL] SGL Carbon, SIGRAFINE R4550 technical data sheet, TDS R4550.00, 04 2019 ("typical average values";
  SGL asks designers to confirm values with its technical sales team).
- [JANAF] M. W. Chase, NIST-JANAF Thermochemical Tables, 4th ed., J. Phys. Chem. Ref. Data Monograph 9 (1998),
  table C-002, C1(ref), janaf.nist.gov/tables/C-002.txt.
- [Incropera] F. P. Incropera, D. P. DeWitt, T. L. Bergman, A. S. Lavine, Fundamentals of Heat and Mass
  Transfer, Table 1.1.

## Known limit: 780 K in graphite at s = 0.5 mm

With constant properties, graphite at s = 0.5 mm reaches about 780 K at the end of the pulse (the demo prints
782 K). JANAF gives Cp = 19.83 J/(mol K) at 800 K, against 8.517 at 298.15 K: c is 2.3 times higher at 780 K.
Integrating c(T) for the same deposit gives a temperature rise of 319 K (612 K) instead of 491 K. Constant c
overestimates that row by 54%. The rows with s >= 1 mm (dT <= 123 K) are less affected but not exact. alpha
and E of graphite also change with T. Temperature-dependent properties are out of scope.
