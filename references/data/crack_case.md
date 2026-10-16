# Cracked plate (stretch goal S1)

Source: M. R. Kabiri, ENSAM Meknes, TP 3, pp. 9-12 (Abaqus).

Known:
- Center crack of length 2a in a strip of width 2b under remote tension sigma (the source
  highlights the finite-width formula K_I = sigma sqrt(pi a) / sqrt(cos(pi a / (2b)))).
- Source reports K_I analytic 840 MPa sqrt(m) and numerical 922.25 MPa sqrt(m) (about 10% apart),
  from a single-point extraction of sigma_yy near the tip.
- Numerical stresses (MPa) at points A (tip), B (crack face), C (ahead):
  A: sxx 203, syy 474.6, sxy 0; B: sxx -202, syy -12.54, sxy -1.12; C: sxx -16.56, syy 199.52, sxy 1.03.

Issues:
- The near-tip formula is printed as sigma_yy = K_I sqrt(2 pi r) f(theta); it must be
  K_I / sqrt(2 pi r) f(theta).
- The analytic K_I is written as sigma sqrt(2 pi a); for a crack of total length 2a the
  infinite-plate value is sigma sqrt(pi a). Check which length the source calls a.
- B lies on the traction-free crack face, so analytically sigma_yy = 0 there (the source table
  lists infinity).
- Geometry (a, b, sigma, units) is missing: ask before reproducing the 840 / 922 numbers.

Independent of the missing data, S1 can be validated against the Feddersen formula for any
chosen a/b.
