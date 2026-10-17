// Quarter square [0, L]^2 for 2D heat (P4-e, P4-h Gaussian; P4-f i slab, iii uniform plate;
// docs/learning/phase4.md sections 4 and 10). Straight edges, so T6 holds quadratics exactly.
// Coordinates in m. Mesh with: gmsh square.geo -2 -order {1|2} -clscale s -format msh41
L = 0.012;
lc = 0.1 * L;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, L, 0, lc};
Point(4) = {0, L, 0, lc};

Line(1) = {1, 2};   // y = 0
Line(2) = {2, 3};   // x = L
Line(3) = {3, 4};   // y = L
Line(4) = {4, 1};   // x = 0

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Physical Surface("domain") = {1};
