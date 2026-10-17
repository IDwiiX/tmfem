// Thin solid disk, quarter model (P4-c, P4-d, P4-g, P4-f ii; docs/learning/phase4.md section 10).
// Radius b, symmetry on x = 0 and y = 0, free curved rim. Coordinates in m.
// Mesh with: gmsh disk.geo -2 -order {1|2} -clscale s -format msh41
b = 0.010;
lc = 0.2 * b;

Point(1) = {0, 0, 0, lc};
Point(2) = {b, 0, 0, lc};
Point(3) = {0, b, 0, lc};

Line(1) = {1, 2};         // y = 0
Circle(2) = {2, 1, 3};    // rim, r = b
Line(3) = {3, 1};         // x = 0

Curve Loop(1) = {1, 2, 3};  // counter-clockwise
Plane Surface(1) = {1};

Physical Curve("sym_y") = {1};   // u_y = 0
Physical Curve("rim") = {2};     // traction-free; Robin in P4-f (ii)
Physical Curve("sym_x") = {3};   // u_x = 0
Physical Point("centre") = {1};
Physical Surface("domain") = {1};
