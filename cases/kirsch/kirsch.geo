// Kirsch plate with a hole, quarter model using symmetry (docs/learning/phase2.md section 9).
// Full plate: square of side 10a with a central hole of radius a; this is the quarter [0, 5a]^2.
// Remote stress along y. Exact Kirsch tractions go on outer_right and outer_top.
// Coordinates in m. Mesh with: gmsh kirsch.geo -2 -order {1|2} -clscale s -format msh41
a = 0.010;        // hole radius, m (Kirsch stresses are scale-free)
W = 5 * a;        // half side of the full plate
lc_hole = 0.2 * a;  // element size at the hole
lc_out = 1.0 * a;   // element size at the outer corners

Point(1) = {0, 0, 0, lc_hole};  // arc centre (not meshed)
Point(2) = {a, 0, 0, lc_hole};
Point(3) = {W, 0, 0, lc_out};
Point(4) = {W, W, 0, lc_out};
Point(5) = {0, W, 0, lc_out};
Point(6) = {0, a, 0, lc_hole};

Line(1) = {2, 3};         // y = 0
Line(2) = {3, 4};         // x = W
Line(3) = {4, 5};         // y = W
Line(4) = {5, 6};         // x = 0
Circle(5) = {6, 1, 2};    // hole, r = a

Curve Loop(1) = {1, 2, 3, 4, 5};  // counter-clockwise
Plane Surface(1) = {1};

Physical Curve("sym_y") = {1};        // u_y = 0
Physical Curve("outer_right") = {2};  // traction sigma_exact . (1, 0)
Physical Curve("outer_top") = {3};    // traction sigma_exact . (0, 1)
Physical Curve("sym_x") = {4};        // u_x = 0
Physical Curve("hole") = {5};         // traction-free
Physical Surface("domain") = {1};

Mesh.Algorithm = 6;  // Frontal-Delaunay
