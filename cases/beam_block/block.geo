// Beam-heated block, quarter of the 40 mm x 40 mm cross-section (demo D4, docs/learning/phase4.md
// section 11). The beam runs along z through the origin; symmetry on x = 0 and y = 0.
// Element size s/10 within 3 s of the beam axis, growing to W/10 at the outer edges.
// Coordinates in m. Mesh with: gmsh block.geo -2 -order 2 -setnumber s_mm <s> -format msh41
DefineConstant[ s_mm = 1.0 ];
s = s_mm * 1e-3;
W = 0.020;

Point(1) = {0, 0, 0, s / 10};
Point(2) = {W, 0, 0, W / 10};
Point(3) = {W, W, 0, W / 10};
Point(4) = {0, W, 0, W / 10};

Line(1) = {1, 2};   // y = 0
Line(2) = {2, 3};   // x = W
Line(3) = {3, 4};   // y = W
Line(4) = {4, 1};   // x = 0

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Field[1] = MathEval;
Field[1].F = Sprintf("Min(%g + 0.1 * Max(Sqrt(x*x + y*y) - %g, 0), %g)", s / 10, 3 * s, W / 10);
Background Field = 1;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Physical Curve("sym_y") = {1};   // u_y = 0, insulated
Physical Curve("outer") = {2, 3};  // free; Robin during cooling
Physical Curve("sym_x") = {4};   // u_x = 0, insulated
Physical Point("centre") = {1};
Physical Surface("domain") = {1};
