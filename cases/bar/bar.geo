// Uniform bar in tension (P2-b): 100 mm x 20 mm, unstructured triangles. Coordinates in m.
// Mesh with: gmsh bar.geo -2 -order {1|2} -format msh41
L = 0.100;
H = 0.020;
lc = 0.004;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {0, H, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Curve("left") = {4};      // u_x = 0
Physical Curve("right") = {2};     // traction (sigma, 0)
Physical Point("corner") = {1};    // u_y = 0
Physical Surface("domain") = {1};

Mesh.Algorithm = 6;
