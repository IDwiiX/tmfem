# Cases

- `steel_plate/`: phase 1 steel plate, synthetic sol-air input.
- `bar/`: phase 2 uniform bar (P2-b).
- `kirsch/`: phase 2 Kirsch plate with a hole (P2-c to P2-f), quarter model.

The `.msh` files are not committed. The build generates them from the `.geo` files with
`scripts/make_meshes.sh` (`cmake/Meshes.cmake`), using **Gmsh 4.12.1** (ASCII msh 4.1). The script
refuses any other version so that the validation numbers are reproducible. Another Gmsh
version gives slightly different meshes.
