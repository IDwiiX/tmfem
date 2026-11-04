# Cases

- `steel_plate/`: phase 1 steel plate, synthetic sol-air input.
- `bar/`: phase 2 uniform bar (P2-b).
- `kirsch/`: phase 2 Kirsch plate with a hole (P2-c to P2-f), quarter model.
- `disk/`: phase 4 thin disk with a parabolic temperature (P4-c, P4-d, P4-g, P4-f ii).
- `heat2d/`: phase 4 quarter square for 2D heat (P4-e, P4-h, P4-f i and iii).
- `beam_block/`: phase 4 beam-heated block, demo D4 (not a validation; see its README).

The `.msh` files are not committed. The build generates them from the `.geo` files with
`scripts/make_meshes.sh` (`cmake/Meshes.cmake`), using **Gmsh 4.12.1** (ASCII msh 4.1). The script
refuses any other version so that the validation numbers are reproducible. Another Gmsh
version gives slightly different meshes.
