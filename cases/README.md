# Cases

- `steel_plate/`: phase 1 steel plate, synthetic sol-air input.
- `bar/`: phase 2 uniform bar (P2-b).
- `kirsch/`: phase 2 Kirsch plate with a hole (P2-c to P2-f), quarter model.

The `.msh` files are generated from the `.geo` files by `scripts/make_meshes.sh` with
**Gmsh 4.12.1** (ASCII msh 4.1) and committed. That way CI needs no Gmsh, and the
validation numbers are reproducible. Another Gmsh version gives slightly different meshes.
