# tmfem

A from-scratch C++17 finite-element solver for 2D transient heat conduction, linear elasticity,
von Mises elastoplasticity and one-way thermo-mechanical coupling. Every capability is checked by a
validation test against an analytical solution or a reference result. It is a learning and validation
project; it does not replace ANSYS or Abaqus.

**Status:** skeleton only. No solver capability is validated yet.

## Results

Every number below is written by the test that asserts it (`results/validation/`).

| Benchmark | Quantity | Reference | tmfem | Error | Test |
| --- | --- | --- | --- | --- | --- |
| | | | | | |

## Build and test

Requires a C++17 compiler, CMake >= 3.16 and Eigen 3 (`sudo apt install libeigen3-dev`).
Catch2 is downloaded at a pinned version during configuration.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Run one phase with `ctest --test-dir build -L phase1`, or one reference row with
`ctest --test-dir build -L P2-c`.

## Limitations

To be completed as phases are validated. Out of scope by design: 3D, large deformation, dynamics,
temperature-dependent properties, contact, radiation, real particle-shower energy maps.

## Credits

Benchmark cases adapted from M. R. Kabiri, ENSAM Meknès lab materials (with permission);
analytical solutions from the textbooks cited next to each test.
