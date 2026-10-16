#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include "core/version.hpp"

// Skeleton check: the core library links, Eigen is usable, Catch2 is wired into CTest.
TEST_CASE("phase0 skeleton links tmfem_core and Eigen", "[phase0]") {
  REQUIRE_FALSE(tmfem::version().empty());

  Eigen::Matrix2d A;
  A << 4.0, 1.0,
       1.0, 3.0;
  const Eigen::Vector2d b(1.0, 2.0);
  const Eigen::Vector2d x = A.ldlt().solve(b);

  // Exact solution of the 2x2 system: x = (1/11, 7/11).
  REQUIRE((x - Eigen::Vector2d(1.0 / 11.0, 7.0 / 11.0)).norm() < 1e-15);
}
