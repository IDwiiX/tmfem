#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <map>
#include <string>

#include "mechanics/elastic2d.hpp"
#include "mechanics/stress_recovery.hpp"
#include "phase2/kirsch_common.hpp"
#include "support/validation_record.hpp"

// P2-a patch test (V2.1): a linear displacement field is reproduced on a distorted mesh, T3 and T6.
// Reference: docs/REFERENCE_VALUES.md row P2-a, tolerance 1e-10 relative.
// Mesh: the MacNeal-Harder patch (0.24 x 0.12 m rectangle, interior nodes (0.04, 0.02), (0.18, 0.03),
// (0.16, 0.08), (0.08, 0.08)), each quadrilateral split into two triangles. In the T6 version every
// interior edge is curved: its mid-side node is moved off the chord by 5% of the edge length.
// docs/learning/phase2.md section 5 explains why an isoparametric element must still pass.

using namespace tmfem;

namespace {

constexpr double kTolPatch = 1e-10;  // P2-a: 1e-10 relative
constexpr double kBow = 0.05;        // interior mid-side offset / edge length (T6)

// Exact linear field u = c + A x and its constant strain [eps_xx, eps_yy, gamma_xy].
Eigen::Vector2d u_exact(const Eigen::Vector2d& x) {
  return {1e-3 + 2e-3 * x(0) + 1e-3 * x(1), -1e-3 + 1e-3 * x(0) + 3e-3 * x(1)};
}
const Eigen::Vector3d kStrain(2e-3, 3e-3, 2e-3);

enum Tag { kDomain = 1, kBottom, kRight, kTop, kLeft, kP0, kP1 };

Mesh patch_mesh(bool quadratic) {
  Mesh m;
  std::vector<Eigen::Vector2d> x{{0.0, 0.0}, {0.24, 0.0}, {0.24, 0.12}, {0.0, 0.12},
                                 {0.04, 0.02}, {0.18, 0.03}, {0.16, 0.08}, {0.08, 0.08}};
  const std::vector<std::array<int, 4>> quads{{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
  std::vector<std::array<int, 3>> tris;
  for (const auto& q : quads) {
    tris.push_back({q[0], q[1], q[2]});
    tris.push_back({q[0], q[2], q[3]});
  }
  const std::vector<std::array<int, 2>> sides{{0, 1}, {1, 2}, {2, 3}, {3, 0}};  // bottom, right, top, left
  const auto is_boundary = [](int i, int j) { return i < 4 && j < 4; };

  std::map<std::pair<int, int>, int> mid;  // edge (lo, hi) -> mid-side node
  const auto mid_node = [&](int i, int j) {
    const auto key = std::minmax(i, j);
    const auto it = mid.find(key);
    if (it != mid.end()) return it->second;
    const Eigen::Vector2d d = x[static_cast<std::size_t>(key.second)] - x[static_cast<std::size_t>(key.first)];
    Eigen::Vector2d p = 0.5 * (x[static_cast<std::size_t>(i)] + x[static_cast<std::size_t>(j)]);
    if (!is_boundary(i, j)) p += kBow * Eigen::Vector2d(-d(1), d(0));  // fixed side for each edge
    x.push_back(p);
    mid[key] = static_cast<int>(x.size()) - 1;
    return mid[key];
  };

  const int npc = quadratic ? 6 : 3;
  ConnMatrix tri_conn(static_cast<Eigen::Index>(tris.size()), npc);
  for (std::size_t e = 0; e < tris.size(); ++e) {
    const auto& t = tris[e];
    const auto r = static_cast<Eigen::Index>(e);
    tri_conn(r, 0) = t[0];
    tri_conn(r, 1) = t[1];
    tri_conn(r, 2) = t[2];
    if (quadratic) {
      tri_conn(r, 3) = mid_node(t[0], t[1]);
      tri_conn(r, 4) = mid_node(t[1], t[2]);
      tri_conn(r, 5) = mid_node(t[2], t[0]);
    }
  }
  m.blocks.push_back({quadratic ? CellType::Tri6 : CellType::Tri3, kDomain, tri_conn});
  const char* side_names[4] = {"bottom", "right", "top", "left"};
  for (int s = 0; s < 4; ++s) {
    const auto& sd = sides[static_cast<std::size_t>(s)];
    ConnMatrix e(1, quadratic ? 3 : 2);
    e(0, 0) = sd[0];
    e(0, 1) = sd[1];
    if (quadratic) e(0, 2) = mid_node(sd[0], sd[1]);
    m.blocks.push_back({quadratic ? CellType::Line3 : CellType::Line2, kBottom + s, e});
    m.physical[side_names[s]] = {1, kBottom + s};
  }
  m.blocks.push_back({CellType::Point1, kP0, ConnMatrix::Constant(1, 1, 0)});
  m.blocks.push_back({CellType::Point1, kP1, ConnMatrix::Constant(1, 1, 1)});
  m.physical["domain"] = {2, kDomain};
  m.physical["p0"] = {0, kP0};
  m.physical["p1"] = {0, kP1};
  m.x_m.resize(static_cast<Eigen::Index>(x.size()), 2);
  for (std::size_t i = 0; i < x.size(); ++i) m.x_m.row(static_cast<Eigen::Index>(i)) = x[i].transpose();
  return m;
}

DisplacementBC exact_bc(const std::string& physical, int comp) {
  return {physical, comp, [comp](const Eigen::Vector2d& x) { return u_exact(x)(comp); }};
}

struct PatchErrors {
  double u_rel, sigma_rel;
};

PatchErrors patch_errors(const Elastic2DProblem& p) {
  const Elastic2DSolution sol = solve_elastic2d(p);
  double du = 0, umax = 0;
  for (int i = 0; i < p.mesh.n_nodes(); ++i) {
    const Eigen::Vector2d ue = u_exact(p.mesh.x_m.row(i).transpose());
    du = std::max(du, (sol.u_m.segment<2>(2 * i) - ue).cwiseAbs().maxCoeff());
    umax = std::max(umax, ue.cwiseAbs().maxCoeff());
  }
  const Eigen::Vector3d s_exact = p.D() * kStrain;
  double ds = 0;
  for (const auto& s : gauss_point_stresses(p, sol.u_m)) ds = std::max(ds, (s - s_exact).cwiseAbs().maxCoeff());
  return {du / umax, ds / s_exact.cwiseAbs().maxCoeff()};
}

// (a) Dirichlet version: exact u on every boundary node; the interior nodes are free.
Elastic2DProblem dirichlet_patch(bool quadratic, PlaneCondition plane) {
  Elastic2DProblem p{patch_mesh(quadratic), test::phase2_material(), plane, "domain", {}, {}};
  for (const char* s : {"bottom", "right", "top", "left"}) {
    p.displacements.push_back(exact_bc(s, 0));
    p.displacements.push_back(exact_bc(s, 1));
  }
  return p;
}

// (b) Traction version: exact tractions sigma n on all four sides, minimal supports (p0 in x and y, p1 in y).
Elastic2DProblem traction_patch(bool quadratic, PlaneCondition plane) {
  Elastic2DProblem p{patch_mesh(quadratic), test::phase2_material(), plane, "domain", {}, {}};
  const Eigen::Vector3d s = p.D() * kStrain;
  const Eigen::Matrix2d sig = (Eigen::Matrix2d() << s(0), s(2), s(2), s(1)).finished();
  const std::pair<const char*, Eigen::Vector2d> normals[4] = {
      {"bottom", {0, -1}}, {"right", {1, 0}}, {"top", {0, 1}}, {"left", {-1, 0}}};
  for (const auto& [name, n] : normals) {
    const Eigen::Vector2d t = sig * n;
    p.tractions.push_back({name, [t](const Eigen::Vector2d&) { return t; }});
  }
  p.displacements = {exact_bc("p0", 0), exact_bc("p0", 1), exact_bc("p1", 1)};
  return p;
}

void run_patch(bool quadratic, const std::string& test_name) {
  double worst = 0;
  for (const auto plane : {PlaneCondition::PlaneStress, PlaneCondition::PlaneStrain}) {
    const char* pname = plane == PlaneCondition::PlaneStress ? "plane stress" : "plane strain";
    const PatchErrors d = patch_errors(dirichlet_patch(quadratic, plane));
    const PatchErrors t = patch_errors(traction_patch(quadratic, plane));
    std::printf("P2-a %s %s: Dirichlet u %.2e sigma %.2e | traction u %.2e sigma %.2e\n", quadratic ? "T6" : "T3",
                pname, d.u_rel, d.sigma_rel, t.u_rel, t.sigma_rel);
    INFO(pname);
    CHECK(d.u_rel < kTolPatch);
    CHECK(d.sigma_rel < kTolPatch);
    CHECK(t.u_rel < kTolPatch);
    CHECK(t.sigma_rel < kTolPatch);
    worst = std::max({worst, d.u_rel, d.sigma_rel, t.u_rel, t.sigma_rel});
  }
  test::record("P2-a", std::string(quadratic ? "T6" : "T3") + " max relative error (u, sigma)", 0.0, worst,
               "1e-10 relative", test_name);
}

}  // namespace

TEST_CASE("P2-a patch test T3 distorted mesh", "[phase2][P2-a]") {
  run_patch(false, "P2-a patch test T3 distorted mesh");
}

TEST_CASE("P2-a patch test T6 distorted mesh", "[phase2][P2-a]") {
  run_patch(true, "P2-a patch test T6 distorted mesh");
}
