#include "mechanics/error_norms.hpp"

#include <cmath>

#include "fem/assembly.hpp"
#include "fem/element_geometry.hpp"
#include "fem/shape_functions.hpp"

namespace tmfem {

ElasticErrors elastic_errors(const Elastic2DProblem& p, const Eigen::VectorXd& u_m, const DisplacementField& u_exact,
                             const StressField& sigma_exact, int n_gauss) {
  const Eigen::Matrix3d D = p.D();
  const Eigen::Matrix3d C = D.inverse();  // compliance
  const QuadRuleTri rule = triangle_collapsed_gauss(n_gauss);
  double e_L2 = 0, n_L2 = 0, e_E = 0, n_E = 0, area = 0;
  int n_cells = 0;
  for (const CellBlock* b : p.mesh.blocks_of(p.region)) {
    for (int e = 0; e < b->nodes.rows(); ++e) {
      const Eigen::MatrixXd x_e = p.mesh.cell_coordinates(*b, e);
      const std::vector<int> dofs = cell_dofs(*b, e, 2);
      Eigen::VectorXd u_e(static_cast<Eigen::Index>(dofs.size()));
      for (std::size_t i = 0; i < dofs.size(); ++i) u_e(static_cast<Eigen::Index>(i)) = u_m(dofs[i]);
      for (std::size_t q = 0; q < rule.w.size(); ++q) {
        const ShapeEval s = shape(b->type, rule.xi[q]);
        const CellGeometry g = cell_geometry(s, x_e);
        const double dA = g.detJ * rule.w[q];
        Eigen::Vector2d u_h = Eigen::Vector2d::Zero();
        for (Eigen::Index a = 0; a < s.N.size(); ++a) u_h += s.N(a) * u_e.segment<2>(2 * a);
        const Eigen::Vector3d s_h = D * strain_displacement(g.dN_dx) * u_e;
        const Eigen::Vector2d u = u_exact(g.x);
        const Eigen::Vector3d sig = sigma_exact(g.x);
        e_L2 += (u - u_h).squaredNorm() * dA;
        n_L2 += u.squaredNorm() * dA;
        e_E += (sig - s_h).dot(C * (sig - s_h)) * dA;
        n_E += sig.dot(C * sig) * dA;
        area += dA;
      }
      ++n_cells;
    }
  }
  return {std::sqrt(e_L2), std::sqrt(n_L2), std::sqrt(e_E), std::sqrt(n_E), std::sqrt(area / n_cells)};
}

}  // namespace tmfem
