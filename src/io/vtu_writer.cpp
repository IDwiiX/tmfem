#include "io/vtu_writer.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace tmfem {

namespace {

int vtk_cell_type(CellType type) {
  switch (type) {
    case CellType::Tri3: return 5;
    case CellType::Tri6: return 22;
    default: throw std::invalid_argument("write_vtu: only Tri3 and Tri6 cells are written");
  }
}

}  // namespace

void write_vtu(const std::filesystem::path& path, const Mesh& mesh, const std::vector<const CellBlock*>& cells,
               const std::vector<VtuPointField>& point_data) {
  for (const auto& f : point_data) {
    if (f.values.rows() != mesh.n_nodes() || f.values.cols() < 1 || f.values.cols() > 3) {
      throw std::invalid_argument("write_vtu: field '" + f.name + "' must have n_nodes rows and 1-3 columns");
    }
  }
  if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("write_vtu: cannot open " + path.string());
  out << std::setprecision(12);

  Eigen::Index n_cells = 0;
  for (const CellBlock* b : cells) n_cells += b->nodes.rows();

  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "<UnstructuredGrid>\n"
      << "<Piece NumberOfPoints=\"" << mesh.n_nodes() << "\" NumberOfCells=\"" << n_cells << "\">\n";

  out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (int i = 0; i < mesh.n_nodes(); ++i) out << mesh.x_m(i, 0) << ' ' << mesh.x_m(i, 1) << " 0\n";
  out << "</DataArray>\n</Points>\n";

  out << "<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (const CellBlock* b : cells) {
    for (Eigen::Index e = 0; e < b->nodes.rows(); ++e) {
      for (Eigen::Index a = 0; a < b->nodes.cols(); ++a) out << b->nodes(e, a) << (a + 1 < b->nodes.cols() ? ' ' : '\n');
    }
  }
  out << "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  long long offset = 0;
  for (const CellBlock* b : cells) {
    for (Eigen::Index e = 0; e < b->nodes.rows(); ++e) {
      offset += b->nodes.cols();
      out << offset << '\n';
    }
  }
  out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (const CellBlock* b : cells) {
    const int t = vtk_cell_type(b->type);
    for (Eigen::Index e = 0; e < b->nodes.rows(); ++e) out << t << '\n';
  }
  out << "</DataArray>\n</Cells>\n";

  out << "<PointData>\n";
  for (const auto& f : point_data) {
    out << "<DataArray type=\"Float64\" Name=\"" << f.name << "\" NumberOfComponents=\"" << f.values.cols()
        << "\" format=\"ascii\">\n";
    for (Eigen::Index i = 0; i < f.values.rows(); ++i) {
      for (Eigen::Index c = 0; c < f.values.cols(); ++c) out << f.values(i, c) << (c + 1 < f.values.cols() ? ' ' : '\n');
    }
    out << "</DataArray>\n";
  }
  out << "</PointData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace tmfem
