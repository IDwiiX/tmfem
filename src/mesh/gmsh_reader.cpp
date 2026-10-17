#include "mesh/gmsh_reader.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace tmfem {

namespace {

using EntityKey = std::pair<int, int>;  // (dim, entity tag)

struct EntityKeyHash {
  std::size_t operator()(const EntityKey& k) const {
    return std::hash<long long>()((static_cast<long long>(k.first) << 32) ^ static_cast<unsigned>(k.second));
  }
};

using EntityPhysicals = std::unordered_map<EntityKey, std::vector<int>, EntityKeyHash>;

[[noreturn]] void fail(const std::string& what) { throw std::runtime_error("read_gmsh41: " + what); }

void expect(std::istream& in, const std::string& token) {
  std::string s;
  if (!(in >> s) || s != token) fail("expected " + token + ", found '" + s + "'");
}

void read_format(std::istream& in) {
  double version = 0;
  int file_type = 0, data_size = 0;
  in >> version >> file_type >> data_size;
  if (std::abs(version - 4.1) > 1e-12 || file_type != 0) fail("only ASCII msh 4.1 is supported");
  expect(in, "$EndMeshFormat");
}

void read_physical_names(std::istream& in, Mesh& mesh) {
  int n = 0;
  in >> n;
  for (int i = 0; i < n; ++i) {
    int dim = 0, tag = 0;
    in >> dim >> tag;
    std::string rest;
    std::getline(in, rest);
    const auto q0 = rest.find('"');
    const auto q1 = rest.rfind('"');
    if (q0 == std::string::npos || q1 <= q0) fail("bad physical name");
    mesh.physical[rest.substr(q0 + 1, q1 - q0 - 1)] = {dim, tag};
  }
  expect(in, "$EndPhysicalNames");
}

std::vector<int> read_tags(std::istream& in) {
  std::size_t n = 0;
  in >> n;
  std::vector<int> tags(n);
  for (auto& t : tags) in >> t;
  return tags;
}

void read_entities(std::istream& in, EntityPhysicals& physicals) {
  std::array<int, 4> count{};
  for (auto& c : count) in >> c;
  for (int dim = 0; dim < 4; ++dim) {
    for (int i = 0; i < count[static_cast<std::size_t>(dim)]; ++i) {
      int tag = 0;
      in >> tag;
      double box = 0;
      const int n_box = dim == 0 ? 3 : 6;
      for (int k = 0; k < n_box; ++k) in >> box;
      physicals[{dim, tag}] = read_tags(in);
      if (dim > 0) read_tags(in);  // bounding entities
    }
  }
  expect(in, "$EndEntities");
}

void read_nodes(std::istream& in, std::vector<std::size_t>& tags, std::vector<std::array<double, 3>>& xyz) {
  std::size_t n_blocks = 0, n_nodes = 0, min_tag = 0, max_tag = 0;
  in >> n_blocks >> n_nodes >> min_tag >> max_tag;
  tags.reserve(n_nodes);
  xyz.reserve(n_nodes);
  for (std::size_t b = 0; b < n_blocks; ++b) {
    int dim = 0, entity = 0, parametric = 0;
    std::size_t n = 0;
    in >> dim >> entity >> parametric >> n;
    if (parametric != 0) fail("parametric node coordinates are not supported");
    const std::size_t first = tags.size();
    for (std::size_t i = 0; i < n; ++i) {
      std::size_t t = 0;
      in >> t;
      tags.push_back(t);
    }
    for (std::size_t i = 0; i < n; ++i) {
      std::array<double, 3> p{};
      in >> p[0] >> p[1] >> p[2];
      xyz.push_back(p);
    }
    if (tags.size() != first + n) fail("node block size mismatch");
  }
  expect(in, "$EndNodes");
}

CellType cell_type_from_gmsh(int gmsh_type) {
  switch (gmsh_type) {
    case 15: return CellType::Point1;
    case 1: return CellType::Line2;
    case 8: return CellType::Line3;
    case 2: return CellType::Tri3;
    case 9: return CellType::Tri6;
    default: fail("unsupported element type " + std::to_string(gmsh_type));
  }
}

// Signed area of the corner triangle; positive when counter-clockwise.
double signed_area(const Mesh& mesh, const ConnMatrix& c, Eigen::Index row) {
  const Eigen::RowVector2d p0 = mesh.x_m.row(c(row, 0));
  const Eigen::RowVector2d p1 = mesh.x_m.row(c(row, 1));
  const Eigen::RowVector2d p2 = mesh.x_m.row(c(row, 2));
  return 0.5 * ((p1(0) - p0(0)) * (p2(1) - p0(1)) - (p2(0) - p0(0)) * (p1(1) - p0(1)));
}

// Reorders clockwise triangles: swap corners 1 and 2 and, for Tri6, mid-sides 3 (edge 01) and 5 (edge 20).
void orient_counter_clockwise(const Mesh& mesh, CellBlock& block) {
  if (cell_dim(block.type) != 2) return;
  for (Eigen::Index e = 0; e < block.nodes.rows(); ++e) {
    if (signed_area(mesh, block.nodes, e) > 0.0) continue;
    std::swap(block.nodes(e, 1), block.nodes(e, 2));
    if (block.type == CellType::Tri6) std::swap(block.nodes(e, 3), block.nodes(e, 5));
  }
}

void read_elements(std::istream& in, const EntityPhysicals& physicals,
                   const std::unordered_map<std::size_t, int>& node_index, Mesh& mesh) {
  std::size_t n_blocks = 0, n_elements = 0, min_tag = 0, max_tag = 0;
  in >> n_blocks >> n_elements >> min_tag >> max_tag;
  for (std::size_t b = 0; b < n_blocks; ++b) {
    int dim = 0, entity = 0, gmsh_type = 0;
    std::size_t n = 0;
    in >> dim >> entity >> gmsh_type >> n;
    const CellType type = cell_type_from_gmsh(gmsh_type);
    const int npc = nodes_per_cell(type);
    ConnMatrix conn(static_cast<Eigen::Index>(n), npc);
    for (std::size_t e = 0; e < n; ++e) {
      std::size_t tag = 0;
      in >> tag;
      for (int a = 0; a < npc; ++a) {
        std::size_t node_tag = 0;
        in >> node_tag;
        const auto it = node_index.find(node_tag);
        if (it == node_index.end()) fail("element references unknown node " + std::to_string(node_tag));
        conn(static_cast<Eigen::Index>(e), a) = it->second;
      }
    }
    const auto ph = physicals.find({dim, entity});
    if (ph == physicals.end()) continue;
    for (const int tag : ph->second) {
      CellBlock block{type, tag, conn};
      orient_counter_clockwise(mesh, block);
      mesh.blocks.push_back(std::move(block));
    }
  }
  expect(in, "$EndElements");
}

void skip_section(std::istream& in, const std::string& name) {
  const std::string end = "$End" + name.substr(1);
  std::string s;
  while (in >> s) {
    if (s == end) return;
  }
  fail("unterminated section " + name);
}

}  // namespace

Mesh read_gmsh41(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) fail("cannot open " + path.string());
  Mesh mesh;
  EntityPhysicals physicals;
  std::vector<std::size_t> tags;
  std::vector<std::array<double, 3>> xyz;
  bool have_nodes = false;
  std::string section;
  while (in >> section) {
    if (section == "$MeshFormat") {
      read_format(in);
    } else if (section == "$PhysicalNames") {
      read_physical_names(in, mesh);
    } else if (section == "$Entities") {
      read_entities(in, physicals);
    } else if (section == "$Nodes") {
      read_nodes(in, tags, xyz);
      // Contiguous indices in increasing tag order.
      std::vector<std::size_t> order(tags.size());
      for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
      std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return tags[a] < tags[b]; });
      mesh.x_m.resize(static_cast<Eigen::Index>(tags.size()), 2);
      for (std::size_t k = 0; k < order.size(); ++k) {
        const auto& p = xyz[order[k]];
        if (p[2] != 0.0) fail("node with z != 0 in a 2D mesh");
        mesh.x_m(static_cast<Eigen::Index>(k), 0) = p[0];
        mesh.x_m(static_cast<Eigen::Index>(k), 1) = p[1];
      }
      std::vector<std::size_t> sorted_tags(tags.size());
      for (std::size_t k = 0; k < order.size(); ++k) sorted_tags[k] = tags[order[k]];
      tags = std::move(sorted_tags);
      have_nodes = true;
    } else if (section == "$Elements") {
      if (!have_nodes) fail("$Elements before $Nodes");
      std::unordered_map<std::size_t, int> node_index;
      for (std::size_t k = 0; k < tags.size(); ++k) node_index[tags[k]] = static_cast<int>(k);
      read_elements(in, physicals, node_index, mesh);
    } else if (!section.empty() && section[0] == '$') {
      skip_section(in, section);
    } else {
      fail("unexpected token '" + section + "'");
    }
  }
  if (!have_nodes) fail("no $Nodes section in " + path.string());
  return mesh;
}

}  // namespace tmfem
