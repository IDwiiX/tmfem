#include "io/timeseries.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace tmfem {

TimeSeries::TimeSeries(std::vector<double> t_s, std::vector<double> values)
    : t_(std::move(t_s)), v_(std::move(values)) {
  if (t_.size() < 2 || t_.size() != v_.size()) {
    throw std::invalid_argument("TimeSeries: need at least two (t, v) pairs of equal length");
  }
  for (std::size_t i = 1; i < t_.size(); ++i) {
    if (!(t_[i] > t_[i - 1])) throw std::invalid_argument("TimeSeries: times must increase strictly");
  }
}

TimeSeries TimeSeries::read_csv(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("TimeSeries: cannot open " + path.string());
  std::vector<double> t;
  std::vector<double> v;
  std::string line;
  bool header_seen = false;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (!header_seen) {
      header_seen = true;
      continue;
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    double ti = 0.0;
    double vi = 0.0;
    if (!(row >> ti >> vi)) throw std::runtime_error("TimeSeries: bad row '" + line + "' in " + path.string());
    t.push_back(ti);
    v.push_back(vi);
  }
  return TimeSeries(std::move(t), std::move(v));
}

double TimeSeries::operator()(double t_s) const {
  const double slack = 1e-9 * (t_.back() - t_.front());
  if (t_s < t_.front() - slack || t_s > t_.back() + slack) {
    throw std::out_of_range("TimeSeries: t = " + std::to_string(t_s) + " s outside the data range");
  }
  const double t = std::clamp(t_s, t_.front(), t_.back());
  auto upper = std::upper_bound(t_.begin(), t_.end(), t);
  if (upper == t_.end()) return v_.back();
  const auto i = static_cast<std::size_t>(upper - t_.begin()) - 1;
  const double s = (t - t_[i]) / (t_[i + 1] - t_[i]);
  return (1.0 - s) * v_[i] + s * v_[i + 1];
}

}  // namespace tmfem
