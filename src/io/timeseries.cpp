#include "io/timeseries.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

#include "io/csv_reader.hpp"

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
  auto cols = read_numeric_csv(path, 2);
  return TimeSeries(std::move(cols[0]), std::move(cols[1]));
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
