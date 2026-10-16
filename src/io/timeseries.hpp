#pragma once

#include <filesystem>
#include <vector>

namespace tmfem {

// Piecewise-linear time series v(t). Evaluating outside [t_first, t_last] throws:
// there is no silent extrapolation.
class TimeSeries {
 public:
  TimeSeries(std::vector<double> t_s, std::vector<double> values);

  // Two-column CSV `t_s,value`: lines starting with '#' are comments, the first other line is a header.
  static TimeSeries read_csv(const std::filesystem::path& path);

  double operator()(double t_s) const;
  double t_first() const { return t_.front(); }
  double t_last() const { return t_.back(); }
  const std::vector<double>& times() const { return t_; }
  const std::vector<double>& values() const { return v_; }

 private:
  std::vector<double> t_;
  std::vector<double> v_;
};

}  // namespace tmfem
