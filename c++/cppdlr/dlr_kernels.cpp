#include "dlr_kernels.hpp"
#include <numbers>

using namespace std;
using namespace std::numbers;

namespace cppdlr {

  double k_it(double t, double om) {

    if (t >= 0) {
      return k_it_abs(t, om);
    } else {
      return k_it_abs(-t, -om);
    }
  }

  double k_it_abs(double t, double om) {

    if (om >= 0.0) {
      return exp(-t * om) / (1.0 + exp(-om));
    } else {
      return exp((1.0 - t) * om) / (1.0 + exp(om));
    }
  }

  std::complex<double> k_if(int n, double om) { return 1.0 / (n * pi * 1i + om); }

} // namespace cppdlr
