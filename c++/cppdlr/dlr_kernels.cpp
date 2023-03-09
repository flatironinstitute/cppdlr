#include "dlr_kernels.hpp"
#include <numbers>

using namespace std;
using namespace std::numbers;

namespace cppdlr {

  /** Evaluate Lehmann kernel with relative format for time variable */

  double kfun(double t, double om) {

    if (t >= 0) {
      return kfun_abs(t, om);
    } else {
      return kfun_abs(-t, -om);
    }
  }

  /** Evaluate Lehmann kernel with absolute format for time variable */

  double kfun_abs(double t, double om) {

    if (om >= 0.0) {
      return exp(-t * om) / (1.0 + exp(-om));
    } else {
      return exp((1.0 - t) * om) / (1.0 + exp(om));
    }
  }

  /** Evaluate Lehmann kernel in imaginary frequency */

  std::complex<double> kfun_if(int n, double om) {

    return 1.0/(n*pi*1i + om);

  }

} // namespace cppdlr
