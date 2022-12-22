#include "dlr_kernels.hpp"
#include <cmath>

using namespace std;

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

} // namespace cppdlr