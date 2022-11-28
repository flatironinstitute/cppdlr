#include "dlr_kernels.hpp"
#include <cmath>

using namespace std;

namespace cppdlr {

  double kfun(double t, double om) {

    double val;

    if (om>=0.0) {
      val = exp(-t*om)/(1.0+exp(-om));
    } else {
      val = exp((1.0-t)*om)/(1.0+exp(om));
    }

    return val;

  }
}
