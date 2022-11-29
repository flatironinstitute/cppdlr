#include "dlr_kernels.hpp"
#include <cmath>

using namespace std;

namespace cppdlr {

  double kfun(double t, double om) {

    if (om>=0.0) {
      return exp(-t*om)/(1.0+exp(-om));
    } else {
      return exp((1.0-t)*om)/(1.0+exp(om));
    }
  }
}
