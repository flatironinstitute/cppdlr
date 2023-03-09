#pragma once
#include <complex>

namespace cppdlr {

  /** Evaluate Lehmann kernel in imaginary time with relative format for time variable */

  double kfun(double t, double om);

  /** Evaluate Lehmann kernel in imaginary time with absolute format for time variable */

  double kfun_abs(double t, double om);

  /** Evaluate Lehmann kernel in imaginary frequency */

  std::complex<double> kfun_if(int n, double om);

} // namespace cppdlr
