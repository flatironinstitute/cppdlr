#pragma once

namespace cppdlr {

  /** Evaluate Lehmann kernel with relative format for time variable */

  double kfun(double t, double om);

  /** Evaluate Lehmann kernel with absolute format for time variable */

  double kfun_abs(double t, double om);

} // namespace cppdlr
