#pragma once
#include <nda/nda.hpp>

namespace cppdlr {

  /**
   * Construct DLR basis for a given accuracy and cutoff parameter by getting DLR frequencies
   */

   nda::vector<double> dlr_freq(double lambda, double eps);

   /**
   * Construct DLR imaginary time grid
   */ 

   nda::vector<double> dlr_it(nda::vector<double> dlr_rf);

} // namespace cppdlr