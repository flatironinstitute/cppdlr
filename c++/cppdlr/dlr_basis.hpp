#pragma once
#include <nda/nda.hpp>

namespace cppdlr {

  /**
   * DLR basis for a given accuracy and cutoff parameter
   */

  class dlr_basis {

    public:
    dlr_basis(double lambda, double eps);

    private:
    int r;                     /// DLR rank
    nda::vector<double> omega; /// DLR frequencies
  };

} // namespace cppdlr