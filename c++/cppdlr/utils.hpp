#pragma once
#include <vector>
#include <nda/nda.hpp>

namespace cppdlr {

  /** 
   * Initialize subroutine barycheb for barycentric Lagrange interpolation at
   * Chebyshev nodes.
   *
   * @param   n   number of Chebyshev nodes
   * @param   xc  n Chebyshev nodes of the first kind 
   * @param   wc  barycentric interpolation weights at Chebyshev nodes of the
   * first kind
   */

  class barycheb {

    public:

      barycheb(int n);

    private:

      nda::vector<double> xc; /// Chebshev nodes
      nda::vector<double> wc; /// Chebshev weights

  };
}
