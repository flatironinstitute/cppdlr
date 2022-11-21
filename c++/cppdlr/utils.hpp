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

    // [Q] Should this be const or ref? How to do this?
    nda::vector_view<double> getnodes();

    // [Q] Is &x and &f correct?
    // [Q] Should f be a view?
    double interp(double &x, nda::vector<double> &f);

    private:
    nda::vector<double> xc; /// Chebshev nodes
    nda::vector<double> wc; /// Chebshev weights
  };
} // namespace cppdlr
