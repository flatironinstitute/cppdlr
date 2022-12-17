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

    nda::vector<double> const &getnodes();

    double interp(double x, nda::vector_const_view<double> f);

    private:
    nda::vector<double> xc; /// Chebshev nodes
    nda::vector<double> wc; /// Chebshev weights
    const int n;            /// Chebyshev n
  };

  std::tuple<nda::matrix<double>, nda::vector<double>, nda::vector<int>> pivrgs(nda::matrix<double> a, double eps);

} // namespace cppdlr