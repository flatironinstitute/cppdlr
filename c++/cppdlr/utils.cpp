#include <numbers>
#include "utils.hpp"
#include <nda/blas.hpp>

using namespace nda;
using std::numbers::pi;

auto _ = range::all;

namespace cppdlr {

  barycheb::barycheb(int n) : xc(n), wc(n), n(n) {

    //auto ns = arange(1,n+1); // This isn't implemented yet

    // vector<int> ns(n);
    // for (int i = 0; i < n; ++i) { ns(i) = i; }

    // auto c = (2 * ns + 1) * pi / (2 * n);
    // xc     = cos(c(range(n-1,-1,-1))); // This line isn't working
    // wc     = sin(c);
    // for (int i = 0; i < n; ++i) {wc(i) = pow(-1,i+1)*wc(i); };

    for (int i = 0; i < n; i++) {
      auto c        = (2 * i + 1) * pi / (2 * n);
      xc(n - i - 1) = std::cos(c);
      wc(n - i - 1) = (1 - 2 * (i % 2)) * std::sin(c);
    }
  }

  nda::vector<double> const &barycheb::getnodes() { return xc; }

  double barycheb::interp(double x, nda::vector_const_view<double> f) {

    for (int i = 0; i < n; ++i) {
      if (x == xc(i)) { return f(i); }
    }

    double num = 0.0, den = 0.0;
    double dif, q;

    for (int i = 0; i < n; ++i) {
      dif = x - xc(i);
      q   = wc(i) / dif;
      num = num + q * f(i);
      den = den + q;
    }

    return num / den;
  }

  /** 
   * Rank-revealing pivoted reorthogonalized Gram-Schmidt
   *
   * This is a translation of the fortran subroutine "qrdgrm" by V. Rokhlin.
   */

  std::tuple<nda::matrix<double>, nda::vector<double>, nda::vector<int>> pivrgs(nda::matrix<double> a, double eps) {

    // Get matrix dimensions

    auto [m, n] = a.shape();

    // Compute norms of columns of input matrix

    nda::vector<double> norms(n);

    double epsscal = 0.0; // Scaled eps threshold parameter

    // Get norms of columns of matrix, and Frobenius norm of matrix

    for (int j = 0; j < n; ++j) {
      norms(j) = blas::dot(a(_, j), a(_, j));
      epsscal += norms(j);
    }

    epsscal *= eps * eps;

    // Begin pivoted double Gram-Schmidt procedure

    int jpiv, jj;
    double nrm;
    auto piv = nda::arange(0, n);
    nda::vector<double> tmp(m);

    for (int j = 0; j < n; ++j) {

      // Find next pivot

      // TODO: make this more concise. Just need index of max of norms(j:end).
      jpiv = j;
      nrm  = norms(j);
      for (int k = j + 1; k < n; ++k) {
        if (norms(k) > nrm) {
          jpiv = k;
          nrm  = norms(k);
        }
      }

      // Swap current column with chosen pivot column

      tmp        = a(_, j);
      a(_, j)    = a(_, jpiv);
      a(_, jpiv) = tmp;

      nrm         = norms(j);
      norms(j)    = norms(jpiv);
      norms(jpiv) = nrm;

      jj        = piv(j);
      piv(j)    = piv(jpiv);
      piv(jpiv) = jj;

      // Orthogonalize current column (now the chosen pivot column) against all
      // previously chosen columns

      for (int k = 0; k < j; ++k) { a(_, j) = a(_, j) - a(_, k) * blas::dot(a(_, j), a(_, k)); }

      // Get norm of current column

      nrm = blas::dot(a(_, j), a(_, j));

      // Terminate if sufficiently small, and return previously selected columns
      // (not including current column)

      if (nrm <= epsscal) { return {a(_, range(0, j)), norms(range(0, j)), piv(range(0, j))}; };

      // Normalize current column

      a(_, j) = a(_, j) * (1 / sqrt(nrm));

      // Orthogonalize remaining columns against current column

      for (int k = j + 1; k < n; ++k) {

        if (norms(k) <= epsscal) { continue; } // Can skip columns with norm less than tolerance

        // [Q] Below could be made more efficient using for loop; is this advisable?

        a(_, k)  = a(_, k) - a(_, j) * blas::dot(a(_, k), a(_, j));
        norms(k) = blas::dot(a(_, k), a(_, k));
      }
    }

    return {a, norms, piv};
  }

} // namespace cppdlr