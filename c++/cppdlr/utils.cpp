#include <numbers>
#include "utils.hpp"
#include <nda/blas.hpp>

using namespace nda;
using std::numbers::pi;

static constexpr auto _ = range::all;

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

    double num = 0, den = 0, dif = 0, q = 0;

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

  // TODO: Change this to a row-pivoted code in order to make things contiguous
  // in memory

  std::tuple<nda::matrix<double>, nda::vector<double>, nda::vector<int>> pivrgs(nda::matrix<double> a, double eps) {

    // Get matrix dimensions

    auto [m, n] = a.shape();

    // Compute norms of rows of input matrix

    auto norms = nda::vector<double>(m);

    double epsscal = 0; // Scaled eps threshold parameter

    // Get norms of rows of matrix, and Frobenius norm of matrix

    for (int j = 0; j < m; ++j) {
      norms(j) = blas::dot(a(j, _), a(j, _));
      epsscal += norms(j);
    }

    epsscal *= eps * eps;

    // Begin pivoted double Gram-Schmidt procedure

    int jpiv = 0, jj = 0;
    double nrm = 0;
    auto piv   = nda::arange(0, m);
    auto tmp   = nda::vector<double>(n);

    for (int j = 0; j < m; ++j) {

      // Find next pivot

      jpiv = j;
      nrm  = norms(j);
      for (int k = j + 1; k < m; ++k) {
        if (norms(k) > nrm) {
          jpiv = k;
          nrm  = norms(k);
        }
      }

      // Swap current row with chosen pivot row

      tmp        = a(j, _);
      a(j, _)    = a(jpiv, _);
      a(jpiv, _) = tmp;

      nrm         = norms(j);
      norms(j)    = norms(jpiv);
      norms(jpiv) = nrm;

      jj        = piv(j);
      piv(j)    = piv(jpiv);
      piv(jpiv) = jj;

      // Orthogonalize current rows (now the chosen pivot row) against all
      // previously chosen rows

      for (int k = 0; k < j; ++k) { a(j, _) = a(j, _) - a(k, _) * blas::dot(a(j, _), a(k, _)); }

      // Get norm of current row

      nrm = blas::dot(a(j, _), a(j, _));

      // Terminate if sufficiently small, and return previously selected rows
      // (not including current row)

      if (nrm <= epsscal) { return {a(range(0, j), _), norms(range(0, j)), piv(range(0, j))}; };

      // Normalize current row

      a(j, _) = a(j, _) * (1 / sqrt(nrm));

      // Orthogonalize remaining rows against current row

      for (int k = j + 1; k < m; ++k) {

        if (norms(k) <= epsscal) { continue; } // Can skip rows with norm less than tolerance

        a(k, _)  = a(k, _) - a(j, _) * blas::dot(a(k, _), a(j, _));
        norms(k) = blas::dot(a(k, _), a(k, _));
      }
    }

    return {a, norms, piv};
  }

  // TODO: write test program for this
  nda::vector<double> eqptsrel(int n) {

    auto t = nda::vector<double>(n);

     for (int i = 0; i < n - 1; ++i) {
       if (i <= (n - 1) / 2) {
         t(i) = i * 1.0 / (n - 1);
       } else {
         t(i) = -(n - 1 - i) * 1.0 / (n - 1);
       }
     }
     t(n - 1) = 1;

    return t;
  }

} // namespace cppdlr
