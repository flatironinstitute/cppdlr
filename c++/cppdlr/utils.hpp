#pragma once
#include <vector>
#include <nda/nda.hpp>
#include <nda/blas.hpp>

using namespace nda;

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

  /** 
   * @brief Rank-revealing pivoted reorthogonalized Gram-Schmidt
   *
   * Determine the epsilon-rank of a matrix and return an orthogonal basis of
   * its epsilon-row space.
   *
   * This is a translation of the Fortran subroutine "qrdgrm" by V.  Rokhlin.
   *
   * @param a   Matrix to be orthogonalized
   * @param eps Rank cutoff tolerance
   *
   * @return Tuple of (1) matrix containing whose rows form orthogonal basis of
   * row space of @p a to @p eps tolerance, (2) vector with entry n given by the
   * squared l2 norm of the orthogonal complement of nth selected row with
   * respect to subspace spanned by first n-1 selected rows, (3) vector of
   * pivots
   */
  // Type T must be scalar-valued rank 2 array/array_view or matrix/matrix_view
  template <nda::MemoryArrayOfRank<2> T, typename S = get_value_t<T>>
    requires(nda::is_scalar_v<S>)
  std::tuple<typename T::regular_type, nda::vector<double>, nda::vector<int>> pivrgs(T const &a, double eps) {

    auto _ = nda::range::all;

    // Copy input data
    auto aa = make_regular(a);

    // Get matrix dimensions
    auto [m, n] = aa.shape();

    // Compute norms of rows of input matrix, and rescale eps tolerance
    auto norms     = nda::vector<double>(m);
    double epsscal = 0; // Scaled eps threshold parameter
    for (int j = 0; j < m; ++j) {
      norms(j) = real(blas::dotc(aa(j, _), aa(j, _)));
      // TODO: Need to choose between these; this choice is consistent w/ libdlr
      //epsscal += norms(j);
      epsscal = std::max(epsscal, norms(j));
    }
    epsscal *= eps * eps;

    // Begin pivoted double Gram-Schmidt procedure
    int jpiv = 0, jj = 0;
    double nrm = 0;
    auto piv   = nda::arange(0, m);
    auto tmp   = nda::vector<S>(n);

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
      tmp        = aa(j, _);
      aa(j, _)    = aa(jpiv, _);
      aa(jpiv, _) = tmp;

      nrm         = norms(j);
      norms(j)    = norms(jpiv);
      norms(jpiv) = nrm;

      jj        = piv(j);
      piv(j)    = piv(jpiv);
      piv(jpiv) = jj;

      // Orthogonalize current rows (now the chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j; ++k) { aa(j, _) = aa(j, _) - aa(k, _) * blas::dotc(aa(k, _), aa(j, _)); }

      // Get norm of current row
      nrm = real(blas::dotc(aa(j, _), aa(j, _)));
      //nrm = nda::norm(aa(j, _));

      // Terminate if sufficiently small, and return previously selected rows
      // (not including current row)
      if (nrm <= epsscal) { return {aa(range(0, j), _), norms(range(0, j)), piv(range(0, j))}; };

      // Normalize current row
      aa(j, _) = aa(j, _) * (1 / sqrt(nrm));

      // Orthogonalize remaining rows against current row
      for (int k = j + 1; k < m; ++k) {
        if (norms(k) <= epsscal) { continue; } // Can skip rows with norm less than tolerance
        aa(k, _)  = aa(k, _) - aa(j, _) * blas::dotc(aa(j, _), aa(k, _));
        norms(k) = real(blas::dotc(aa(k, _), aa(k, _)));
      }
    }

    return {aa, norms, piv};
  }

  nda::vector<double> eqptsrel(int n);

  template <nda::MemoryArray T>
  using make_real_t = decltype(make_regular(real(std::declval<T>())));

  template <nda::MemoryArray T>
  using make_cplx_t = decltype(make_regular(imag(std::declval<T>())));

} // namespace cppdlr
