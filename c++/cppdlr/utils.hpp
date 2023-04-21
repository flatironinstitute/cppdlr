#pragma once
#include <vector>
#include <nda/nda.hpp>
#include <nda/blas.hpp>
#include <type_traits>

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
      tmp         = aa(j, _);
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
        aa(k, _) = aa(k, _) - aa(j, _) * blas::dotc(aa(j, _), aa(k, _));
        norms(k) = real(blas::dotc(aa(k, _), aa(k, _)));
      }
    }

    return {aa, norms, piv};
  }

  nda::vector<double> eqptsrel(int n);

  template <nda::MemoryArray T> using make_real_t = decltype(make_regular(real(std::declval<T>())));

  template <nda::MemoryArray T> using make_cplx_t = decltype(make_regular(std::declval<T>() * 1i));

  /**
  * @brief Contract the last dimension of an array a with the first dimension of
  * an array b
  *
  * @param a  An array/matrix/vector or array/matrix/vector view of rank at least 2
  * @param b  An array/matrix/vector or array/matrix/vector view of rank at least 2
  *
  * @return Contraction of the inner dimensions of \p a and \p b
  */
  template <nda::MemoryArray Ta, nda::MemoryArray Tb, typename Sa = nda::get_value_t<Ta>, typename Sb = nda::get_value_t<Tb>,
            typename S = typename std::common_type<Sa, Sb>::type>
    requires(nda::is_scalar_v<Sa> and nda::is_scalar_v<Sb>)
  nda::array<S, Ta::rank + Tb::rank - 2> arraymult(Ta const &a, Tb const &b) {

    // Get ranks of input arrays
    constexpr int ra = Ta::rank;
    constexpr int rb = Tb::rank;

    // Get inner dimensions of input arrays
    int p = a.shape(ra - 1);
    if (b.shape(0) != p) throw std::runtime_error("last dim of a != first dim of b");

    // Get product of outer dimensions of input arrays
    int m = a.size() / p;
    int n = b.size() / p;

    // Reshape input arrays to 2D arrays
    auto a_reshaped = nda::reshaped_view(a, std::array<int, 2>({m, p}));
    auto b_reshaped = nda::reshaped_view(b, std::array<int, 2>({p, n}));

    // Get shape of output array
    auto c_shape = std::array<int, ra + rb - 2>();
    for (int i = 0; i < ra - 1; ++i) { c_shape[i] = a.shape(i); }
    for (int i = ra - 1; i < ra + rb - 2; ++i) { c_shape[i] = b.shape(i - ra + 2); }

    // Compute the contraction, reshape, and return
    return reshape(matmul(a_reshaped, b_reshaped), c_shape);
  }

} // namespace cppdlr
