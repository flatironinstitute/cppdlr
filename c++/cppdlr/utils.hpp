// Copyright (c) 2022-2023 Simons Foundation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0.txt
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors: Jason Kaye

#pragma once
#include "nda/concepts.hpp"
#include <nda/nda.hpp>
#include <nda/blas.hpp>


namespace cppdlr {
  using dcomplex = std::complex<double>;

  /**
   * Calculate the squared norm of a vector
   *
   * @param v The input vector
   * @return x The squared norm of the vector
   */
  double normsq(nda::MemoryVector auto const &v) { return nda::real(nda::blas::dotc(v, v)); }

  /** 
   * Class constructor for barycheb: barycentric Lagrange interpolation at
   * Chebyshev nodes.
   *
   * @param   n   number of Chebyshev nodes
   * @param   x   n Chebyshev nodes of the first kind 
   * @param   w   barycentric interpolation weights at Chebyshev nodes of the
   * first kind
   */

  class barycheb {

    public:
    barycheb(int n);

    nda::vector<double> const &getnodes();

    double interp(double xeval, nda::vector_const_view<double> f);

    private:
    nda::vector<double> x; /// Chebshev nodes
    nda::vector<double> w; /// Chebshev weights
  };

  /** 
   * Class constructor for baryleg: barycentric Lagrange interpolation at
   * Legendre nodes.
   *
   * @param   n   number of Legendre nodes
   * @param   x  n Legendre nodes 
   * @param   w  barycentric interpolation weights at Legendre nodes
   */

  class baryleg {

    public:
    baryleg(int n);

    nda::vector<double> const &getnodes();

    double interp(double xeval, nda::vector_const_view<double> f);

    private:
    nda::vector<double> x; /// Legendre nodes
    nda::vector<double> w; /// Legendre barycentric weights
  };

  /**
  * @brief Barycentric Lagrange interpolation
  *
  * @param[in] x       Interpolation nodes
  * @param[in] w       Barycentric interpolation weights
  * @param[in] f       Function values at interpolation nodes
  * @param[in] xeval   Point at which to evaluate interpolated function
  *
  * @return Interpolated function value at xeval
  */
  double baryinterp(nda::vector_const_view<double> x, nda::vector_const_view<double> w, nda::vector_const_view<double> f, double xeval);

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
   * @return Tuple of (1) matrix whose rows form orthogonal basis of
   * row space of @p a to @p eps tolerance, (2) vector with entry n given by the
   * squared l2 norm of the orthogonal complement of nth selected row with
   * respect to subspace spanned by first n-1 selected rows, (3) vector of
   * pivots
   */

  // Type T must be scalar-valued rank 2 array/array_view or matrix/matrix_view
  template <nda::MemoryArrayOfRank<2> T, nda::Scalar S = nda::get_value_t<T>>
  std::tuple<typename T::regular_type, nda::vector<double>, nda::vector<int>> pivrgs(T const &a, double eps) {

    auto _ = nda::range::all;

    // Copy input data
    auto aa = make_regular(a);

    // Get matrix dimensions
    auto [m, n] = aa.shape();
    int maxrnk  = std::min(m, n);

    // Compute norms of rows of input matrix, and rescale eps tolerance
    auto norms   = nda::vector<double>(m);
    double epssq = eps * eps;
    for (int j = 0; j < m; ++j) { norms(j) = normsq(aa(j, _)); }

    // Begin pivoted double Gram-Schmidt procedure
    int jpiv   = 0;
    double nrm = 0;
    auto piv   = nda::arange(m);
    auto tmp   = nda::vector<S>(n);

    for (int j = 0; j < maxrnk; ++j) {

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
      deep_swap(aa(j, _), aa(jpiv, _));
      std::swap(norms(j), norms(jpiv));
      std::swap(piv(j), piv(jpiv));

      // Orthogonalize current rows (now the chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j; ++k) { aa(j, _) = aa(j, _) - aa(k, _) * nda::blas::dotc(aa(k, _), aa(j, _)); }

      // Get norm of current row
      nrm = normsq(aa(j, _));

      // Terminate if sufficiently small, and return previously selected rows
      // (not including current row)
      if (nrm <= epssq) { return {aa(nda::range(0, j), _), norms(nda::range(0, j)), piv(nda::range(0, j))}; };

      // Normalize current row
      aa(j, _) /= sqrt(nrm);

      // Orthogonalize remaining rows against current row
      for (int k = j + 1; k < m; ++k) {
        if (norms(k) <= epssq) { continue; } // Can skip rows with norm less than tolerance
        aa(k, _) = aa(k, _) - aa(j, _) * nda::blas::dotc(aa(j, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }
    }

    return {aa(nda::range(maxrnk), _), norms(nda::range(maxrnk)), piv(nda::range(maxrnk))};
  }

  /** 
   * @brief Symmetrized rank-revealing pivoted reorthogonalized Gram-Schmidt
   *
   * Determine the epsilon-rank of a matrix and return an orthogonal basis of
   * its epsilon-row space, enforcing symmetrization of pivots.
   *
   * This is a translation of the Fortran subroutine "qrdgrm" by V.  Rokhlin,
   * with symmetrization added.
   *
   * @param a   Matrix to be orthogonalized
   * @param eps Rank cutoff tolerance
   *
   * @return Tuple of (1) matrix whose rows form orthogonal basis of
   * row space of @p a to @p eps tolerance, (2) vector with entry n given by the
   * squared l2 norm of the orthogonal complement of nth selected row with
   * respect to subspace spanned by first n-1 selected rows, (3) vector of
   * pivots
   *
   * \note The symmetrization condition is that if A(i,:), the ith row of A, is
   * selected as a pivot, then A(m-i-1,:) is also selected as a pivot. Here, m
   * is the row dimension of A, and A is zero-indexed. m must be even.
   */

  // Type T must be scalar-valued rank 2 array/array_view or matrix/matrix_view
  template <nda::MemoryArrayOfRank<2> T, nda::Scalar S = nda::get_value_t<T>>
  std::tuple<typename T::regular_type, nda::vector<double>, nda::vector<int>> pivrgs_sym(T const &a, double eps) {

    auto _ = nda::range::all;

    // Get matrix dimensions
    auto [m, n] = a.shape();
    int maxrnk  = std::min(m, n);

    if (m % 2 != 0) { throw std::runtime_error("Input matrix must have even number of rows."); }

    // Copy input data, re-ordering rows to make symmetric rows adjacent.
    auto aa                    = typename T::regular_type(m, n);
    aa(nda::range(0, m, 2), _) = a(nda::range(0, m / 2), _);
    aa(nda::range(1, m, 2), _) = a(nda::range(m - 1, m / 2 - 1, -1), _);

    // Compute norms of rows of input matrix, and rescale eps tolerance
    auto norms   = nda::vector<double>(m);
    double epssq = eps * eps;
    for (int j = 0; j < m; ++j) { norms(j) = normsq(aa(j, _)); }

    // Begin pivoted double Gram-Schmidt procedure
    int jpiv                 = 0;
    double nrm               = 0;
    auto piv                 = nda::arange(0, m);
    piv(nda::range(0, m, 2)) = nda::arange(0, m / 2); // Re-order pivots to match re-ordered input matrix
    piv(nda::range(1, m, 2)) = nda::arange(m - 1, m / 2 - 1, -1);

    if (maxrnk % 2 != 0) { // If n < m and n is odd, decrease maxrnk to maintain symmetry
      maxrnk -= 1;
    }

    for (int j = 0; j < maxrnk; j += 2) {

      // Find next pair of pivots
      jpiv = j;
      nrm  = norms(j) + norms(j + 1);
      for (int k = j + 2; k < m; k += 2) {
        if (norms(k) + norms(k + 1) > nrm) {
          jpiv = k;
          nrm  = norms(k) + norms(k + 1);
        }
      }

      // Swap current row pair with chosen pivot row pair
      deep_swap(aa(j, _), aa(jpiv, _));
      deep_swap(aa(j + 1, _), aa(jpiv + 1, _));
      std::swap(norms(j), norms(jpiv));
      std::swap(norms(j + 1), norms(jpiv + 1));
      std::swap(piv(j), piv(jpiv));
      std::swap(piv(j + 1), piv(jpiv + 1));

      // Orthogonalize current row (now the first chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j; ++k) { aa(j, _) = aa(j, _) - aa(k, _) * nda::blas::dotc(aa(k, _), aa(j, _)); }

      // Get norm of current row
      nrm = normsq(aa(j, _));

      // Terminate if sufficiently small, and return previously selected rows
      // (not including current row)
      if (nrm <= epssq) { return {aa(nda::range(0, j), _), norms(nda::range(0, j)), piv(nda::range(0, j))}; };

      // Normalize current row
      aa(j, _) /= sqrt(nrm);

      // Orthogonalize remaining rows against current row
      for (int k = j + 1; k < m; ++k) {
        if (norms(k) <= epssq) { continue; } // Can skip rows with norm less than tolerance
        aa(k, _) = aa(k, _) - aa(j, _) * nda::blas::dotc(aa(j, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }

      // Orthogonalize current row (now the second chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j + 1; ++k) { aa(j + 1, _) = aa(j + 1, _) - aa(k, _) * nda::blas::dotc(aa(k, _), aa(j + 1, _)); }

      // Normalize current row
      aa(j + 1, _) /= sqrt(normsq(aa(j + 1, _)));

      // Orthogonalize remaining rows against current row
      for (int k = j + 2; k < m; ++k) {
        if (norms(k) <= epssq) { continue; } // Can skip rows with norm less than tolerance
        aa(k, _) = aa(k, _) - aa(j + 1, _) * nda::blas::dotc(aa(j + 1, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }
    }

    return {aa(nda::range(maxrnk), _), norms(nda::range(maxrnk)), piv(nda::range(maxrnk))};
  }

  /** 
   * @brief Symmetrized pivoted reorthogonalized Gram-Schmidt with specified rank
   *
   * Return the leading r vectors of an orthogonal basis of the row space,
   * enforcing symmetrization of pivots.
   *
   * This is a translation of the Fortran subroutine "qrdgrm" by V.  Rokhlin,
   * with symmetrization added.
   *
   * @param a   Matrix to be orthogonalized
   * @param r   Rank cutoff
   *
   * @return Tuple of (1) matrix whose rows are leading @p r vectors in
   * orthogonal basis of row space of @p a, (2) vector with entry n given by the
   * squared l2 norm of the orthogonal complement of nth selected row with
   * respect to subspace spanned by first n-1 selected rows, (3) vector of
   * pivots
   *
   * \note The symmetrization condition is that if A(i,:), the ith row of A, is
   * selected as a pivot, then A(m-i-1,:) is also selected as a pivot. Here, m
   * is the row dimension of A, and A is zero-indexed. A can have an odd number of
   * rows if and only if @p r is odd, and in this case the middle row (index
   * (m-1)/2) of A is automatically selected as a pivot.
   */

  // Type T must be scalar-valued rank 2 array/array_view or matrix/matrix_view
  template <nda::MemoryArrayOfRank<2> T, nda::Scalar S = nda::get_value_t<T>>
  std::tuple<typename T::regular_type, nda::vector<double>, nda::vector<int>> pivrgs_sym(T const &a, int r) {

    auto _ = nda::range::all;

    // Get matrix dimensions
    auto [m, n] = a.shape();

    if (m % 2 == 1 && r % 2 == 0) { throw std::runtime_error("If input matrix has odd number of rows, r must be odd."); }
    if (r % 2 == 1 && m % 2 == 0) { throw std::runtime_error("If r is odd, input matrix must have odd number of rows."); }
    if (r > m || r > n + 1) { throw std::runtime_error("r must be less than or equal to min(m,n+1)."); }
    if (r == n + 1 && (n % 2 == 1 || n > m)) { throw std::runtime_error("If r = n+1, n must be even and less than or equal to m."); }

    // Copy input data, re-ordering rows to make symmetric rows adjacent. If m
    // odd, put middle row at the end.
    auto aa = typename T::regular_type(m, n);
    if (m % 2 == 0) {
      aa(nda::range(0, m, 2), _) = a(nda::range(0, m / 2), _);
      aa(nda::range(1, m, 2), _) = a(nda::range(m - 1, m / 2 - 1, -1), _);
    } else {
      aa(0, _)                   = a((m - 1) / 2, _);
      aa(nda::range(1, m, 2), _) = a(nda::range(0, (m - 1) / 2), _);
      aa(nda::range(2, m, 2), _) = a(nda::range(m - 1, (m - 1) / 2, -1), _);
      //aa(m - 1, _)                   = a((m - 1) / 2, _);
    }

    // Compute norms of rows of input matrix
    auto norms = nda::vector<double>(m);
    for (int j = 0; j < m; ++j) { norms(j) = normsq(aa(j, _)); }

    // Begin pivoted double Gram-Schmidt procedure
    int jpiv   = 0;
    double nrm = 0;
    auto piv   = nda::arange(0, m);
    if (m % 2 == 0) {
      piv(nda::range(0, m, 2)) = nda::arange(0, m / 2); // Re-order pivots to match re-ordered input matrix
      piv(nda::range(1, m, 2)) = nda::arange(m - 1, m / 2 - 1, -1);
    } else {
      piv(0)                   = (m - 1) / 2;
      piv(nda::range(1, m, 2)) = nda::arange(0, (m - 1) / 2);
      piv(nda::range(2, m, 2)) = nda::arange(m - 1, (m - 1) / 2, -1);
      //piv(m - 1)                   = (m - 1) / 2;
    }

    // If m odd, first choose middle row (now last row) as first pivot

    if (m % 2 == 1) {
      // Normalize
      aa(0, _) /= sqrt(normsq(aa(0, _)));

      // Orthogonalize remaining rows against current row
      for (int k = 1; k < m; ++k) {
        aa(k, _) = aa(k, _) - aa(0, _) * nda::blas::dotc(aa(0, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }
    }

    int jstrt = (m % 2 == 1) ? 1 : 0;
    // Then proceed with pivoted GS algorithm as normal
    for (int j = jstrt; j < r; j += 2) {

      // Find next pair of pivots
      jpiv = j;
      nrm  = norms(j) + norms(j + 1);
      for (int k = j + 2; k < m; k += 2) {
        if (norms(k) + norms(k + 1) > nrm) {
          jpiv = k;
          nrm  = norms(k) + norms(k + 1);
        }
      }

      // Swap current row pair with chosen pivot row pair
      deep_swap(aa(j, _), aa(jpiv, _));
      deep_swap(aa(j + 1, _), aa(jpiv + 1, _));
      std::swap(norms(j), norms(jpiv));
      std::swap(norms(j + 1), norms(jpiv + 1));
      std::swap(piv(j), piv(jpiv));
      std::swap(piv(j + 1), piv(jpiv + 1));

      // Orthogonalize current row (now the first chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j; ++k) { aa(j, _) = aa(j, _) - aa(k, _) * nda::blas::dotc(aa(k, _), aa(j, _)); }

      // Normalize current row
      aa(j, _) /= sqrt(normsq(aa(j, _)));

      // Orthogonalize remaining rows against current row
      for (int k = j + 1; k < m; ++k) {
        aa(k, _) = aa(k, _) - aa(j, _) * nda::blas::dotc(aa(j, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }

      // Orthogonalize current row (now the second chosen pivot row) against all
      // previously chosen rows
      for (int k = 0; k < j + 1; ++k) { aa(j + 1, _) = aa(j + 1, _) - aa(k, _) * nda::blas::dotc(aa(k, _), aa(j + 1, _)); }

      // Normalize current row
      aa(j + 1, _) /= sqrt(normsq(aa(j + 1, _)));

      // Orthogonalize remaining rows against current row
      for (int k = j + 2; k < m; ++k) {
        aa(k, _) = aa(k, _) - aa(j + 1, _) * nda::blas::dotc(aa(j + 1, _), aa(k, _));
        norms(k) = normsq(aa(k, _));
      }
    }

    return {aa(nda::range(r), _), norms(nda::range(r)), piv(nda::range(r))};
  }

  /**
  * @brief Get grid of equispaced points on [0,1] in relative time format
  *
  * @param n  Number of points
  *
  * @return Vector of equispaced points on [0,1] in relative time format
  */
  nda::vector<double> eqptsrel(int n);

  /**
  * @brief Convert points on [0,1] from relative to absolute time format
  *
  * @param t  Vector of points on [0,1] in relative time format
  *
  * @return Vector of points on [0,1] in absolute time format
  *
  * @note Converting a point from relative to absolute time format will, in
  * general, result in a loss of relative accuracy in the location of the point
  * if the point is close to t = 1. For example, in three-digit arithmetic, the
  * point t = 0.999111 could be stored in relative format as t^* = -0.889e-3,
  * but only as t = 0.999 in absolute format. For more information on the
  * relative time format used in cppdlr, please see the Background section of
  * the cppdlr documentation.
  */
  nda::vector<double> rel2abs(nda::vector_const_view<double> t);

  /**
  * @copydoc rel2abs(nda::vector_const_view<double>)
  */
  double rel2abs(double t);

  /**
  * @brief Convert points on [0,1] from absolute to relative time format
  *
  * @param t_abs Vector of points on [0,1] in absolute time format
  *
  * @return Vector of points on [0,1] in relative time format
  *
  * @note cppdlr uses the relative time format to describe imaginary time
  * points. Therefore, if you wish to specify an imaginary time point in the
  * standard absolute time format, for example to specify a point at which to
  * evaluate a DLR expansion, you must first convert the point to the relative
  * time format using this function. However, in order to maintain full relative
  * precision in all calculations, you must conform to the cppdlr standard and
  * specify point in the relative time format from the beginning. On the other
  * hand, in most cases only a mild loss of accuracy will result from beginning
  * with the absolute format and then converting to the relative format. For
  * more information on the relative time format used in cppdlr, please see the
  * Background section of the cppdlr documentation. 
  */
  nda::vector<double> abs2rel(nda::vector_const_view<double> t_abs);

  /**
  * @copydoc rel2abs(nda::vector_const_view<double>)
  */
  double abs2rel(double t_abs);

  /**
  * @brief Get real-valued type corresponding to type of given nda MemoryArray
  */
  template <nda::MemoryArray T> using make_real_t = decltype(make_regular(nda::real(std::declval<T>())));

  /**
  * @brief Get complex-valued type corresponding to type of given nda MemoryArray
  */
  template <nda::MemoryArray T> using make_cplx_t = decltype(make_regular(std::declval<T>() * 1i));

  /**
  * @brief Get type of given nda MemoryArray with scalar value type replaced by
  * common type of two given types (real if both are real, complex otherwise)
  */
  template <nda::MemoryArray T, nda::Scalar S1, nda::Scalar S2> struct make_common_helper {
    using S              = std::common_type_t<S1, S2>;
    static constexpr S x = 0;
    using type           = decltype(make_regular(std::declval<T>() * x));
  };
  template <nda::MemoryArray T, nda::Scalar S1, nda::Scalar S2> using make_common_t = typename make_common_helper<T, S1, S2>::type;

  /**
  * @brief Contract the last dimension of an array a with the first dimension of
  * an array b
  *
  * @param a  An array/matrix/vector or array/matrix/vector view of rank at least 2
  * @param b  An array/matrix/vector or array/matrix/vector view of rank at least 2
  *
  * @return Contraction of the inner dimensions of \p a and \p b
  */
  template <nda::MemoryArray Ta, nda::MemoryArray Tb, nda::Scalar Sa = nda::get_value_t<Ta>, nda::Scalar Sb = nda::get_value_t<Tb>,
            nda::Scalar S = std::common_type_t<Sa, Sb>>
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
    auto a_reshaped = nda::reshape(a, m, p);
    auto b_reshaped = nda::reshape(b, p, n);

    // Get shape of output array
    auto c_shape = std::array<int, ra + rb - 2>();
    for (int i = 0; i < ra - 1; ++i) { c_shape[i] = a.shape(i); }
    for (int i = ra - 1; i < ra + rb - 2; ++i) { c_shape[i] = b.shape(i - ra + 2); }

    // Compute the contraction, reshape, and return
    return reshape(matmul(a_reshaped, b_reshaped), c_shape);
  }

  /**
  * @brief Quick and dirty adaptive Gauss quadrature
  *
  * This function implements adaptive Gauss-Legendre quadrature with local error
  * estimation only, using a stack.
  * 
  * @param[in] f  Function to be integrated
  * @param[in] a  Lower integration limit
  * @param[in] b  Upper integration limit
  * @param[in] tol  Absolute error tolerance
  * @param[in] xgl Gauss-Legendre nodes
  * @param[in] wgl Gauss-Legendre weights
  *
  * @return Integral of \p f from \p a to \p b
  *
  * \note This is a quick and dirty adaptive integration function, which tries
  * to achieve an error tolerance \p tol but doesn't guarantee it. A more robust
  * implementation would use global error estimation. Nevertheless, this works
  * quite well most of the time.
  */

  template <typename S>
  // TODO: require S is scalar
  S adapgl(std::function<nda::array<S, 1>(nda::array<double, 1>)> f, double a, double b, double tol, nda::vector<double> xgl,
           nda::vector<double> wgl) {

    int maxnde = 100000; // Max nodes in integration tree
    int maxstk = 1000;   // Max stack size

    S sum = 0; // Integral

    // Initialize stack: there is a stack of interval endpoints associated w/
    // each node in the integration tree, and a corresponding stack of values
    // associated w/ each such interval

    auto endpt = nda::array<double, 2>(maxstk, 2); // Stack of interval endpoints
    auto sums  = nda::array<S, 1>(maxstk);         // Stack of sums associated w/ each interval

    int istk = 0; // Stack pointer
    int nnde = 1; // Number of nodes in integration tree

    endpt(istk, 0) = a; // Root node interval is [a,b]
    endpt(istk, 1) = b;

    double h   = (b - a) / 2;                                        // Interval half-width
    sums(istk) = h * nda::blas::dotc(wgl, f(xgl * h + (a + b) / 2)); // Integral of f on root node interval

    double aa = 0, bb = 0, cc = 0; // Interval endpoints, midpoint
    S sum1 = 0, sum2 = 0;          // Integrals on subintervals

    // Main loop
    while (true) {

      aa = endpt(istk, 0); // Get interval endpoints
      bb = endpt(istk, 1);
      cc = (aa + bb) / 2; // Get midpoint

      // Get integral on each subinterval
      h    = (cc - aa) / 2;
      sum1 = h * nda::blas::dotc(wgl, f(xgl * h + (aa + cc) / 2));
      sum2 = h * nda::blas::dotc(wgl, f(xgl * h + (cc + bb) / 2));

      // If error on subintervals is less than tolerance, add to integral and
      // move to next node

      if (std::abs(sum1 + sum2 - sums(istk)) < tol) {
        sum += sum1 + sum2;           // Add to subintervals to integral
        istk--;                       // Remove current interval from stack
        if (istk < 0) { return sum; } // If stack is empty, return integral
      } else {                        // Integral on current interval not converged

        // Remove current interval from stack, and add subintervals
        endpt(istk, 1)     = cc;
        endpt(istk + 1, 0) = cc;
        endpt(istk + 1, 1) = bb;
        sums(istk)         = sum1;
        sums(istk + 1)     = sum2;

        istk++;
        nnde += 2;

        if (istk > maxstk) { throw std::runtime_error("integration stack too large"); }
        if (nnde > maxnde) { throw std::runtime_error("integration tree too large"); }
      }
    }
  }

  /**
  * @brief Gauss-Legendre nodes and weights
  *
  * Uses Newton iteration to obtain the Gauss-Legendre nodes and weights
  *
  * @param[in] n  Number of nodes
  *
  * @return Tuple of nodes and weights
  */
  std::tuple<nda::vector<double>, nda::vector<double>> gaussquad(int n);

  /**
  * @brief Evaluate Legendre polynomial of degree n and its derivative
  *
  * Uses Legendre three-term recurrence
  *
  * @param[in] n  Degree of polynomial Pn(x)
  * @param[in] x  Point at which to evaluate polynomial
  *
  * @return Tuple of polynomial value Pn(x) and derivative Pn'(x)
  */
  std::tuple<double, double> leg_eval(int n, double x);

} // namespace cppdlr
