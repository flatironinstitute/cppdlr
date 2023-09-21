// Copyright (c) 2022-2023 Simons Foundation
// Copyright (c) 2023 Hugo Strand
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
// Authors: Hugo U. R. Strand, Nils Wentzell, Jason Kaye

#pragma once
#include <nda/nda.hpp>
#include "cppdlr/dlr_kernels.hpp"
#include "cppdlr/dlr_build.hpp"
#include "cppdlr/utils.hpp"
#include "nda/clef/clef.hpp"

#include <h5/h5.hpp>
#include <nda/h5.hpp>

namespace cppdlr {

  static constexpr auto _ = nda::range::all;

  /**
  * Option for ordinary or time-ordered convolution 
  */
   static constexpr bool ORDINARY = false, TIME_ORDERED = true;

  /**
  * @class imtime_ops
  * @brief Class responsible for all DLR imaginary time operations, including
  * building imaginary time grid and transformations.
  *
  * \note First dimension of all Green's function and coefficient arrays must be
  * DLR rank r.
  */

  class imtime_ops {

    public:
    /** 
    * @brief Constructor for imtime_ops
    * 
    * @param[in] lambda DLR cutoff parameter
    * @param[in] dlr_rf DLR frequencies
    * @param[in] statistic Particle statistics: Boson or Fermion
    * @param[in] symmetrize Non-symmetrized (false, default) or symmetrized DLR (true)
    *
    * @note In case symmetrize and Boson options are selected, we enforce that
    * tau = 0 is chosen as a DLR frequency. The reasoning is as follows: to
    * maintain symmetry of the bosonic DLR imaginary frequency nodes, we must
    * either forbid the zero imaginary frequency from being chosen, or always
    * choose it. We decide to always choose it, which means the number of
    * symmetrized imaginary frequency nodes will be odd. To maintain the same
    * number r of DLR imaginary time and imaginary frequency nodes, we must
    * therefore have an odd number of DLR imaginary time nodes. To obtain a
    * symmetric DLR imaginary time grid with an odd number of points, we must
    * choose tau = 1/2 as a DLR frequency.
    */
    imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf, statistic_t statistic, bool symmetrize);

    /** 
    * @brief Constructor for imtime_ops
    * 
    * @param[in] lambda DLR cutoff parameter
    * @param[in] dlr_rf DLR frequencies
    */
    imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf);

    imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf, nda::vector_const_view<double> dlr_it, nda::matrix_const_view<double> cf2it,
               nda::matrix_const_view<double> it2cf_lu, nda::vector_const_view<int> it2cf_piv)
       : lambda_(lambda), r(dlr_rf.size()), dlr_rf(dlr_rf), dlr_it(dlr_it), cf2it(cf2it), it2cf{it2cf_lu, it2cf_piv} {};

    imtime_ops() = default;

    /** 
    * @brief Transform values of Green's function G on DLR imaginary time grid to
    * DLR coefficients 
    *
    * @param[in] g Values of G on DLR imaginary time grid
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> typename T::regular_type vals2coefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape g to matrix w/ first dimension r
      auto g_rs = nda::reshape(g, r, g.size() / r);
      auto gct  = nda::matrix<S>(transpose(g_rs));

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs (we transpose because LAPACK requires index into RHS # to be
      // slowest)

      if constexpr (nda::have_same_value_type_v<T, decltype(it2cf.lu)>) {
        nda::lapack::getrs(it2cf.lu, gct, it2cf.piv);
      } else {
        // getrs requires matrix and rhs to have same value type
        nda::lapack::getrs(nda::matrix<S>(it2cf.lu), gct, it2cf.piv);
      }

      // Reshape to original dimensions and return
      auto gc = nda::matrix<S>(transpose(gct));
      return nda::reshape(gc, g.shape());
    }

    /** 
    * @brief Transform DLR coefficients of Green's function G to values on DLR
    * imaginary time grid
    *
    * @param[in] gc DLR coefficients of G
    *
    * @return Values of G on DLR imaginary time grid
    * */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> typename T::regular_type coefs2vals(T const &gc) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape gc to a matrix w/ first dimension r
      auto gc_rs = nda::reshape(gc, r, gc.size() / r);

      // Apply coeffs -> vals matrix
      auto g = cf2it * nda::matrix_const_view<S>(gc_rs);

      // Reshape to original dimensions and return
      return nda::reshape(g, gc.shape());
    }

    /** 
    * @brief Evaluate DLR expansion of G, given by its DLR coefficients, at imaginary
    * time point 
    *
    * @param[in] gc DLR coefficients of G
    * @param[in] t  Evaluation point, in relative format
    *
    * @return Value of G at @p t
    *
    * @note The given evaluation point must be scaled to the interval [0, 1]
    * (rather than [0, beta]) and then given in the relative time format. Please
    * see the "Imaginary time point format" section in the Background page of
    * the documentation for more information.
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> auto coefs2eval(T const &gc, double t) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions

      if constexpr (T::rank == 1) {

        // TODO: can be further optimized to reduce # exponential evals.
        // Evaluate DLR expansion
        S g = 0;
        if (t >= 0) {
          for (int l = 0; l < r; ++l) { g += k_it_abs(t, dlr_rf(l)) * gc(l); }
        } else {
          for (int l = 0; l < r; ++l) { g += k_it_abs(-t, -dlr_rf(l)) * gc(l); }
        }

        return g;
      } else {

        // Reshape gc to matrix w/ first dimension r
        auto gc_rs = nda::reshape(gc, r, gc.size() / r);

        // Get output shape
        std::array<long, T::rank - 1> shape_out;
        for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i + 1); }

        // Get vector of evaluation of DLR expansion at a point
        auto kvec = build_evalvec(t);

        // Evaluate DLR expansion, reshape to original dimensions (with first
        // dimension summed out), and return
        return nda::reshape(transpose(nda::matrix_const_view<S>(gc_rs)) * kvec, shape_out);
      }
    }

    /** 
    * @brief Build vector of evaluation of DLR expansion at an imaginary time point 
    *
    * @param[in] t  Evaluation point
    *
    * @return Vector of evaluation at @p t
    *
    * @note The given evaluation point must be scaled to the interval [0, 1]
    * (rather than [0, beta]) and then given in the relative time format. Please
    * see the "Imaginary time point format" section in the Background page of
    * the documentation for more information.
    **/
    nda::vector<double> build_evalvec(double t) const;

    /** 
    * @brief Obtain DLR coefficients of a Green's function G from scattered
    * imaginary time grid by least squares fitting
    *
    * @param[in] t Imaginary time grid points at which G is sampled, in relative format
    * @param[in] g Values of G on grid t
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>>
    typename T::regular_type fitvals2coefs(nda::vector_const_view<double> t, T const &g) const {

      // Check first dimension of g equal to length of t
      int n = t.size();
      if (n != g.shape(0)) throw std::runtime_error("First dim of g must be equal to length of t.");

      // Get matrix for least squares fitting: columns are DLR basis functions
      // evaluating at data points t. Must built in Fortran layout for
      // compatibility with LAPACK.
      auto kmat = nda::matrix<S, F_layout>(n, r); // Make sure matrix has same scalar type as g
      for (int j = 0; j < r; ++j) {
        for (int i = 0; i < n; ++i) { kmat(i, j) = k_it(t(i), dlr_rf(j)); }
      }

      // Reshape g to matrix w/ first dimension n, and put in Fortran layout for
      // compatibility w/ LAPACK
      auto g_rs = nda::matrix<S, F_layout>(nda::reshape(g, n, g.size() / n));

      // Solve least squares problem to obtain DLR coefficients
      auto s   = nda::vector<double>(r); // Singular values (not needed)
      int rank = 0;                      // Rank of system matrix (not needed)
      nda::lapack::gelss(kmat, g_rs, s, 0.0, rank);

      auto gc_shape                = g.shape();                          // Output shape is same as g...
      gc_shape[0]                  = r;                                  // ...with first dimension n replaced by r
      auto gc                      = typename T::regular_type(gc_shape); // Output array: output type is same as g
      reshape(gc, r, g.size() / n) = g_rs(nda::range(r), _);             // Place result in output array

      return gc;
    }

    /**
    * @brief Compute reflection of imaginary time Green's function
    *
    * The reflection of g(t) is g(beta - t). This method takes the values of a
    * Green's function g at the DLR imaginary time nodes to the values of its
    * reflection at the DLR imaginary time nodes.
    *
    * @param[in] g Values of g at DLR imaginary time nodes
    *
    * @return Values of g(beta - t) at DLR imaginary time nodes
    */
    template <nda::MemoryArray T> typename T::regular_type reflect(T const &g) const {

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Initialize reflection matrix, if it hasn't been done already
      if (refl.empty()) { reflect_init(); }

      if constexpr (T::rank == 1) { // Scalar-valued Green's function
        return matmul(refl, g);
      } else {
        auto gr                      = typename T::regular_type(g.shape());       // Output has same type/shape as g
        reshape(gr, r, g.size() / r) = matmul(refl, reshape(g, r, g.size() / r)); // Reshape, matrix multiply, reshape back
        return gr;
      }
    }


    /** 
    * @brief Compute convolution of two imaginary time Green's functions
    *
    * The convolution of f and g is defined as h(t) = (f * g)(t) = int_0^beta
    * f(t-t') g(t') dt', where fermionic/bosonic antiperiodicity/periodicity
    * are used to define the Green's functions on (-beta, 0).  This method takes
    * the DLR coefficients of f and g as input and returns the values of h on
    * the DLR imaginary time grid.
    *
    * By specifying the @p time_order flag, this method can be used to compute
    * the time-ordered convolution of f and g, defined as h(t) = (f * g)(t) =
    * int_0^tau f(t-t') g(t') dt'.
    *
    * The convolution is computed on-the-fly in O(r^2) operations using the
    * method described in Appendix A of 
    *
    * J. Kaye, H. U. R. Strand, D. Golez, "Decomposing imaginary time Feynman
    * diagrams using separable basis functions: Anderson impurity model strong
    * coupling expansion," arXiv:2307.08566 (2023).
    *
    * @param[in] beta Inverse temperature
    * @param[in] statistic Fermionic ("Fermion" or 0) or bosonic ("Boson" or 1)
    * @param[in] fc DLR coefficients of f
    * @param[in] gc DLR coefficients of g
    * @param[in] time_order Flag for ordinary (false or ORDINARY, default) or
    * time-ordered (true or TIME_ORDERED) convolution
    *
    * @return Values of h = f * g on DLR imaginary time grid
    * */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>>
    typename T::regular_type convolve(double beta, statistic_t statistic, T const &fc, T const &gc, bool time_order = false) const {

      if (r != fc.shape(0) || r != gc.shape(0)) throw std::runtime_error("First dim of input arrays must be equal to DLR rank r.");
      if (fc.shape() != gc.shape()) throw std::runtime_error("Input arrays must have the same shape.");

      // TODO: implement bosonic case and remove
      if (statistic == 0) throw std::runtime_error("imtime_ops::convolve not yet implemented for bosonic Green's functions.");

      // Initialize convolution, if it hasn't been done already
      if (!time_order & hilb.empty()) { convolve_init(); }
      if (time_order & thilb.empty()) { tconvolve_init(); }

      // Get view of helper matrices based on time_order flag
      auto hilb_v   = (time_order ? nda::matrix_const_view<double>(thilb) : nda::matrix_const_view<double>(hilb));
      auto tcf2it_v = (time_order ? nda::matrix_const_view<double>(ttcf2it) : nda::matrix_const_view<double>(tcf2it));

      if constexpr (T::rank == 1) { // Scalar-valued Green's function

        // Take array view of fc and gc
        auto fca = nda::array_const_view<S, 1>(fc);
        auto gca = nda::array_const_view<S, 1>(gc);

        // Diagonal contribution
        auto h = arraymult(tcf2it_v, make_regular(fca * gca));

        // Off-diagonal contribution
        auto tmp = fca * arraymult(hilb_v, gca) + gca * arraymult(hilb_v, fca);
        return beta * (h + arraymult(cf2it, make_regular(tmp)));

      } else if (T::rank == 3) { // Matrix-valued Green's function

        // Diagonal contribution
        auto fcgc = nda::array<S, 3>(fc.shape()); // Product of coefficients of f and g
        for (int i = 0; i < r; ++i) { fcgc(i, _, _) = arraymult(fc(i, _, _), gc(i, _, _)); }
        auto h = arraymult(tcf2it_v, fcgc);

        // Off-diagonal contribution
        auto tmp1 = arraymult(hilb_v, fc);
        auto tmp2 = arraymult(hilb_v, gc);
        for (int i = 0; i < r; ++i) { tmp1(i, _, _) = arraymult(tmp1(i, _, _), gc(i, _, _)) + arraymult(fc(i, _, _), tmp2(i, _, _)); }

        return beta * (h + arraymult(cf2it, tmp1));

      } else {
        throw std::runtime_error("Input arrays must be rank 1 (scalar-valued Green's function) or 3 (matrix-valued Green's function).");
      }
    }

    /** 
    * @brief Compute convolution of two imaginary time Green's functions,
    * given matrix of convolution by one of them
    *
    * The convolution of f and g is defined as h(t) = (f * g)(t) = int_0^beta
    * f(t-t') g(t') dt', where fermionic/bosonic antiperiodicity/periodicity are
    * used to define the Green's functions on (-beta, 0).  This method takes the
    * matrix of convolution by f (computed by the convmat method) and the values
    * of g on the DLR imaginary time grid, and returns the values of h on the
    * DLR imaginary time grid
    *
    * By passing in the matrix of time-ordered convolution by f (computed by the
    * convmat method with the @p time_order flag set to true, or TIME_ORDERED),
    * this method can be used to compute the time-ordered convolution of f and
    * g, defined as h(t) = (f * g)(t) = int_0^tau f(t-t') g(t') dt'.
    *
    * @param[in] fconv Matrix of convolution by f
    * @param[in] g Values of g on the DLR imaginary time grid
    *
    * @return Values of h = f * g on DLR imaginary time grid
    * */
    template <nda::MemoryMatrix Tf, nda::MemoryArray Tg, nda::Scalar Sf = nda::get_value_t<Tf>, nda::Scalar Sg = nda::get_value_t<Tg>,
              nda::Scalar S = typename std::common_type<Sf, Sg>::type>
    typename Tg::regular_type convolve(Tf const &fconv, Tg const &g) const {

      if (r != g.shape(0)) throw std::runtime_error("First dim of input g must be equal to DLR rank r.");

      if constexpr (Tg::rank == 1) { // Scalar-valued Green's function

        if (fconv.shape(1) != g.shape(0)) throw std::runtime_error("Input array dimensions incompatible.");

        return arraymult(fconv, g);

      } else if (Tg::rank == 3) { // Matrix-valued Green's function

        if (fconv.shape(1) != g.shape(0) * g.shape(1)) throw std::runtime_error("Input array dimensions incompatible.");

        int norb1 = g.shape(1);
        int norb2 = g.shape(2);

        auto h                       = nda::array<S, 3>(r, norb1, norb2);
        reshape(h, r * norb1, norb2) = arraymult(fconv, reshape(g, r * norb1, norb2));

        return h;

      } else {
        throw std::runtime_error("Input arrays must be rank 1 (scalar-valued Green's function) or 3 (matrix-valued Green's function).");
      }
    }

    /** 
    * @brief Compute matrix of convolution by an imaginary time Green's function
    *
    * The convolution of f and g is defined as h(t) = (f * g)(t) = int_0^beta
    * f(t-t') g(t') dt', where fermionic/bosonic antiperiodicity/periodicity are
    * used to define the Green's functions on (-beta, 0).  This method takes the
    * DLR coefficients of f as input and returns the matrix of convolution by f.
    * This matrix can be applied to the values of g on the DLR imaginary time
    * grid, to produce the values of h on the DLR imaginary time grid.
    *
    * By specifying the @p time_order flag, this method can be used to compute
    * the time-ordered convolution of f and g, defined as h(t) = (f * g)(t) =
    * int_0^tau f(t-t') g(t') dt'.
    *
    * The convolution matrix is constructed using the method described in
    * Appendix A of 
    *
    * J. Kaye, H. U. R. Strand, D. Golez, "Decomposing imaginary time Feynman
    * diagrams using separable basis functions: Anderson impurity model strong
    * coupling expansion," arXiv:2307.08566 (2023).
    *
    * @param[in] beta Inverse temperature
    * @param[in] statistic Fermionic ("Fermion" or 0) or bosonic ("Boson" or 1)
    * @param[in] fc DLR coefficients of f
    * @param[in] time_order Flag for ordinary (false or ORDINARY, default) or
    * time-ordered (true or TIME_ORDERED) convolution
    *
    * @return Matrix of convolution by f
    *
    * \note Whereas the method imtime_ops::convolve takes the DLR coefficients
    * of f and g as input and computes their convolution h directly, this method
    * returns a matrix which should be applied to the DLR imaginary time grid values
    * of g, rather than its DLR coefficients, in to order to obtain the
    * convolution h. The purpose of this is to make the input and output
    * representations of the convolution matrix equal, which is often convenient
    * in practice.
    *
    * \note In the case of matrix-valued Green's functions, we think of the
    * matrix of convolution by f as an r*norb x r*norb matrix, or a block r x r
    * matrix of norb x norb blocks. Here r is the DLR rank and norb is the
    * number of orbital indices. This matrix would then be applied to a Green's
    * function g, represented as an r*norb x norb matrix, or a block r x 1
    * matrix of norb x norb blocks.
    * */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>>
    nda::matrix<S> convmat(double beta, statistic_t statistic, T const &fc, bool time_order = false) const {

      if (r != fc.shape(0)) throw std::runtime_error("First dim of input array must be equal to DLR rank r.");

      // TODO: implement bosonic case and remove
      if (statistic == 0) throw std::runtime_error("imtime_ops::convmat not yet implemented for bosonic Green's functions.");

      // Initialize convolution, if it hasn't been done already
      if (!time_order & hilb.empty()) { convolve_init(); }
      if (time_order & thilb.empty()) { tconvolve_init(); }

      // Get view of helper matrices based on time_order flag
      auto hilb_v   = (time_order ? nda::matrix_const_view<double>(thilb) : nda::matrix_const_view<double>(hilb));
      auto tcf2it_v = (time_order ? nda::matrix_const_view<double>(ttcf2it) : nda::matrix_const_view<double>(tcf2it));

      if constexpr (T::rank == 1) { // Scalar-valued Green's function

        // First construct convolution matrix from DLR coefficients to DLR grid
        // values
        auto fconv = nda::matrix<S>(r, r); // Matrix of convolution by f

        // Diagonal contribution (given by diag(tau_k) * K(tau_k, om_l) * diag(fc_l))
        for (int k = 0; k < r; ++k) {
          for (int l = 0; l < r; ++l) { fconv(k, l) = tcf2it_v(k, l) * fc(l); }
        }

        // Off-diagonal contribution (given by K * (diag(hilb*fc) +
        // (diag(fc)*hilb)), where K is the matrix K(dlr_it(k), dlr_rf(l))
        auto tmp1 = arraymult(hilb_v, fc); // hilb * fc
        auto tmp2 = nda::matrix<S>(r, r);
        for (int k = 0; k < r; ++k) {
          for (int l = 0; l < r; ++l) {
            tmp2(k, l) = fc(k) * hilb_v(k, l); // diag(fc) * hilb
          }
          tmp2(k, k) += tmp1(k); // diag(fc)*hilb + diag(hilb*fc)
        }
        fconv += arraymult(cf2it, tmp2);

        // Then precompose with DLR grid values to DLR coefficients matrix
        if constexpr (nda::have_same_value_type_v<T, decltype(it2cf.lu)>) {
          nda::lapack::getrs(transpose(it2cf.lu), fconv, it2cf.piv); // Note: lapack effectively tranposes fconv by fortran reordering here
        } else {
          // getrs requires matrix and rhs to have same value type
          nda::lapack::getrs(transpose(nda::matrix<S>(it2cf.lu)), fconv, it2cf.piv);
        }

        return beta * fconv;

      } else if (T::rank == 3) { // Matrix-valued Green's function

        int norb1 = fc.shape(1);
        int norb2 = fc.shape(2);

        // First construct convolution matrix from DLR coefficients to DLR grid
        // values
        auto fconv    = nda::matrix<S>(r * norb1, r * norb2);    // Matrix of convolution by f
        auto fconv_rs = nda::reshape(fconv, r, norb1, r, norb2); // Array view to index into fconv for conevenience

        // Diagonal contribution (given by diag(tau_k) * K(tau_k, om_l) * diag(fc_l))
        for (int k = 0; k < r; ++k) {
          for (int l = 0; l < r; ++l) { fconv_rs(k, _, l, _) = tcf2it_v(k, l) * fc(l, _, _); }
        }

        // Off-diagonal contribution (given by K * (diag(hilb*fc) +
        // (diag(fc)*hilb)), where K is the matrix K(dlr_it(k), dlr_rf(l))
        auto tmp1 = arraymult(hilb_v, fc); // hilb * fc
        auto tmp2 = nda::array<S, 4>(r, norb1, r, norb2);
        for (int k = 0; k < r; ++k) {
          for (int l = 0; l < r; ++l) {
            tmp2(k, _, l, _) = fc(k, _, _) * hilb_v(k, l); // diag(fc) * hilb
          }
          tmp2(k, _, k, _) += tmp1(k, _, _); // diag(fc)*hilb + diag(hilb*fc)
        }
        fconv_rs += arraymult(cf2it, tmp2);

        // Then precompose with DLR grid values to DLR coefficients matrix

        // Transpose last two indices to put make column DLR index the last
        // index
        auto fconvtmp    = nda::matrix<S>(r * norb1 * norb2, r);
        auto fconvtmp_rs = nda::reshape(fconvtmp, r, norb1, norb2, r);
        for (int i = 0; i < norb2; ++i) {
          for (int k = 0; k < r; ++k) { fconvtmp_rs(_, _, i, k) = fconv_rs(_, _, k, i); }
        }

        // Do the solve
        if constexpr (nda::have_same_value_type_v<T, decltype(it2cf.lu)>) {
          nda::lapack::getrs(transpose(it2cf.lu), fconvtmp, it2cf.piv); // Note: lapack effectively tranposes fconv by fortran reordering here
        } else {
          // getrs requires matrix and rhs to have same value type
          nda::lapack::getrs(transpose(nda::matrix<S>(it2cf.lu)), fconvtmp, it2cf.piv);
        }

        // Transpose back
        for (int i = 0; i < norb2; ++i) {
          for (int k = 0; k < r; ++k) { fconv_rs(_, _, k, i) = fconvtmp_rs(_, _, i, k); }
        }

        return beta * fconv;

      } else {
        throw std::runtime_error("Input arrays must be rank 1 (scalar-valued Green's function) or 3 (matrix-valued Green's function).");
      }
    }

    /** 
    * @brief Get DLR imaginary time nodes
    *
    * @return DLR imaginary time nodes
    */
    nda::vector_const_view<double> get_itnodes() const { return dlr_it; };
    double get_itnodes(int i) const { return dlr_it(i); };

    /** Access DLR imaginary real frequency nodes*/
    /**
    * @brief Get DLR real frequency nodes
    *
    * @return DLR real frequency nodes
    */
    nda::vector_const_view<double> get_rfnodes() const { return dlr_rf; };
    double get_rfnodes(int i) const { return dlr_rf(i); };

    /**
    * @brief Get transformation matrix from DLR coefficients to values at DLR imaginary time nodes
    *
    * @return Transformation matrix
    */
    nda::matrix_const_view<double> get_cf2it() const { return cf2it; };

    /**
    * @brief Get LU factors of transformation matrix from DLR imaginary time values to coefficients
    *
    * @return LU factors
    */
    nda::matrix_const_view<double> get_it2cf_lu() const { return it2cf.lu; };

    /**
    * @brief Get LU pivots of transformation matrix from DLR imaginary time values to coefficients
    *
    * @return LU pivots
    */
    nda::vector_const_view<int> get_it2cf_piv() const { return it2cf.piv; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }
    double lambda() const { return lambda_; }

    /**
    * @brief Initialization for convolution methods
    *
    * Initialize matrices required for the convolution methods. This method is
    * called automatically the first time one of the relevant convolution
    * methods is called, but it may also be called manually to avoid the
    * additional overhead in the first convolution call.
    */
    void convolve_init() const {

      hilb   = nda::matrix<double>(r, r);
      tcf2it = nda::matrix<double>(r, r);

      // "Discrete Hilbert transform" matrix -(1-delta_jk)/(dlr_rf(j) -
      // dlr_rf(k)), scaled by beta
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (j == k) {
            hilb(j, k) = 0;
          } else {
            hilb(j, k) = 1.0 / (dlr_rf(j) - dlr_rf(k));
          }
        }
      }

      // Matrix which applies DLR coefficients to imaginary time grid values
      // transformation matrix, and then multiplies the result by tau, the
      // imaginary time variable
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (dlr_it(j) > 0) {
            tcf2it(j, k) = -(dlr_it(j) + k_it(1.0, dlr_rf(k))) * cf2it(j, k);
          } else {
            tcf2it(j, k) = -(dlr_it(j) - k_it(0.0, dlr_rf(k))) * cf2it(j, k);
          }
        }
      }
    }

    /**
    * @brief Initialization for time-ordered convolution methods
    *
    * Initialize matrices required for the time-ordered convolution methods.
    * This method is called automatically the first time one of the relevant
    * time-ordered convolution methods is called, but it may also be called
    * manually to avoid the additional overhead in the first time-ordered
    * convolution call. 
    */
    void tconvolve_init() const {

      thilb   = nda::matrix<double>(r, r);
      ttcf2it = nda::matrix<double>(r, r);

      // "Discrete Hilbert transform" matrix
      // -(1-delta_jk)*K(0,dlr_rf(k))/(dlr_rf(j) - dlr_rf(k)) for time-ordered
      // convolution, scaled by beta
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (j == k) {
            thilb(j, k) = 0;
          } else {
            thilb(j, k) = k_it(0.0, dlr_rf(k)) / (dlr_rf(k) - dlr_rf(j));
          }
        }
      }

      // Matrix which applies K(0,dlr_rf(j)) multiplication, then DLR
      // coefficients to imaginary time grid values transformation matrix, and
      // then multiplies the result by tau, the imaginary time variable
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (dlr_it(j) > 0) {
            ttcf2it(j, k) = dlr_it(j) * cf2it(j, k) * k_it(0.0, dlr_rf(k));
          } else {
            ttcf2it(j, k) = (1 + dlr_it(j)) * cf2it(j, k) * k_it(0.0, dlr_rf(k));
          }
        }
      }
    }

    /**
    * @brief Initialization for reflection method
    *
    * Initialize matrix of the reflection G(tau) -> G(beta - tau) acting on the
    * values of the Green's function at the DLR imaginary time nodes.  This
    * matrix is required for the reflect method. It is called automatically the
    * first time the reflect method is called, but it may also be called
    * manually to avoid the additional overhead of the first call to reflect.
    */
    void reflect_init() const {

      refl = nda::matrix<double>(r, r);

      // Matrix of reflection acting on DLR coefficients and returning values at
      // DLR nodes
      for (int i = 0; i < r; ++i) {
        for (int j = 0; j < r; ++j) { refl(i, j) = k_it(-dlr_it(i), dlr_rf(j)); }
      }

      // Precompose with DLR values to coefficients matrix
      nda::lapack::getrs(transpose(it2cf.lu), refl, it2cf.piv); // Lapack effectively transposes refl by fortran reordering here
    }

    private:
    double lambda_;
    int r;                      ///< DLR rank
    nda::vector<double> dlr_rf; ///< DLR frequencies
    nda::vector<double> dlr_it; ///< DLR imaginary time nodes
    nda::matrix<double> cf2it;  ///< Transformation matrix from DLR coefficients to values at DLR imaginary time nodes

    /**
    * @brief Struct for transformation from DLR imaginary time values to coefficients
    */
    struct {
      nda::matrix<double> lu; ///< LU factors (LAPACK format) of imaginary time vals -> coefs matrix
      nda::vector<int> piv;   ///< LU pivots (LAPACK format) of imaginary time vals -> coefs matrix
    } it2cf;

    // Arrays used for dlr_imtime::convolve
    mutable nda::matrix<double> hilb;   ///< "Discrete Hilbert transform" matrix
    mutable nda::matrix<double> tcf2it; ///< A matrix required for convolution

    // Arrays used for dlr_imtime::tconvolve
    mutable nda::matrix<double> thilb;   ///< "Discrete Hilbert transform" matrix, modified for time-ordered convolution
    mutable nda::matrix<double> ttcf2it; ///< A matrix required for time-ordered convolution

    // Arrays used for dlr_imtime::reflect
    mutable nda::matrix<double> refl; ///< Matrix of reflection

    // -------------------- hdf5 -------------------

    public:
    static std::string hdf5_format() { return "cppdlr::imtime_ops"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, imtime_ops const &m) {

      h5::group gr = fg.create_group(subgroup_name);
      write_hdf5_format(gr, m);

      h5::write(gr, "lambda", m.lambda());
      h5::write(gr, "rf", m.get_rfnodes());
      h5::write(gr, "it", m.get_itnodes());
      h5::write(gr, "cf2it", m.get_cf2it());
      h5::write(gr, "it2cf_lu", m.get_it2cf_lu());
      h5::write(gr, "it2cf_piv", m.get_it2cf_piv());
    }

    friend void h5_read(h5::group fg, std::string const &subgroup_name, imtime_ops &m) {

      h5::group gr = fg.open_group(subgroup_name);
      assert_hdf5_format(gr, m);

      auto lambda    = h5::read<double>(gr, "lambda");
      auto rf        = h5::read<nda::vector<double>>(gr, "rf");
      auto it        = h5::read<nda::vector<double>>(gr, "it");
      auto cf2it_     = h5::read<nda::matrix<double>>(gr, "cf2it");
      auto it2cf_lu  = h5::read<nda::matrix<double>>(gr, "it2cf_lu");
      auto it2cf_piv = h5::read<nda::vector<int>>(gr, "it2cf_piv");

      m = imtime_ops(lambda, rf, it, cf2it_, it2cf_lu, it2cf_piv);
    }
  };

} // namespace cppdlr
