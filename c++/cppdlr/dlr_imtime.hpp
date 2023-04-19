#pragma once
#include <nda/nda.hpp>
#include "cppdlr/dlr_kernels.hpp"
#include "cppdlr/dlr_build.hpp"
#include "cppdlr/utils.hpp"

#include <h5/h5.hpp>
#include <nda/h5.hpp>

namespace cppdlr {

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
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
      requires(nda::is_scalar_v<S>)
    typename T::regular_type vals2coefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape g to matrix w/ first dimension r
      auto g_rs = nda::reshaped_view(g, std::array<long, 2>{r, g.size() / r});
      auto gct  = nda::matrix<S>(transpose(g_rs));

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs (we transpose because LAPACK requires index into RHS # to be
      // slowest)

      if constexpr (nda::have_same_value_type_v<T, decltype(it2cf.lu)>) {
        nda::lapack::getrs(it2cf.lu, gct, it2cf.piv);
      } else {
        // NOTE: getrs require the first and second matrix to have the same value type
        // So if it2cf.lu does not have the same value type as gct we need a cast.
        nda::lapack::getrs(nda::matrix<S>(it2cf.lu), gct, it2cf.piv);
      }

      // Reshape to original dimensions and return
      auto gc = nda::matrix<S>(transpose(gct));
      return nda::reshaped_view(gc, g.shape());
    }

    /** 
    * @brief Transform DLR coefficients of Green's function G to values on DLR
    * imaginary time grid
    *
    * @param[in] gc DLR coefficients of G
    *
    * @return Values of G on DLR imaginary time grid
    * */
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
      requires(nda::is_scalar_v<S>)
    typename T::regular_type coefs2vals(T const &gc) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape gc to a matrix w/ first dimension r
      auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

      // Apply coeffs -> vals matrix
      auto g = cf2it * nda::matrix_const_view<S>(gc_rs);

      // [Q] Why reshaped_view instead of reshape? But this works, and reshape doesn't compile.
      // Reshape to original dimensions and return
      return nda::reshaped_view(g, gc.shape());
    }

    /** 
    * @brief Evaluate DLR expansion of G, given by its DLR coefficients, at imaginary
    * time point 
    *
    * @param[in] gc DLR coefficients of G
    * @param[in] t  Evaluation point
    *
    * @return Value of G at @p t
    */
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
      requires(nda::is_scalar_v<S>)
    auto coefs2eval(T const &gc, double t) const {

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
        auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

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
    **/
    nda::vector<double> build_evalvec(double t) const;

    /** 
    * @brief Compute convolution of two imaginary time Green's functions
    *
    * The convolution of f and g is defined as h(t) = (f * g)(t) int_0^beta
    * f(t-t') g(t') dt', where fermionic/bosonic antiperiodicity/periodicity
    * are used to define the Green's functions on (-beta, 0).  This method takes
    * the DLR coefficients of f and g as input and returns the values of h on
    * the DLR imaginary time grid.
    *
    * The convolution is computed on-the-fly in O(r^2) operations using Eq. (to
    * be added) in
    *
    * J. Kaye, H. U. R. Strand, D. Golez. Manuscript in preparation (2023).
    *
    * @param[in] beta Inverse temperature
    * @param[in] statistic Fermionic ("Fermion" or 0) or bosonic ("Boson" or 1)
    * @param[in] fc DLR coefficients of f
    * @param[in] gc DLR coefficients of g
    *
    * @return Values of h = f * g on DLR imaginary time grid
    * */
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
      requires(nda::is_scalar_v<S>)
    typename T::regular_type convolve(double beta, statistic_t statistic, T const &fc, T const &gc) {

      static constexpr auto _ = nda::range::all;

      if (r != fc.shape(0) || r != gc.shape(0)) throw std::runtime_error("First dim of input arrays must be equal to DLR rank r.");
      if (fc.shape() != gc.shape()) throw std::runtime_error("Input arrays must have the same shape.");

      // TODO: implement bosonic case and remove
      if (statistic == 0) throw std::runtime_error("imtime_ops::convolve not yet implemented for bosonic Green's functions.");

      // TODO: do this more cleanly
      // Initialize convolution, if it hasn't been done already
      if (hilb(0, 0) == -1.0) { convolve_init(); }

      if constexpr (T::rank == 1) { // Scalar-valued Green's function

        // Take array view of fc and gc
        auto fca = nda::array_const_view<S, 1>(fc);
        auto gca = nda::array_const_view<S, 1>(gc);

        // Diagonal contribution
        auto h = arraymult(convtmp, make_regular(fca * gca));

        // Off-diagonal contribution
        auto tmp = fca * arraymult(hilb, gca) + gca * arraymult(hilb, fca);
        return h + arraymult(cf2it, make_regular(tmp));

      } else if (T::rank == 3) { // Matrix-valued Green's function

        // Diagonal contribution
        auto fcgc = nda::array<S, 3>(fc.shape()); // Product of coefficients of f and g
        for (int i = 0; i < r; ++i) { fcgc(i, _, _) = arraymult(fc(i, _, _), gc(i, _, _)); }
        auto h = arraymult(convtmp, fcgc);

        // Off-diagonal contribution
        auto tmp1 = arraymult(hilb, fc);
        auto tmp2 = arraymult(hilb, gc);
        for (int i = 0; i < r; ++i) { tmp1(i, _, _) = arraymult(tmp1(i, _, _), gc(i, _, _)) + arraymult(fc(i, _, _), tmp2(i, _, _)); }

        return h + arraymult(cf2it, tmp1);

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

    /** Access DLR imaginary real frequency nodes*/
    /**
    * @brief Get DLR real frequency nodes
    *
    * @return DLR real frequency nodes
    */
    nda::vector_const_view<double> get_rfnodes() const { return dlr_rf; };

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

    private:
    /**
    * @brief Initialization for convolution methods
    *
    * This method pre-builds some matrices required for the convolution methods.
    */
    void convolve_init() {
      // "Discrete Hilbert transform" matrix -(1-delta_jk)/(dlr_rf(j) -
      // dlr_rf(k)), scaled by beta
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (j == k) {
            hilb(j, k) = 0;
          } else {
            hilb(j, k) = 1.0 / (dlr_rf(k) - dlr_rf(j));
          }
        }
      }

      // Matrix which applies DLR coefficients to imaginary time grid values
      // transformation matrix, and then multiplies the result by tau, the
      // imaginary time variable
      for (int j = 0; j < r; ++j) {
        for (int k = 0; k < r; ++k) {
          if (dlr_it(j) > 0) {
            tcf2it(j, k) = (dlr_it(j) - k_it(1.0, dlr_rf(k))) * cf2it(j, k);
          } else {
            tcf2it(j, k) = (dlr_it(j) + k_it(0.0, dlr_rf(k))) * cf2it(j, k);
          }
        }
      }
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
    nda::matrix<double> hilb;   ///< "Discrete Hilbert transform" matrix
    nda::matrix<double> tcf2it; ///< A matrix required for convolution

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
      auto cf2it     = h5::read<nda::matrix<double>>(gr, "cf2it");
      auto it2cf_lu  = h5::read<nda::matrix<double>>(gr, "it2cf_lu");
      auto it2cf_piv = h5::read<nda::vector<int>>(gr, "it2cf_piv");

      m = imtime_ops(lambda, rf, it, cf2it, it2cf_lu, it2cf_piv);
    }
  };

} // namespace cppdlr
