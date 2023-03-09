#pragma once
#include <nda/nda.hpp>
#include "cppdlr/dlr_kernels.hpp"

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

      // Reshape g to matrix w/ second dimension r
      auto g_rs = nda::reshaped_view(g, std::array<long, 2>{r, g.size() / r});
      auto gct  = nda::matrix<S>(transpose(g_rs));

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs (we transpose because LAPACK requires index into RHS # to be
      // slowest)
      nda::lapack::getrs(it2cf.lu, gct, it2cf.piv);

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

      // Reshape gc to a matrix w/ second dimension r
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
          for (int l = 0; l < r; ++l) { g += kfun_abs(t, dlr_rf(l)) * gc(l); }
        } else {
          for (int l = 0; l < r; ++l) { g += kfun_abs(-t, -dlr_rf(l)) * gc(l); }
        }

        return g;
      } else {

        // Reshape gc to matrix w/ second dimension r
        auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

        // Get output shape
        std::array<long, T::rank - 1> shape_out;
        for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i + 1); }

        // Get vector of evaluation of DLR expansion at a point
        auto kvec = get_kevalvec(t);

        // Evaluate DLR expansion, reshape to original dimensions (with first
        // dimension summed out), and return
        return nda::reshape(transpose(nda::matrix_const_view<S>(gc_rs)) * kvec, shape_out);
      }
    }

    /** 
    * @brief Get vector of evaluation of DLR expansion at an imaginary time point 
    *
    * @param[in] t  Evaluation point
    *
    * @return Vector of evaluation at @p t
    **/
    nda::vector<double> get_kevalvec(double t) const;

    /** 
    * @brief Get DLR imaginary time nodes
    *
    * @return DLR imaginary time nodes
    */
    nda::vector_const_view<double> get_itnodes() const { return dlr_it; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }

    private:
    int r;                      ///< DLR rank
    nda::vector<double> dlr_rf; ///< DLR frequencies
    nda::vector<double> dlr_it; ///< DLR imaginary time nodes
    nda::matrix<double> cf2it;  ///< Transformation matrix from DLR coefficients to values at DLR imaginary time nodes

    /**
    * Struct for transformation from DLR imaginary time values to coefficients
    */
    struct {
      nda::matrix<double> lu;   ///< LU factors (LAPACK format) of imaginary time vals -> coefs matrix
      nda::vector<int> piv;     ///< LU pivots (LAPACK format) of imaginary time vals -> coefs matrix
    } it2cf; 
  };

} // namespace cppdlr
