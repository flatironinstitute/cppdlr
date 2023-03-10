#pragma once
#include <nda/nda.hpp>
#include "cppdlr/dlr_kernels.hpp"
#include "cppdlr/utils.hpp"

namespace cppdlr {

  /**
  * @class imfreq_ops
  * @brief Class responsible for all DLR imaginary frequency operations, including
  * building imaginary frequency grid and transformations.
  *
  * \note First dimension of all Green's function and coefficient arrays must be
  * DLR rank r.
  */

  class imfreq_ops {

    public:
    /** 
    * @brief Constructor for imfreq_ops
    * 
    * @param[in] lambda DLR cutoff parameter
    * @param[in] dlr_rf DLR frequencies
    * @param[in] xi     Fermionic (xi = -1) or bosonic (xi = 1) Matsubara frequencies
    */
    imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, int xi);

    /** 
    * @brief Transform values of Green's function G on DLR imaginary frequency grid to
    * real-valued DLR coefficients
    *
    * Use this vals2coefs method when the DLR coefficients of G should be
    * real-valued; in other words, when G is real-valued in imaginary time 
    *
    * @param[in] g Values of G on DLR imaginary frequency grid
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T> make_real_t<T> vals2realcoefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape g to matrix w/ second dimension r
      auto g_rs = nda::reshaped_view(g, std::array<long, 2>{r, g.size() / r});
      auto gct  = nda::matrix<nda::dcomplex>(transpose(g_rs));

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs (we transpose because LAPACK requires index into RHS # to be
      // slowest)
      nda::lapack::getrs(if2cf.lu, gct, if2cf.piv);

      // Take real part (imaginary part of DLR coefficients should be zero),
      // reshape to original dimensions and return
      auto gc = nda::matrix<double>(transpose(real(gct)));
      return nda::reshaped_view(gc, g.shape());
    }

    /** 
    * @brief Transform values of Green's function G on DLR imaginary frequency grid to
    * complex-valued DLR coefficients
    *
    * Use this vals2coefs method when the DLR coefficients of G should be
    * complex-valued; in other words, when G is complex-valued in imaginary time
    *
    * @param[in] g Values of G on DLR imaginary frequency grid
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T> typename T::regular_type vals2cplxcoefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape g to matrix w/ second dimension r
      auto g_rs = nda::reshaped_view(g, std::array<long, 2>{r, g.size() / r});
      auto gct  = nda::matrix<nda::dcomplex>(transpose(g_rs));

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs (we transpose because LAPACK requires index into RHS # to be
      // slowest)
      nda::lapack::getrs(if2cf.lu, gct, if2cf.piv);

      // Reshape to original dimensions and return
      auto gc = nda::matrix<nda::dcomplex>(transpose(gct));
      return nda::reshaped_view(gc, g.shape());
    }

    /** 
* @brief Transform DLR coefficients of Green's function G to values on DLR
* imaginary frequency grid
*
* @param[in] gc DLR coefficients of G
*
* @return Values of G on DLR imaginary frequency grid
* */
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
    requires(nda::is_scalar_v<S>) make_cplx_t<T> coefs2vals(T const &gc)
    const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Reshape gc to a matrix w/ second dimension r
      auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

      // Apply coeffs -> vals matrix
      auto g = cf2if * nda::matrix_const_view<S>(gc_rs);

      // [Q] Why reshaped_view instead of reshape? But this works, and reshape doesn't compile.
      // Reshape to original dimensions and return
      return nda::reshaped_view(g, gc.shape());
    }

    /** 
    * @brief Evaluate DLR expansion of G, given by its DLR coefficients, at imaginary
    * frequency point 
    *
    * @param[in] gc   DLR coefficients of G
    * @param[in] iom  Evaluation point
    *
    * @return Value of G at @p iom
    */
    template <nda::MemoryArray T, typename S = nda::get_value_t<T>>
    requires(nda::is_scalar_v<S>) auto coefs2eval(T const &gc, int n) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions

      if constexpr (T::rank == 1) {

        // Evaluate DLR expansion
        std::complex<double> g = 0;
        for (int l = 0; l < r; ++l) { g += kfun_if(2*n+(1-xi)/2, dlr_rf(l)) * gc(l); }

        return g;
      } else {

        // Reshape gc to matrix w/ second dimension r
        auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

        // Get output shape
        std::array<long, T::rank - 1> shape_out;
        for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i + 1); }

        // Get vector of evaluation of DLR expansion at a point
        auto kvec = get_kevalvec(n);

        // Evaluate DLR expansion, reshape to original dimensions (with first
        // dimension summed out), and return
        return nda::reshape(transpose(nda::matrix_const_view<S>(gc_rs)) * kvec, shape_out);
      }
    }

    /** 
    * @brief Get vector of evaluation of DLR expansion at an imaginary frequency point 
    *
    * @param[in] n  Evaluation point index
    *
    * @return Vector of evaluation at Matsubara frequency with index @p n
    **/
    nda::vector<nda::dcomplex> get_kevalvec(int n) const;

    /** 
    * @brief Get DLR imaginary frequency nodes
    *
    * @return DLR imaginary frequency nodes
    */
    nda::vector_const_view<int> get_ifnodes() const { return dlr_if; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }

    private:
    int xi;                           ///< Fermionic (xi = -1) or bosonic (xi = 1) Matsubara frequencies
    int r;                            ///< DLR rank
    nda::vector<double> dlr_rf;       ///< DLR frequencies
    nda::vector<int> dlr_if;          ///< DLR imaginary frequency nodes
    nda::matrix<nda::dcomplex> cf2if; /// Transformation matrix from DLR coefficients to values at DLR imaginary frequency nodes

    /**
    * Struct for transformation from DLR imaginary frequency values to coefficients
    */
    struct {
      nda::matrix<nda::dcomplex> lu; ///< LU factors (LAPACK format) of imaginary frequency vals -> coefs matrix
      nda::vector<int> piv;          ///< LU pivots (LAPACK format) of imaginary frequency vals -> coefs matrix
    } if2cf;
  };

} // namespace cppdlr
