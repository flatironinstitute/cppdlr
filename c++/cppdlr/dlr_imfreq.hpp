#pragma once
#include <nda/nda.hpp>
#include "cppdlr/dlr_kernels.hpp"
#include "cppdlr/utils.hpp"

#include <h5/h5.hpp>
#include <nda/h5.hpp>

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

    imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, int xi, nda::vector_const_view<int> dlr_if,
               nda::matrix_const_view<nda::dcomplex> cf2if, nda::matrix_const_view<nda::dcomplex> if2cf_lu, nda::vector_const_view<int> if2cf_piv)
       : lambda_(lambda), xi(xi), r(dlr_rf.size()), dlr_rf(dlr_rf), dlr_if(dlr_if), cf2if(cf2if), if2cf{if2cf_lu, if2cf_piv} {};

    imfreq_ops() = default;

    /** 
    * @brief Transform values of Green's function G on DLR imaginary frequency grid to
    * DLR coefficients
    *
    * @param[in] g Values of G on DLR imaginary frequency grid
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T> typename T::regular_type vals2coefs(T const &g) const {

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
      requires(nda::is_scalar_v<S>)
    make_cplx_t<T> coefs2vals(T const &gc) const {

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
      requires(nda::is_scalar_v<S>)
    auto coefs2eval(T const &gc, int n) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions

      if constexpr (T::rank == 1) {

        // Evaluate DLR expansion
        std::complex<double> g = 0;
        for (int l = 0; l < r; ++l) { g += kfun_if(2 * n + (1 - xi) / 2, dlr_rf(l)) * gc(l); }

        return g;
      } else {

        // Reshape gc to matrix w/ second dimension r
        auto gc_rs = nda::reshaped_view(gc, std::array<long, 2>{r, gc.size() / r});

        // Get output shape
        std::array<long, T::rank - 1> shape_out;
        for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i + 1); }

        // Get vector of evaluation of DLR expansion at a point
        auto kvec = build_evalvec(n);

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
    nda::vector<nda::dcomplex> build_evalvec(int n) const;

    /** 
    * @brief Get DLR imaginary frequency nodes
    *
    * @return DLR imaginary frequency nodes
    */
    nda::vector_const_view<int> get_ifnodes() const { return dlr_if; };

    /** Accessors */
    nda::vector_const_view<double> get_rfnodes() const { return dlr_rf; };
    nda::matrix_const_view<nda::dcomplex> get_cf2if() const { return cf2if; };
    nda::matrix_const_view<nda::dcomplex> get_if2cf_lu() const { return if2cf.lu; };
    nda::vector_const_view<int> get_if2cf_piv() const { return if2cf.piv; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }
    double lambda() const { return lambda_; }
    int get_xi() const { return xi; }

    private:
    double lambda_;
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

    // -------------------- hdf5 -------------------

    static std::string hdf5_format() { return "cppdlr::imfreq_ops"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, imfreq_ops const &m) {

      h5::group gr = fg.create_group(subgroup_name);
      write_hdf5_format_as_string(gr, "cppdlr::imfreq_ops");

      h5_write(gr, "lambda", m.lambda());
      h5_write(gr, "xi", m.get_xi());
      h5_write(gr, "rf", m.get_rfnodes());
      h5_write(gr, "if", m.get_ifnodes());
      h5_write(gr, "cf2if", m.get_cf2if());
      h5_write(gr, "if2cf_lu", m.get_if2cf_lu());
      h5_write(gr, "if2cf_piv", m.get_if2cf_piv());
    }

    friend void h5_read(h5::group fg, std::string const &subgroup_name, imfreq_ops &m) {

      h5::group gr = fg.open_group(subgroup_name);
      assert_hdf5_format_as_string(gr, "cppdlr::imfreq_ops", true);

      double lambda;
      int xi;
      nda::vector<double> rf;
      nda::vector<int> if_;
      nda::matrix<nda::dcomplex> cf2if;
      nda::matrix<nda::dcomplex> if2cf_lu;
      nda::vector<int> if2cf_piv;

      h5_read(gr, "lambda", lambda);
      h5_read(gr, "xi", xi);
      h5_read(gr, "rf", rf);
      h5_read(gr, "if", if_);
      h5_read(gr, "cf2if", cf2if);
      h5_read(gr, "if2cf_lu", if2cf_lu);
      h5_read(gr, "if2cf_piv", if2cf_piv);

      m = imfreq_ops(lambda, rf, xi, if_, cf2if, if2cf_lu, if2cf_piv);
    }
  };

} // namespace cppdlr
