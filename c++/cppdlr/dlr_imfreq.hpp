// Copyright (c) 2023 Simons Foundation
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
#include "dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

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
    * @param[in] statistic Particle statistic: Fermion or Boson
    * @param[in] symmetrize NONSYM or false for non-symmetrized DLR frequencies,
    * SYM or true for symmetrized
    *
    * @note In case Boson and SYM options are selected, we enforce that i*nu_n=0 is
    * chosen as a DLR imaginary frequency node.
    */
    imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, statistic_t statistic, bool symmetrize = false);

    imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, statistic_t statistic, //
               nda::vector_const_view<int> dlr_if,                                          //
               nda::matrix_const_view<nda::dcomplex> cf2if,                                 //
               nda::matrix_const_view<nda::dcomplex> if2cf_lu,                              //
               nda::vector_const_view<int> if2cf_piv)
       : lambda_(lambda), statistic(statistic), r(dlr_rf.size()), dlr_rf(dlr_rf), dlr_if(dlr_if), cf2if(cf2if), if2cf{if2cf_lu, if2cf_piv} {};

    imfreq_ops() = default;

    /** 
    * @brief Transform values of Green's function G on DLR imaginary frequency grid to
    * DLR coefficients
    *
    * @param[in] beta Inverse temperature
    * @param[in] g Values of G on DLR imaginary frequency grid
    *
    * @return DLR coefficients of G
    */
    template <nda::MemoryArray T> typename T::regular_type vals2coefs(double beta, T const &g) const { return vals2coefs(make_regular(g / beta)); }

    template <nda::MemoryArray T> typename T::regular_type vals2coefs(T const &g) const {

      // MemoryArray type can be nda vector, matrix, array, or view of any of
      // these; taking a regular_type converts, for example, a matrix view to a
      // matrix.

      if (niom != g.shape(0)) throw std::runtime_error("First dim of g != # DLR imaginary frequency nodes.");

      // Make a copy of the data in Fortran Layout as required by getrs
      auto gf = nda::array<get_value_t<T>, get_rank<T>, F_layout>(g);

      // Reshape as matrix_view with r rows
      auto gfv = nda::reshape(gf, niom, g.size() / niom);

      // Solve linear system (multiple right hand sides) to convert vals ->
      // coeffs
      if (niom == r) {
        nda::lapack::getrs(if2cf.lu, gfv, if2cf.piv);
      } else {                                 // Non-square system---use least squares solver
        auto s       = nda::vector<double>(r); // Not needed
        double rcond = 0;                      // Not needed
        int rank     = 0;                      // Not needed
        nda::lapack::gelss(nda::matrix<dcomplex, F_layout>(cf2if), gfv, s, rcond, rank);
      }

      return gf(nda::range(r), nda::ellipsis());
    }

    /** 
    * @brief Transform DLR coefficients of Green's function G to values on DLR
    * imaginary frequency grid
    *
    * @param[in] beta Inverse temperature
    * @param[in] gc DLR coefficients of G
    *
    * @return Values of G on DLR imaginary frequency grid
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> make_cplx_t<T> coefs2vals(double beta, T const &gc) const {
      return coefs2vals(make_regular(beta * gc));
    }

    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> make_cplx_t<T> coefs2vals(T const &gc) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of gc != DLR rank r.");

      // Reshape gc to a matrix w/ first dimension r
      auto gc_rs = nda::reshape(gc, r, gc.size() / r);

      // Apply coeffs -> vals matrix
      auto g = cf2if * nda::matrix_const_view<S>(gc_rs);

      // Get output shape
      std::array<long, T::rank> shape_out;
      shape_out[0] = niom;
      for (int i = 1; i < T::rank; ++i) { shape_out[i] = gc.shape(i); }

      // Reshape to original dimensions and return
      return nda::reshape(g, shape_out);
    }

    /** 
    * @brief Evaluate DLR expansion of G, given by its DLR coefficients, at imaginary
    * frequency point 
    *
    * @param[in] beta Inverse temperature
    * @param[in] gc   DLR coefficients of G
    * @param[in] iom  Evaluation point
    *
    * @return Value of G at @p iom
    */
    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> auto coefs2eval(double beta, T const &gc, int n) const {
      return beta * coefs2eval(gc, n);
    }

    template <nda::MemoryArray T, nda::Scalar S = nda::get_value_t<T>> auto coefs2eval(T const &gc, int n) const {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions

      if constexpr (T::rank == 1) {

        // Evaluate DLR expansion
        std::complex<double> g = 0;
        for (int l = 0; l < r; ++l) { g += k_if(n, dlr_rf(l), statistic) * gc(l); }

        return g;
      } else {

        // Reshape gc to matrix w/ first dimension r
        auto gc_rs = nda::reshape(gc, r, gc.size() / r);

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
    * @param[in] beta Inverse temperature
    * @param[in] n  Evaluation point index
    *
    * @return Vector of evaluation at Matsubara frequency with index @p n
    **/
    nda::vector<nda::dcomplex> build_evalvec(double beta, int n) const;
    nda::vector<nda::dcomplex> build_evalvec(int n) const;

    /** 
    * @brief Get DLR imaginary frequency nodes
    *
    * @return DLR imaginary frequency nodes
    */
    nda::vector_const_view<int> get_ifnodes() const { return dlr_if; };
    int get_ifnodes(int i) const { return dlr_if(i); };

    /**
    * @brief Get DLR real frequency nodes
    *
    * @return DLR real frequency nodes
    */
    nda::vector_const_view<double> get_rfnodes() const { return dlr_rf; };
    double get_rfnodes(int i) const { return dlr_rf(i); };

    /**
    * @brief Get transformation matrix from DLR coefficients to values at DLR imaginary frequency nodes
    *
    * @return Transformation matrix
    */
    nda::matrix_const_view<nda::dcomplex> get_cf2if() const { return cf2if; };

    /**
    * @brief Get LU factors of transformation matrix from DLR imaginary frequency values to coefficients
    *
    * @return LU factors
    */
    nda::matrix_const_view<nda::dcomplex> get_if2cf_lu() const { return if2cf.lu; };

    /**
    * @brief Get LU pivots of transformation matrix from DLR imaginary frequency values to coefficients
    *
    * @return LU pivots
    */
    nda::vector_const_view<int> get_if2cf_piv() const { return if2cf.piv; };

    /** 
    * @brief Get DLR rank
    *
    * @return DLR rank
    */
    int rank() const { return r; }
    double lambda() const { return lambda_; }
    statistic_t get_statistic() const { return statistic; }

    private:
    double lambda_;                   ///< Energy cutoff divided by temperature
    statistic_t statistic;            ///< Particle statistic: Fermion or Boson
    int r;                            ///< DLR rank
    int niom;                         ///< # DLR imaginary freq nodes (different from r in symmetrized bosonic case)
    nda::vector<double> dlr_rf;       ///< DLR frequencies
    nda::vector<int> dlr_if;          ///< DLR imaginary frequency nodes
    nda::matrix<nda::dcomplex> cf2if; /// Transformation matrix from DLR coefficients to values at DLR imaginary frequency nodes

    /**
    * @brief Struct for transformation from DLR imaginary frequency values to coefficients
    */
    struct {
      nda::matrix<nda::dcomplex> lu; ///< LU factors (LAPACK format) of imaginary frequency vals -> coefs matrix
      nda::vector<int> piv;          ///< LU pivots (LAPACK format) of imaginary frequency vals -> coefs matrix
    } if2cf;

    // -------------------- hdf5 -------------------

    public:
    static std::string hdf5_format() { return "cppdlr::imfreq_ops"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, imfreq_ops const &m) {

      h5::group gr = fg.create_group(subgroup_name);
      write_hdf5_format_as_string(gr, "cppdlr::imfreq_ops");

      h5::write(gr, "lambda", m.lambda());
      h5::write<int>(gr, "statistic", m.get_statistic());
      h5::write(gr, "rf", m.get_rfnodes());
      h5::write(gr, "if", m.get_ifnodes());
      h5::write(gr, "cf2if", m.get_cf2if());
      h5::write(gr, "if2cf_lu", m.get_if2cf_lu());
      h5::write(gr, "if2cf_piv", m.get_if2cf_piv());
    }

    friend void h5_read(h5::group fg, std::string const &subgroup_name, imfreq_ops &m) {

      h5::group gr = fg.open_group(subgroup_name);
      assert_hdf5_format_as_string(gr, "cppdlr::imfreq_ops", true);

      double lambda;
      statistic_t statistic_;
      nda::vector<double> rf;
      nda::vector<int> if_;
      nda::matrix<nda::dcomplex> cf2if_;
      nda::matrix<nda::dcomplex> if2cf_lu;
      nda::vector<int> if2cf_piv;

      h5::read(gr, "lambda", lambda);
      statistic_ = statistic_t(h5::read<int>(gr, "statistic"));
      h5::read(gr, "rf", rf);
      h5::read(gr, "if", if_);
      h5::read(gr, "cf2if", cf2if_);
      h5::read(gr, "if2cf_lu", if2cf_lu);
      h5::read(gr, "if2cf_piv", if2cf_piv);

      m = imfreq_ops(lambda, rf, statistic_, if_, cf2if_, if2cf_lu, if2cf_piv);
    }
  };

} // namespace cppdlr
