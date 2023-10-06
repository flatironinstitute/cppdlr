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
// Authors: Nils Wentzell, Jason Kaye

#pragma once
#include <nda/nda.hpp>
#include "dlr_kernels.hpp"

using namespace nda;

namespace cppdlr {

  /** 
  * @class fineparams
  * @brief Class containing parameters for fine composite Chebyshev grid
  * discretizations of imaginary time and frequency
  *
  * @param[in] lambda DLR cutoff parameter
  * @param[in] p Order of composite Chebyshev grid discretization
  *
  * \note Values have been chosen empirically to give fine discretization of
  * Lehmann kernel accurate to double machine precision
  */
  class fineparams {

    public:
    fineparams(double lambda, int p = 24);

    const double lambda; ///< DLR cutoff parameter
    const int p;         ///< Order of composite Chebyshev grid
    const int nmax;      ///< Imaginary frequency cutoff
    const int npom;      ///< # fine real frequency grid panels
    const int npt;       ///< # fine imaginary time grid panels
    const int nom;       ///< Total # fine real frequency grid points
    const int nt;        ///< Total # fine imaginary time grid points
  };

  /**
  * @brief Build fine composite Chebyshev grid in real frequency 
  *
  * @param[in] fine Fine grid parameters
  *
  * @return Fine real frequency grid
  */
  nda::vector<double> build_rf_fine(fineparams &fine);

  /**
  * @brief Get fine composite Legendre grid in imaginary time and corresponding
  * square root quadrature weights
  *
  * @param[in] fine Fine grid parameters
  *
  * @return Tuple containing fine imaginary time grid and corresponding square
  * root quadrature weights
  *
  * \note Fine imaginary time grid is given in relative format 
  */
  std::tuple<nda::vector<double>, nda::vector<double>> build_it_fine(fineparams &fine);

  /**
  * @brief Get imaginary time discretization of analytic continuation kernel
  *
  * @param[in] t Imaginary time grid in relative format
  * @param[in] om Frequency real grid
  *
  * @return Discretization of analytic continuation kernel on given grid
  */
  nda::matrix<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> om);

  /**
  * @brief Get imaginary time discretization of analytic continuation kernel
  * with L^2 weighting
  *
  * @param[in] t Imaginary time grid in relative format
  * @param[in] w Square root quadrature weights
  * @param[in] om Frequency real grid
  *
  * @return Discretization of analytic continuation kernel on given grid
  */
  nda::matrix<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> w, nda::vector_const_view<double> om);

  /** 
  * @copydoc build_k_it(nda::vector_const_view<double>, nda::vector_const_view<double>)
  */
  nda::vector<double> build_k_it(double t, nda::vector_const_view<double> om);

  /** 
  * @copydoc build_k_it(nda::vector_const_view<double>, nda::vector_const_view<double>)
  */
  nda::vector<double> build_k_it(nda::vector_const_view<double> t, double om);

  /** 
  * @copydoc build_k_it(nda::vector_const_view<double>, nda::vector_const_view<double>, nda::vector_const_view<double>)
  */
  nda::vector<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> w, double om);

  /**
  * @brief Get imaginary frequency discretization of analytic continuation kernel
  *
  * @param[in] nmax Imaginary frequency cutoff
  * @param[in] om Real frequency grid
  * @param[in] statistic Particle Statistic: Boson or Fermion
  *
  * @return Discretization of analytic continuation kernel on given grid
  */
  nda::matrix<dcomplex> build_k_if(int nmax, nda::vector_const_view<double> om, statistic_t statistic);

  /**
  * @brief Get error of fine composite Chebyshev discretization of analytic
  * continuation kernel in imaginary time
  * 
  * @param[in] fine Fine grid parameters
  * @param[in] t Imaginary time grid in relative format
  * @param[in] om Real frequency grid
  * @param[in] kmat Discretization of analytic continuation kernel on given grid
  *
  * @return Error of given fine discretization of analytic continuation kernel
  * in imaginary time and real frequency
  *
  * \note Error is given as an estimate of the maximum absolute value of the
  * difference between the given discretization and the exact analytic
  * continuation kernel
  *
  * \note \p kmat should be computed using the function get_kfine with composite
  * Chebyshev grids produced by get_omfine and get_tfine
  */
  std::tuple<double, double> geterr_k_it(fineparams &fine, nda::vector_const_view<double> t, nda::vector_const_view<double> om,
                                         nda::matrix_const_view<double> kmat);

  /**
  * Option for unsymmetrized or symmetrized DLR
  */
  enum { NONSYM = false, SYM = true };

  /**
  * @brief Construct DLR basis by obtaining DLR frequencies
  *
  * @param[in] lambda DLR cutoff parameter
  * @param[in] eps Accuracy of DLR basis
  * @param[in] statistic Particle statistics: Boson or Fermion
  * @param[in] symmetrize NONSYM or false for non-symmetrized DLR frequencies,
  * SYM or true for symmetrized
  *
  * @return DLR frequencies
  *
  * @note In case Boson and SYM options are selected, we enforce that omega = 0
  * is chosen as a DLR frequency. For an explanation, please see the
  * "Symmetrized DLR grids" subsection in the "Background" section of the
  * documentation.
  */
  nda::vector<double> build_dlr_rf(double lambda, double eps, statistic_t statistic, bool symmetrize);

  /**
  * @brief Construct DLR basis by obtaining DLR frequencies
  *
  * @param[in] lambda DLR cutoff parameter
  * @param[in] eps Accuracy of DLR basis
  *
  * @return DLR frequencies
  *
  * @note With this signature, the DLR frequencies are unsymmetrized.
  */
  nda::vector<double> build_dlr_rf(double lambda, double eps);

} // namespace cppdlr
