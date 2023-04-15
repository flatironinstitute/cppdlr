#pragma once
#include <nda/nda.hpp>

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
    const int npom;      ///< # fine frequency grid panels
    const int npt;       ///< # fine imaginary time grid panels
    const int nom;       ///< Total # fine frequency grid points
    const int nt;        ///< Total # fine imaginary time grid points
  };

  /**
  * @brief Get fine composite Chebyshev grid in frequency 
  *
  * @param[in] fine Fine grid parameters
  *
  * @return Fine frequency grid
  */
  nda::vector<double> get_omfine(fineparams &fine);

  /**
  * @brief Get fine composite Chebyshev grid in imaginary time
  *
  * @param[in] fine Fine grid parameters
  *
  * @return Fine imaginary time grid
  *
  * \note Fine imaginary time grid is given in relative format 
  */
  nda::vector<double> get_tfine(fineparams &fine);

  /**
  * @brief Get imaginary time discretization of analytic continuation kernel
  *
  * @param[in] t Imaginary time grid in relative format
  * @param[in] om Frequency grid
  *
  * @return Discretization of analytic continuation kernel on given grid
  */
  nda::matrix<double> get_kfine(nda::vector_const_view<double> t, nda::vector_const_view<double> om);

  /**
  * @brief Get imaginary frequency discretization of analytic continuation kernel
  *
  * @param[in] nmax Imaginary frequency cutoff
  * @param[in] om Frequency grid
  * @param[in] xi Fermionic (xi = -1) or bosonic (xi = 1) imaginary frequency grid
  *
  * @return Discretization of analytic continuation kernel on given grid
  */
  nda::matrix<dcomplex> get_kif(int nmax, nda::vector_const_view<double> om, int xi);

  /**
  * @brief Get error of fine composite Chebyshev discretization of analytic
  * continuation kernel in imaginary time
  * 
  * @param[in] fine Fine grid parameters
  * @param[in] t Imaginary time grid in relative format
  * @param[in] om Frequency grid
  * @param[in] kmat Discretization of analytic continuation kernel on given grid
  *
  * @return Error of given fine discretization of analytic continuation kernel
  * in imaginary time and in frequency
  *
  * \note Error is given as an estimate of the maximum absolute value of the
  * difference between the given discretization and the exact analytic
  * continuation kernel
  *
  * \note \p kmat should be computed using the function get_kfine with composite
  * Chebyshev grids produced by get_omfine and get_tfine
  */
  std::tuple<double, double> get_kfineerr(fineparams &fine, nda::vector_const_view<double> t, nda::vector_const_view<double> om,
                                          nda::matrix_const_view<double> kmat);

  /**
   * Construct DLR basis for a given accuracy and cutoff parameter by getting DLR frequencies
   */

  /**
  * @brief Construct DLR basis by obtaining DLR frequencies
  *
  * @param[in] lambda DLR cutoff parameter
  * @param[in] eps Accuracy of DLR basis
  *
  * @return DLR frequencies
  */
  nda::vector<double> dlr_freq(double lambda, double eps);

} // namespace cppdlr
