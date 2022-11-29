#pragma once
#include <vector>
#include <nda/nda.hpp>

namespace cppdlr {

  /**
   * Contains empirically-chosen fine grid parameters to use in
   * construction of DLR basis.
   */

  class fineparams {

    public:
    fineparams(double lambda);

    const double lambda; /// DLR cutoff
    const int p; /// Fine grid Chebyshev order
    const int npt; /// Number fine imaginary time grid panels
    const int npo; /// Number fine frequency grid panels
    const int nt; /// Total number fine imaginary time grid points
    const int no; /// Total number fine frequency grid points

  };

  /**
   * Get fine grids (composite Chebyshev grids) in frequency and imaginary time.
   */

  std::tuple<nda::vector<double>, nda::vector<double>>
    get_finegrids(fineparams &fine);


  /**
   * Get fine discretization of Lehmann kernel.
   */

  nda::matrix<double> get_kfine(fineparams &fine, nda::vector<double> &t,
      nda::vector<double> &om);
  
} // namespace cppdlr
