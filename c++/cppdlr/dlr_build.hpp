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

    // [Q] I want to make all of these const, so that they are set in the
    // constructor once and then can't be modified. It isn't really essential,
    // because this is not a user-facing class, but it seems technically
    // correct. But it seems that if I do this, I'm forced to set the values of
    // these variables in the first line of the constructor, which is awkward
    // and ugly. What's the right way to do this?
    
    double lambda; /// DLR cutoff
    int p; /// Fine grid Chebyshev order
    int npt; /// Number fine imaginary time grid panels
    int npo; /// Number fine frequency grid panels
    int nt; /// Total number fine imaginary time grid points
    int no; /// Total number fine frequency grid points

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
