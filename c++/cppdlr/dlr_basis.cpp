#include "dlr_basis.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

namespace cppdlr {

  nda::vector<double> dlr_freq(double lambda, double eps) {

    // Get fine grid parameters

    auto fine = fineparams(lambda);

    // Get fine grids in frequency and imaginary time

    auto t = get_tfine(fine);
    auto om = get_omfine(fine);

    // Get discretization of analytic continuation kernel on fine grid (the K matrix)

    auto kmat = get_kfine(fine, t, om);

    // Pivoted Gram-Schmidt on columns of K matrix to obtain DLR frequencies

    auto [q, norms, piv] = pivrgs(transpose(kmat), eps);
    int r                = norms.size();

    auto omega = nda::vector<double>(r);
    for (int i = 0; i < r; ++i) { omega(i) = om(piv(i)); }

    return omega;
  }

  // nda::vector<double> dlr_it(nda::vector<double> dlr_rf) {



  // }

} // namespace cppdlr