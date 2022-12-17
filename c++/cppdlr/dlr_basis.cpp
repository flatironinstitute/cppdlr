#include "dlr_basis.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

namespace cppdlr {

  dlr_basis::dlr_basis(double lambda, double eps) {

    // Get fine grid parameters

    fineparams fine(lambda);

    // Get fine grids in frequency and imaginary time

    auto [t, om] = get_finegrids(fine);

    // Get discretization of analytic continuation kernel on fine grid (the K matrix)

    auto kmat = get_kfine(fine, t, om);

    // Pivoted Gram-Schmidt on columns of K matrix to obtain DLR frequencies

    auto [q, norms, piv] = pivrgs(kmat, eps);
    r                    = norms.size();

    auto omega = nda::vector<double>(r);
    // [Q] Can this be done more concisely?
    for (int i = 0; i < r; ++i) { omega(i) = om(piv(i)); }
  }

} // namespace cppdlr