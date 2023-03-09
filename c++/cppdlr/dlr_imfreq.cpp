#include "dlr_imfreq.hpp"
#include "cppdlr/dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

using namespace nda;

namespace cppdlr {

  imfreq_ops::imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, int xi) : xi(xi), r(dlr_rf.size()), dlr_rf(dlr_rf) {

    if (xi != 1 && xi != -1) throw std::runtime_error("xi must be -1 (fermionic) or 1 (bosonic).");

    dlr_if    = nda::vector<int>(r);
    cf2if     = nda::matrix<dcomplex>(r, r);
    if2cf.lu  = nda::matrix<dcomplex>(r, r);
    if2cf.piv = nda::vector<int>(r);

    // Get analytic continuation kernel at DLR frequencies, up to imaginary
    // frequency cutoff
    auto nmax = fineparams(lambda).nmax;
    auto kmat = get_kif(nmax, dlr_rf, xi);

    // Pivoted Gram-Schmidt to obtain DLR imaginary time nodes
    auto [q, norms, piv] = pivrgs(kmat, 1e-100);
    for (int i = 0; i < r; ++i) { dlr_if(i) = piv(i) - nmax; }

    // Obtain coefficients to imaginary frequency values transformation matrix
    for (int i = 0; i < r; ++i) {
      for (int j = 0; j < r; ++j) { cf2if(i, j) = kmat(piv(i), j); }
    }

    // Prepare imaginary time values to coefficients transformation by computing
    // LU factors of coefficient to imaginary time matrix
    if2cf.lu = cf2if;
    lapack::getrf(if2cf.lu, if2cf.piv);
  }

  nda::vector<dcomplex> imfreq_ops::get_kevalvec(int n) const {

    auto kvec = nda::vector<dcomplex>(r);
    for (int l = 0; l < r; ++l) { kvec(l) = kfun_if(2*n+(1-xi)/2, dlr_rf(l)); }

    return kvec;
  }

} // namespace cppdlr
