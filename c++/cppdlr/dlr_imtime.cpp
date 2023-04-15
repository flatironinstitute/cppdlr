#include "dlr_imtime.hpp"
#include "cppdlr/dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

using namespace nda;

namespace cppdlr {

  imtime_ops::imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf) : lambda_(lambda), r(dlr_rf.size()), dlr_rf(dlr_rf) {

    dlr_it    = nda::vector<double>(r);
    cf2it     = nda::matrix<double>(r, r);
    it2cf.lu  = nda::matrix<double>(r, r);
    it2cf.piv = nda::vector<int>(r);

    // Get discretization of analytic continuation kernel on fine grid in
    // imaginary time, at DLR frequencies
    auto fine = fineparams(lambda);
    auto t    = build_it_fine(fine);
    auto kmat = build_k_it(t, dlr_rf);

    // Pivoted Gram-Schmidt to obtain DLR imaginary time nodes
    auto [q, norms, piv] = pivrgs(kmat, 1e-100);
    for (int i = 0; i < r; ++i) { dlr_it(i) = t(piv(i)); }

    // Obtain coefficients to imaginary time values transformation matrix
    for (int i = 0; i < r; ++i) {
      for (int j = 0; j < r; ++j) { cf2it(i, j) = kmat(piv(i), j); }
    }

    // Prepare imaginary time values to coefficients transformation by computing
    // LU factors of coefficient to imaginary time matrix
    it2cf.lu = cf2it;
    lapack::getrf(it2cf.lu, it2cf.piv);
  }

  nda::vector<double> imtime_ops::build_evalvec(double t) const {

    // TODO: can be further optimized to reduce # exponential evals.
    auto kvec = nda::vector<double>(r);
    if (t >= 0) {
      for (int l = 0; l < r; ++l) { kvec(l) = kfun_abs(t, dlr_rf(l)); }
    } else {
      for (int l = 0; l < r; ++l) { kvec(l) = kfun_abs(-t, -dlr_rf(l)); }
    }

    return kvec;
  }

} // namespace cppdlr
