#include "dlr_imtime.hpp"
#include "cppdlr/dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

using namespace nda;

static constexpr auto _  = range::all;
static constexpr auto __ = nda::ellipsis();

namespace cppdlr {

  imtime_ops::imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf) : r(dlr_rf.size()), dlr_rf(dlr_rf) {

    dlr_it    = nda::vector<double>(r);
    cf2it     = nda::matrix<double>(r, r);
    it2cf.lu  = nda::matrix<double>(r, r);
    it2cf.piv = nda::vector<int>(r);

    // Get discretization of analytic continuation kernel on fine grid in
    // imaginary time, at DLR frequencies

    auto fine = fineparams(lambda);
    auto t    = get_tfine(fine);
    auto kmat = get_kfine(t, dlr_rf);

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

  nda::matrix<double> imtime_ops::vals2coefs_mat(nda::matrix_const_view<double> g) {

    auto gc = nda::matrix<double>(g);
    lapack::getrs(it2cf.lu, gc, it2cf.piv);

    return gc;
  }

  nda::matrix<double> imtime_ops::coefs2vals_mat(nda::matrix_const_view<double> gc) {

    return transpose(cf2it * transpose(gc));

  }


  nda::vector<double> imtime_ops::coefs2eval_mat(nda::matrix_const_view<double> gc, double t) {

    // TODO: this can be further optimized; for example, can reduce # of exponential evaluations.
    auto kvec = nda::vector<double>(r);
    if (t >= 0) {
      for (int l = 0; l < r; ++l) {
        kvec(l) = kfun_abs(t, dlr_rf(l));
      }
    } else {
      for (int l = 0; l < r; ++l) {
        kvec(l) = kfun_abs(-t, -dlr_rf(l));
      }
    }

    return gc*kvec;
  }

  double imtime_ops::coefs2eval_vec(nda::vector_const_view<double> gc, double t) {

    // TODO: this can be further optimized; for example, can reduce # of exponential evaluations.
    double g = 0;
    if (t >= 0) {
      for (int l = 0; l < r; ++l) {
        g += kfun_abs(t, dlr_rf(l))*gc(l);
      }
    } else {
      for (int l = 0; l < r; ++l) {
        g += kfun_abs(-t, -dlr_rf(l))*gc(l);
      }
    }

    return g;
  }

} // namespace cppdlr
