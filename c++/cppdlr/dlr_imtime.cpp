#include "dlr_imtime.hpp"
#include "cppdlr/dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>

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

    // [Q] Need to transpose first in order to get expected behavior...this should be documented clearly
    it2cf.lu = transpose(cf2it);

    // [Q] How to get rid of compiler warnings about unused info variable?

    int info = lapack::getrf(it2cf.lu, it2cf.piv);

  }

  nda::array<double, 3> imtime_ops::vals2coefs(nda::array_const_view<double, 3> g) {

    // [Q] This is very ugly. Can this be improved? Need to jump through these
    // hoops because (1) getrs must take in rank 2 array, and (2) getrs is using
    // Fortran-ordering instead of C-ordering. Want to be able to just call
    // lapack::getrs(it2cf.lu, g(__,_), it2cf.piv) and have it do the right
    // thing, i.e., reshape as an norb^2 x r array, and do norb^2 rxr solves.

    int norb = g.shape(0);
    auto gc  = nda::array<double, 3>(norb, norb, r);

    auto tmp = nda::matrix<double>(1, r);
    int info = 0;
    for (int i = 0; i < norb; ++i) {
      for (int j = 0; j < norb; ++j) {
        tmp(0, _)   = g(i, j, _);
        info        = lapack::getrs(it2cf.lu, tmp, it2cf.piv);
        gc(i, j, _) = tmp(0, _);
      }
    }

    return gc;
  }

  nda::array<double, 3> imtime_ops::coefs2vals(nda::array_const_view<double, 3> gc) {
    int norb = gc.shape(0);
    auto g   = nda::array<double, 3>(norb, norb, gc.shape(2));

    // [Q] Can I do this in LAPACK style with one line?
    for (int i = 0; i < norb; ++i) {
      for (int j = 0; j < norb; ++j) { g(i, j, _) = cf2it * nda::vector<double>(gc(i, j, _)); }
    }
    return g;
  }

  nda::matrix<double> imtime_ops::coefs2eval(nda::array_const_view<double, 3> gc, double t) {

    int norb = gc.shape(0);

    auto g = nda::matrix<double>(norb, norb);
    g = 0;

    double kval = 0;
    if (t >= 0) {
      for (int l = 0; l < r; ++l) {
        kval = kfun_abs(t, dlr_rf(l));
        g += kval * gc(_, _, l);
      }
    } else {
      for (int l = 0; l < r; ++l) {
        kval = kfun_abs(-t, -dlr_rf(l));
        g += kval * gc(_, _, l);
      }
    }

    return g;
  }

  nda::vector_const_view<double> imtime_ops::get_itnodes() const {return dlr_it;}

} // namespace cppdlr