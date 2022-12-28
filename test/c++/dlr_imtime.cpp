#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_basis.hpp>
#include <cppdlr/utils.hpp>
#include <cppdlr/dlr_kernels.hpp>

using namespace cppdlr;

double gfun(double beta, double t) {

  int n  = 5;
  auto a = nda::vector<double>(n);

  a(0) = -0.804;
  a(1) = -0.443;
  a(2) = 0.093;
  a(3) = 0.915;
  a(4) = 0.929;

  double g = 0;
  for (int i = 0; i < n; ++i) { g += kfun(t, beta * a(i)); }
  return g;
}

TEST(dlr_imtime, interp) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  int norb = 1; // Scalar-valued Green's function

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = dlr_rf.size();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();

  auto g = nda::array<double, 3>(norb, norb, r);

  for (int i = 0; i < r; ++i) { g(0, 0, i) = gfun(beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes

  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-15);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error

  auto gtru  = nda::array<double, 3>(norb, norb, ntst);
  auto gtst  = nda::array<double, 3>(norb, norb, ntst);
  double err = 0;
  auto tmp   = nda::matrix<double>(norb, norb);
  for (int i = 0; i < ntst; ++i) {
    gtru(0, 0, i) = gfun(beta, ttst(i));

    // [Q] Best way to clean this up?
    tmp           = itops.coefs2eval(gc, ttst(i));
    gtst(0, 0, i) = tmp(0, 0);

    err = std::max(err, abs(gtru(0, 0, i) - gtst(0, 0, i)));
  }

  EXPECT_LT(err, 10 * eps);
}