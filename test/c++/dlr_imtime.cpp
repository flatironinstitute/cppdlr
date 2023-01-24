#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_basis.hpp>
#include <cppdlr/utils.hpp>
#include <cppdlr/dlr_kernels.hpp>

using namespace cppdlr;

static constexpr auto _ = nda::range::all;

nda::matrix<double> gfun(int norb, double beta, double t) {

  int npeak  = 5;

  auto g = nda::matrix<double>(norb,norb);
  g = 0;
  double rnd = 0;
  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      for (int l = 0; l < npeak; ++l) {
        rnd = sin(1000.0 * (i+2*j+3*l)); // Quick and dirty random number in [-1,1]
        g(i,j) += kfun(t, beta * rnd);
      }
    }
  }

  return g;
}

// Test DLR interpolation and evaluation for matrix-valued Green's function

TEST(dlr_imtime, interp_matrix) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  int norb = 2; // Orbital dimensions

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = dlr_rf.size();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();

  auto g = nda::array<double, 3>(norb, norb, r);

  for (int i = 0; i < r; ++i) { g(_, _, i) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes

  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error

  auto gtru  = nda::array<double, 3>(norb, norb, ntst);
  auto gtst  = nda::array<double, 3>(norb, norb, ntst);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru(_, _, i) = gfun(norb, beta, ttst(i));

    // [Q] Best way to clean this up?
    gtst(_, _, i) = itops.coefs2eval(gc, ttst(i));

    err = std::max(err, max_element(abs(gtru(_, _, i) - gtst(_, _, i))));
  }

  EXPECT_LT(err, 10 * eps);
}


// Test DLR interpolation and evaluation for scalar-valued Green's function

TEST(dlr_imtime, interp_scalar) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = dlr_rf.size();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();

  auto g = nda::vector<double>(r);

  for (int i = 0; i < r; ++i) { g(i) = gfun(1, beta, dlr_it(i))(0,0); }

  // DLR coefficients of G
  auto gc2 = itops.vals2coefs(g);
  auto gc = nda::array<double, 3>(1,1,r);
  gc(0,0,_) = gc2;

  // Check that G can be recovered at imaginary time nodes

  EXPECT_LT(max_element(abs(itops.coefs2vals(gc2) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error

  auto gtru  = nda::vector<double>(ntst);
  auto gtst  = nda::vector<double>(ntst);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru(i) = gfun(1, beta, ttst(i))(0,0);

    // [Q] Best way to clean this up?
    gtst(i) = itops.coefs2eval(gc, ttst(i))(0,0);

    err = std::max(err, abs(gtru(i) - gtst(i)));
  }

  EXPECT_LT(err, 10 * eps);
}