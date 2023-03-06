#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_basis.hpp>
#include <cppdlr/utils.hpp>
#include <cppdlr/dlr_kernels.hpp>
#include <nda/blas.hpp>

using namespace cppdlr;
using namespace nda;

static constexpr auto _ = nda::range::all;

// Green's function G_ij(tau) = sum_l c_ijl K(tau,om_ijl) with random c_ijl, om_ijl

nda::matrix<double> gfun(int norb, double beta, double t) {

  int npeak = 5;

  auto g    = nda::matrix<double>(norb, norb);
  g         = 0;
  auto c    = nda::vector<double>(npeak);
  double om = 0;
  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      // Get random weights that sum to 1
      for (int l = 0; l < npeak; ++l) {
        c(l) = (sin(1000.0 * (i + 2 * j + 3 * l + 7)) + 1) / 2; // Quick and dirty rand # gen on [0,1]
      }
      c = c / sum(c);

      // Evaluate Green's function
      for (int l = 0; l < npeak; ++l) {
        om = sin(2000.0 * (3 * i + 2 * j + l + 6)); // Rand # on [-1,1]
        g(i, j) += c(l) * kfun(t, beta * om);
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
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();

  auto g = nda::array<double, 3>(r, norb, norb);

  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes

  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error

  auto gtru  = nda::array<double, 3>(ntst, norb, norb);
  auto gtst  = nda::array<double, 3>(ntst, norb, norb);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru(i, _, _) = gfun(norb, beta, ttst(i));
    gtst(i, _, _) = itops.coefs2eval(gc, ttst(i));

    err = std::max(err, max_element(abs(gtru(i, _, _) - gtst(i, _, _))));
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
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();

  auto g = nda::vector<double>(r);

  for (int i = 0; i < r; ++i) { g(i) = gfun(1, beta, dlr_it(i))(0, 0); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes

  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error

  auto gtru  = nda::vector<double>(ntst);
  auto gtst  = nda::vector<double>(ntst);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru(i) = gfun(1, beta, ttst(i))(0, 0);

    gtst(i) = itops.coefs2eval(gc, ttst(i));

    err = std::max(err, abs(gtru(i) - gtst(i)));
  }

  EXPECT_LT(err, 10 * eps);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method

  auto kvec = itops.get_kevalvec(ttst(0));
  EXPECT_LT((abs(blas::dot(gc,kvec) - gtst(0))),1e-14);

  PRINT((abs(blas::dot(gc,kvec) - gtst(0))));

}
