/** 
* @file imtime_ops.cpp
*
* @brief Tests for imtime_ops class.
*/

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>

#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_basis.hpp>
#include <cppdlr/utils.hpp>
#include <cppdlr/dlr_kernels.hpp>
#include <nda/blas.hpp>

using namespace cppdlr;
using namespace nda;

static constexpr auto _ = nda::range::all;

/**
* @brief Green's function which is a random sum of exponentials
*
* G_ij(t) = sum_l c_ijl K(t,om_ijl) with random c_ijl, om_ijl
*
* @param norb Number of orbital indices
* @param beta Inverse temperature
* @param t    Imaginary time evaluation point
*
* @return Green's function evaluated at t
*/
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

/**
* @brief Test DLR interpolation and evaluation for matrix-valued Green's
* function
*/
TEST(imtime_ops, interp_matrix) {

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
  auto gtru  = nda::matrix<double>(norb, norb);
  auto gtst  = nda::matrix<double>(norb, norb);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru = gfun(norb, beta, ttst(i));
    gtst = itops.coefs2eval(gc, ttst(i));
    err = std::max(err, max_element(abs(gtru - gtst)));
  }

  EXPECT_LT(err, 10 * eps);

}

/**
* @brief Test DLR interpolation and evaluation for complex and matrix-valued Green's
* function
*/
TEST(imtime_ops, interp_matrix_complex) {

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
  auto g = nda::array<std::complex<double>, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute L infinity error
  auto gtru  = nda::matrix<std::complex<double>>(norb, norb);
  auto gtst  = nda::matrix<std::complex<double>>(norb, norb);
  double err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru = gfun(norb, beta, ttst(i));
    gtst = itops.coefs2eval(gc, ttst(i));
    err = std::max(err, max_element(abs(gtru - gtst)));
  }

  EXPECT_LT(err, 10 * eps);

}


/**
* @brief Test DLR interpolation and evaluation for scalar-valued Green's
* function
*/
TEST(imtime_ops, interp_scalar) {

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
  double gtru  = 0, gtst = 0, err = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru = gfun(1, beta, ttst(i))(0, 0);
    gtst = itops.coefs2eval(gc, ttst(i));
    err = std::max(err, abs(gtru - gtst));
  }

  EXPECT_LT(err, 10 * eps);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method
  auto kvec = itops.get_kevalvec(ttst(ntst-1));
  EXPECT_LT((abs(blas::dot(gc,kvec) - gtst)),1e-14);

}

TEST(dlr_imtime, h5_rw) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  auto filename = "data_imtime_ops_h5_rw.h5";
  auto name = "itops";

  {
    h5::file file(filename, 'w');
    h5_write(file, name, itops);
  }

  imtime_ops itops_ref;
  {
    h5::file file(filename, 'r');
    h5_read(file, name, itops_ref);
  }

  // Check equal
  EXPECT_EQ(itops.lambda(), itops_ref.lambda());
  EXPECT_EQ(itops.rank(), itops_ref.rank());
  EXPECT_EQ_ARRAY(itops.get_rfnodes(), itops_ref.get_rfnodes());
  EXPECT_EQ_ARRAY(itops.get_itnodes(), itops_ref.get_itnodes());
  EXPECT_EQ_ARRAY(itops.get_cf2it(), itops_ref.get_cf2it());
  EXPECT_EQ_ARRAY(itops.get_it2cf_lu(), itops_ref.get_it2cf_lu());
  EXPECT_EQ_ARRAY(itops.get_it2cf_piv(), itops_ref.get_it2cf_piv());
  
}
