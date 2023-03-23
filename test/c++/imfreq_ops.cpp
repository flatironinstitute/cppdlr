/** 
* @file imfreq_ops.cpp
*
* @brief Tests for imfreq_ops class.
*/

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>

#include <cppdlr/dlr_imfreq.hpp>
#include <cppdlr/dlr_basis.hpp>
#include <cppdlr/utils.hpp>
#include <cppdlr/dlr_kernels.hpp>
#include <nda/blas.hpp>

using namespace cppdlr;
using namespace nda;

static constexpr auto _ = nda::range::all;

/**
* @brief Green's function which is a random sum of poles
*
* G_ij(iom_n) = sum_l c_ijl K(n,om_ijl) with random c_ijl, om_ijl
*
* @param[in] norb Number of orbital indices
* @param[in] beta Inverse temperature
* @param[in] n    Imaginary frequency evaluation point index
* @param[in] xi   Fermionic (xi = -1) or bosonic (xi = 1) Matsubara frequencies
*
* @return Green's function evaluated at iom_n
*/
nda::matrix<dcomplex> gfun(int norb, double beta, int n, int xi) {

  int npeak = 5;

  auto g    = nda::matrix<dcomplex>(norb, norb);
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
        g(i, j) += c(l) * kfun_if(2 * n + (1 - xi) / 2, beta * om);
      }
    }
  }

  return g;
}

/**
* @brief Test DLR interpolation and evaluation for matrix-valued Green's
* function
*/
TEST(imfreq_ops, interp_matrix) {

  double lambda = 1000;   // DLR cutoff
  double eps    = 1e-10;  // DLR tolerance
  int xi        = -1;     // Fermionic Green's function

  double beta = 1000;     // Inverse temperature
  int nmaxtst = 5000;     // # imag time test points

  int norb = 2; // Orbital dimensions

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, xi);

  // Sample Green's function G at DLR imaginary frequency nodes
  int r = ifops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_if = ifops.get_ifnodes();
  auto g = nda::array<dcomplex, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), xi); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(g);

  // Check that G can be recovered at imaginary frequency nodes
  EXPECT_LT(max_element(abs(ifops.coefs2vals(gc) - g)), 1e-14);

  // Compute L infinity error
  auto gtru  = nda::matrix<dcomplex>(norb, norb);
  auto gtst  = nda::matrix<dcomplex>(norb, norb);
  double err = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru = gfun(norb, beta, n, xi);
    gtst = ifops.coefs2eval(gc, n);
    err = std::max(err, max_element(abs(gtru - gtst)));
  }

  EXPECT_LT(err, 10 * eps);
}

/**
* @brief Test DLR interpolation and evaluation for scalar-valued Green's
* function
*/
TEST(imfreq_ops, interp_scalar) {

  double lambda = 1000;   // DLR cutoff
  double eps    = 1e-10;  // DLR tolerance
  int xi        = -1;     // Fermionic Green's function

  double beta = 1000;     // Inverse temperature
  int nmaxtst = 5000;     // # imag time test points

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, xi);

  // Sample Green's function G at DLR imaginary frequency nodes
  int r = ifops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_if = ifops.get_ifnodes();
  auto g = nda::vector<dcomplex>(r);
  for (int i = 0; i < r; ++i) { g(i) = gfun(1, beta, dlr_if(i), xi)(0, 0); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(g);

  // Check that G can be recovered at imaginary frequency nodes
  EXPECT_LT(max_element(abs(ifops.coefs2vals(gc) - g)), 1e-14);

  // Compute L infinity error
  std::complex<double> gtru  = 0, gtst = 0;
  double err = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru = gfun(1, beta, n, xi)(0, 0);
    gtst = ifops.coefs2eval(gc, n);
    err = std::max(err, abs(gtru - gtst));
  }

  EXPECT_LT(err, 10 * eps);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method
  gtst = ifops.coefs2eval(gc, 3);
  auto kvec = ifops.get_kevalvec(3);
  auto zgc = nda::vector<dcomplex>(gc);
  EXPECT_LT((abs(blas::dotc(zgc,kvec) - gtst)),1e-14);

}

TEST(dlr_imfreq, h5_rw) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance
  int xi        = -1;    // Fermionic Green's function

  // Get DLR frequencies
  auto dlr_rf = dlr_freq(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, xi);

  auto filename = "data_imfreq_ops_h5_rw.h5";
  auto name = "ifops";

  {
    h5::file file(filename, 'w');
    h5_write(file, name, ifops);
  }

  imfreq_ops ifops_ref;
  {
    h5::file file(filename, 'r');
    h5_read(file, name, ifops_ref);
  }

  // Check equal
  EXPECT_EQ(ifops.lambda(), ifops_ref.lambda());
  EXPECT_EQ(ifops.rank(), ifops_ref.rank());
  EXPECT_EQ_ARRAY(ifops.get_rfnodes(), ifops_ref.get_rfnodes());
  EXPECT_EQ_ARRAY(ifops.get_ifnodes(), ifops_ref.get_ifnodes());
  EXPECT_EQ_ARRAY(ifops.get_cf2if(), ifops_ref.get_cf2if());
  EXPECT_EQ_ARRAY(ifops.get_if2cf_lu(), ifops_ref.get_if2cf_lu());
  EXPECT_EQ_ARRAY(ifops.get_if2cf_piv(), ifops_ref.get_if2cf_piv());
  
}
