// Copyright (c) 2023 Simons Foundation
// Copyright (c) 2023 Hugo Strand
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0.txt
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors: Hugo U. R. Strand, Nils Wentzell, Jason Kaye

/** 
* @file imfreq_ops.cpp
*
* @brief Tests for imfreq_ops class.
*/

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <cppdlr/cppdlr.hpp>
#include <nda/gtest_tools.hpp>
#include <fmt/format.h>
#include <random>

using namespace cppdlr;
using namespace nda;

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
        g(i, j) += c(l) * k_it(t, om, beta);
      }
    }
  }

  return g;
}

/**
* @brief Green's function which is a random sum of poles
*
* G_ij(iom_n) = sum_l c_ijl K(n,om_ijl) with random c_ijl, om_ijl
*
* @param[in] norb      Number of orbital indices
* @param[in] beta      Inverse temperature
* @param[in] n         Imaginary frequency evaluation point index
* @param[in] statistic Particle Statistic: Fermion or Boson
*
* @return Green's function evaluated at iom_n
*/
nda::matrix<dcomplex> gfun(int norb, double beta, int n, statistic_t statistic) {

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
      c = beta * c;

      // Evaluate Green's function
      for (int l = 0; l < npeak; ++l) {
        om = sin(2000.0 * (3 * i + 2 * j + l + 6)); // Rand # on [-1,1]
        g(i, j) += c(l) * k_if(n, beta * om, statistic);
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

  double lambda  = 1000;    // DLR cutoff
  double eps     = 1e-10;   // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic);

  // Sample Green's function G at DLR imaginary frequency nodes
  int r              = ifops.rank();
  auto const &dlr_if = ifops.get_ifnodes();
  auto g             = nda::array<dcomplex, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), statistic); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(beta, g);

  // Check that G can be recovered at imaginary frequency nodes
  EXPECT_LT(max_element(abs(ifops.coefs2vals(beta, gc) - g)), 5e-14);

  // Compute error in imaginary frequency
  auto gtru      = nda::matrix<dcomplex>(norb, norb);
  auto gtst      = nda::matrix<dcomplex>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru    = gfun(norb, beta, n, statistic);
    gtst    = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary time
  auto itops = imtime_ops(lambda, dlr_rf);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  auto gtru_it = nda::matrix<double>(norb, norb);
  auto gtst_it = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru_it = gfun(norb, beta, ttst(i));
    gtst_it = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru_it - gtst_it)));
    errl2 += pow(frobenius_norm(gtru_it - gtst_it), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 200 * eps);
  EXPECT_LT(errl2, 3 * eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test DLR interpolation and evaluation for scalar-valued Green's
* function
*/
TEST(imfreq_ops, interp_scalar) {

  double lambda  = 1000;    // DLR cutoff
  double eps     = 1e-10;   // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  double beta = 1000; // Inverse temperature
  int nmaxtst = 5000; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic);

  // Sample Green's function G at DLR imaginary frequency nodes
  int r              = ifops.rank();
  auto const &dlr_if = ifops.get_ifnodes();
  auto g             = nda::vector<dcomplex>(r);
  for (int i = 0; i < r; ++i) { g(i) = gfun(1, beta, dlr_if(i), statistic)(0, 0); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(beta, g);

  // Check that G can be recovered at imaginary frequency nodes
  EXPECT_LT(max_element(abs(ifops.coefs2vals(beta, gc) - g)), 1e-13);

  // Compute error in imaginary frequency
  std::complex<double> gtru = 0, gtst = 0;
  double errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru    = gfun(1, beta, n, statistic)(0, 0);
    gtst    = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, abs(gtru - gtst));
    errl2 += pow(abs(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method
  gtst      = ifops.coefs2eval(beta, gc, 3);
  auto kvec = ifops.build_evalvec(beta, 3);
  auto zgc  = nda::vector<dcomplex>(gc);
  EXPECT_LT((abs(blas::dotc(zgc, kvec) - gtst)), 1e-13);
}

/**
* @brief Test symmetrized DLR interpolation and evaluation for matrix-valued
* Green's function
*/
TEST(imfreq_ops, interp_matrix_sym_fer) {

  double lambda  = 1000;    // DLR cutoff
  double eps     = 1e-10;   // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps, SYM);
  int r       = dlr_rf.size();

  // Verify symmetry
  EXPECT_EQ(max_element(abs(dlr_rf(range(r / 2)) + dlr_rf(range(r - 1, r / 2 - 1, -1)))), 0);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic, SYM);

  // Sample Green's function G at DLR imaginary frequency nodes
  auto const &dlr_if = ifops.get_ifnodes();

  // Verify symmetry
  EXPECT_EQ(max_element(abs(2 * dlr_if(range(r / 2)) + 1 + 2 * dlr_if(range(r - 1, r / 2 - 1, -1)) + 1)), 0);

  auto g = nda::array<dcomplex, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), statistic); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(beta, g);

  // Check that G can be recovered at imaginary frequency nodes
  EXPECT_LT(max_element(abs(ifops.coefs2vals(beta, gc) - g)), 3e-13);

  // Compute error in imaginary frequency
  auto gtru      = nda::matrix<dcomplex>(norb, norb);
  auto gtst      = nda::matrix<dcomplex>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru    = gfun(norb, beta, n, statistic);
    gtst    = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary time
  auto itops = imtime_ops(lambda, dlr_rf);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  auto gtru_it = nda::matrix<double>(norb, norb);
  auto gtst_it = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru_it = gfun(norb, beta, ttst(i));
    gtst_it = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru_it - gtst_it)));
    errl2 += pow(frobenius_norm(gtru_it - gtst_it), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test symmetrized DLR interpolation and evaluation for bosonic
* matrix-valued Green's function
*/
TEST(imfreq_ops, interp_matrix_sym_bos) {

  double lambda  = 1000;  // DLR cutoff
  double eps     = 1e-10; // DLR tolerance
  auto statistic = Boson; // Bosonic Green's function

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps, SYM);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic, SYM);

  // Obtain DLR imaginary frequency nodes
  auto const &dlr_if = ifops.get_ifnodes();
  int niom           = dlr_if.size();

  // Verify symmetry
  EXPECT_EQ(max_element(abs(2 * dlr_if(range((niom - 1) / 2)) + 2 * dlr_if(range(niom - 1, (niom - 1) / 2, -1)))), 0);

  // Verify n = 0 was selected
  EXPECT_EQ(dlr_if((niom - 1) / 2), 0);

  // Sample Green's function at DLR nodes
  auto g = nda::array<dcomplex, 3>(niom, norb, norb);
  for (int i = 0; i < niom; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), statistic); }

  // DLR coefficients of G
  auto gc = ifops.vals2coefs(beta, g);

  // Compute error in imaginary frequency
  auto gtru      = nda::matrix<dcomplex>(norb, norb);
  auto gtst      = nda::matrix<dcomplex>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru    = gfun(norb, beta, n, statistic);
    gtst    = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary time
  auto itops = imtime_ops(lambda, dlr_rf);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  auto gtru_it = nda::matrix<double>(norb, norb);
  auto gtst_it = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru_it = gfun(norb, beta, ttst(i));
    gtst_it = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru_it - gtst_it)));
    errl2 += pow(frobenius_norm(gtru_it - gtst_it), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

TEST(dlr_imfreq, h5_rw) {

  double lambda  = 1000;    // DLR cutoff
  double eps     = 1e-10;   // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic);

  auto filename = "data_imfreq_ops_h5_rw.h5";
  auto name     = "ifops";

  {
    h5::file file(filename, 'w');
    h5::write(file, name, ifops);
  }

  imfreq_ops ifops_ref;
  {
    h5::file file(filename, 'r');
    h5::read(file, name, ifops_ref);
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

/**
* @brief Test direct application of if2it matrix against standard vals2coefs -> coefs2vals
* pipeline, for all statistic/symmetrization combinations
*/
TEST(imfreq_ops, if2it) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000; // Inverse temperature
  int norb    = 2;    // Orbital dimensions

  // Test cases: {statistic, symmetrize}
  std::vector<std::pair<statistic_t, bool>> cases = {{Fermion, NONSYM}, {Fermion, SYM}, {Boson, NONSYM}, {Boson, SYM}};

  for (auto const &[statistic, symmetrize] : cases) {

    std::string label = fmt::format("{}, {}", statistic == Fermion ? "Fermion" : "Boson", symmetrize ? "SYM" : "NONSYM");

    auto dlr_rf = build_dlr_rf(lambda, eps, symmetrize);

    auto ifops = imfreq_ops(lambda, dlr_rf, statistic, symmetrize);
    auto itops = imtime_ops(lambda, dlr_rf);

    int r              = ifops.rank();
    int niom           = ifops.get_ifnodes().size();
    auto const &dlr_if = ifops.get_ifnodes();

    std::cout << fmt::format("if2it ({}): r = {}, niom = {}\n", label, r, niom);

    // Sample Green's function at DLR imaginary frequency nodes
    auto g = nda::array<dcomplex, 3>(niom, norb, norb);
    for (int i = 0; i < niom; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), statistic); }

    // Standard approach: vals2coefs then coefs2vals
    auto gc       = ifops.vals2coefs(beta, g);
    auto g_it_ref = itops.coefs2vals(gc);

    // Using if2it matrix directly
    auto if2it                        = build_if2it(ifops, itops);
    auto g_it_new                     = nda::array<dcomplex, 3>(r, norb, norb);
    reshape(g_it_new, r, norb * norb) = matmul(if2it, nda::matrix<dcomplex>(nda::reshape(g, niom, norb * norb))) / beta;

    // Compare
    auto err = max_element(abs(g_it_ref - g_it_new));
    std::cout << fmt::format("build_if2it ({}) error: {:e}\n", label, err);
    EXPECT_LT(err, 5e-13);
  }
}

// FAILING reproducer of the symmetrized-DLR conditioning issue that breaks the
// tprf DLR Eliashberg solver. The symmetric grid places a near-degenerate
// +/-delta pole pair around omega=0, which ill-conditions the DLR fit. A
// DLR-built operator is then no longer linear on non-spectral inputs -- exactly
// the vectors the tprf eigensolver (ARPACK) feeds it. EXPECTED TO FAIL here;
// passes with the self-symmetric-fixed-point fix. See issue_summary.md.
TEST(imfreq_ops, symmetrized_operator_linearity) {

  double beta = 2.0, lambda = 10.0, eps = 1e-8; // tprf DLR Eliashberg parameters

  // Relative additivity residual |M(a+b) - M(a) - M(b)| of one tprf Eliashberg
  // matvec M (eliashberg_product_fft, single k-point + scalar target). M is
  // exactly linear in delta, so a faithful operator gives ~0.
  auto residual = [&](bool symmetrize) {
    auto rf    = build_dlr_rf(lambda, eps, symmetrize);
    int r      = rf.size();
    auto itops = imtime_ops(lambda, rf, symmetrize);
    auto ifops = imfreq_ops(lambda, rf, Fermion, symmetrize);
    int niom   = ifops.get_ifnodes().size();

    // Fixed (delta-independent) operator data: the G*G bubble as DLR coefficients
    // and the dynamic interaction in imaginary time. These keep M linear in delta.
    auto gg = nda::vector<dcomplex>(r), gam = nda::vector<dcomplex>(r);
    for (int i = 0; i < r; ++i) {
      gg(i)  = 0.5 + 0.25 * std::cos(3.0 * i);
      gam(i) = 0.5 + 0.25 * std::cos(2.0 * i);
    }

    auto M = [&](nda::vector<dcomplex> d) {                          // the cppdlr calls tprf makes:
      auto ft = itops.convolve(beta, gg, ifops.vals2coefs(beta, d)); // F = (G*G) (*) delta, in tau
      auto F  = ifops.coefs2vals(beta, itops.vals2coefs(ft));        //   -> back to imfreq
      auto Ft = itops.coefs2vals(ifops.vals2coefs(beta, F));         // fourier_wr_to_tr: imfreq -> imtime
      for (int i = 0; i < Ft.size(); ++i) Ft(i) *= gam(i);           // multiply by Gamma(tau)
      return ifops.coefs2vals(beta, itops.vals2coefs(Ft));           // fourier_tr_to_wr: imtime -> imfreq
    };

    std::mt19937 gen(12345); // fixed seed: deterministic
    std::normal_distribution<double> dist(0.0, 1.0);
    auto a = nda::vector<dcomplex>(niom), b = nda::vector<dcomplex>(niom);
    for (int i = 0; i < niom; ++i) {
      a(i) = dcomplex(dist(gen), dist(gen)); // non-spectral (random) inputs, like the ARPACK iterates
      b(i) = dcomplex(dist(gen), dist(gen));
    }
    auto Ma = M(a), Mb = M(b), Mab = M(make_regular(a + b));
    return max_element(abs(Mab - Ma - Mb)) / (max_element(abs(Ma)) + max_element(abs(Mb)));
  };

  double res_nonsym = residual(NONSYM), res_sym = residual(SYM);
  std::cout << fmt::format("Eliashberg operator linearity residual: nonsym = {:.1e}, sym = {:.1e}\n", res_nonsym, res_sym);

  // On the symmetric grid the operator is ~10000x less linear here (sym ~7e-5 vs
  // nonsym ~6e-9): a manifestly non-linear operator, which derails ARPACK.
  EXPECT_LT(res_sym, 5 * res_nonsym);
}
