// Copyright (c) 2022-2023 Simons Foundation
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
* @file imtime_ops.cpp
*
* @brief Tests for imtime_ops class.
*/

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <cppdlr/cppdlr.hpp>
#include <nda/gtest_tools.hpp>
#include <fmt/format.h>

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
TEST(imtime_ops, interp_matrix) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();
  auto g             = nda::array<double, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<double>(norb, norb);
  auto gtst      = nda::matrix<double>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary frequency
  auto ifops = imfreq_ops(lambda, dlr_rf, Fermion);

  auto gtru_if = nda::matrix<dcomplex>(norb, norb);
  auto gtst_if = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru_if = gfun(norb, beta, n, Fermion);
    gtst_if = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru_if - gtst_if)));
    errl2 += pow(frobenius_norm(gtru_if - gtst_if), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
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

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();
  auto g             = nda::array<dcomplex, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<dcomplex>(norb, norb);
  auto gtst      = nda::matrix<dcomplex>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, 2 * eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
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

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();
  auto g             = nda::vector<double>(r);
  for (int i = 0; i < r; ++i) { g(i) = gfun(1, beta, dlr_it(i))(0, 0); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  double gtru = 0, gtst = 0, errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(1, beta, ttst(i))(0, 0);
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, abs(gtru - gtst));
    errl2 += pow(gtru - gtst, 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method
  auto kvec = itops.build_evalvec(ttst(ntst - 1));
  EXPECT_LT((abs(blas::dot(gc, kvec) - gtst)), 1e-14);
}

/**
* @brief Test DLR fitting for real matrix-valued Green's function
*/
TEST(imtime_ops, fit_matrix) {

  double lambda = 1000; // DLR cutoff
  double eps    = 1e-8; // DLR tolerance

  double beta  = 1000;  // Inverse temperature
  int nsample  = 5000;  // # imag time sampling points
  double noise = 1e-6;  // Noise level
  int ntst     = 10000; // # imag time test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}, noise = {:e}\n", eps, lambda, noise);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at equispaced imaginary time nodes and add
  // uniform random noise
  auto t = eqptsrel(nsample);
  auto g = nda::array<double, 3>(nsample, norb, norb);
  for (int i = 0; i < nsample; ++i) { g(i, _, _) = gfun(norb, beta, t(i)) + noise * (2 * nda::rand() - 1); }

  // DLR coefficients of G
  auto gc = itops.fitvals2coefs(t, g);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<double>(norb, norb);
  auto gtst      = nda::matrix<double>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * noise);
  EXPECT_LT(errl2, 2 * noise);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test DLR fitting for complex matrix-valued Green's function
*/
TEST(imtime_ops, fit_matrix_cmplx) {

  double lambda = 1000; // DLR cutoff
  double eps    = 1e-8; // DLR tolerance

  double beta  = 1000;  // Inverse temperature
  int nsample  = 5000;  // # imag time sampling points
  double noise = 1e-6;  // Noise level
  int ntst     = 10000; // # imag time test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}, noise = {:e}\n", eps, lambda, noise);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at equispaced imaginary time nodes and add
  // uniform random noise
  auto t = eqptsrel(nsample);
  auto g = nda::array<dcomplex, 3>(nsample, norb, norb);
  for (int i = 0; i < nsample; ++i) { g(i, _, _) = gfun(norb, beta, t(i)) + noise * (2 * (nda::rand() + 1i * nda::rand()) - 1); }

  // DLR coefficients of G
  auto gc = itops.fitvals2coefs(t, g);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<dcomplex>(norb, norb);
  auto gtst      = nda::matrix<dcomplex>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * noise);
  EXPECT_LT(errl2, 2 * noise);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test DLR fitting for scalar-valued Green's function
*/
TEST(imtime_ops, fit_scalar) {

  double lambda = 1000; // DLR cutoff
  double eps    = 1e-8; // DLR tolerance

  double beta  = 1000;  // Inverse temperature
  int nsample  = 5000;  // # imag time sampling points
  double noise = 1e-6;  // Noise level
  int ntst     = 10000; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}, noise = {:e}\n", eps, lambda, noise);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at equispaced imaginary time nodes and add
  // uniform random noise
  auto t = eqptsrel(nsample);
  auto g = nda::vector<double>(nsample);
  for (int i = 0; i < nsample; ++i) { g(i) = gfun(1, beta, t(i))(0, 0) + noise * (2 * nda::rand() - 1); }

  // DLR coefficients of G
  auto gc = itops.fitvals2coefs(t, g);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  double gtru = 0, gtst = 0, errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(1, beta, ttst(i))(0, 0);
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, abs(gtru - gtst));
    errl2 += pow(gtru - gtst, 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * noise);
  EXPECT_LT(errl2, 2 * noise);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Test that constructing vector of evaluation at a point and then applying to
  // coefficients gives same result as direct evaluation method
  auto kvec = itops.build_evalvec(ttst(ntst - 1));
  EXPECT_LT((abs(blas::dot(gc, kvec) - gtst)), 1e-14);
}

/**
* @brief Test convolution and time-ordered convolution of two real-valued
* Green's functions
*
* We use Green's functions f and g given by a single exponential, so that the
* result of the convolution is easy to compute analytically.
*/
TEST(imtime_ops, convolve_scalar_real) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-12; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Specify frequency of single exponentials to be used for f and g
  double omf = 0.1234, omg = -0.5678;

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it = itops.get_itnodes();
  auto f             = nda::array<double, 1>(r);
  auto g             = nda::array<double, 1>(r);
  for (int i = 0; i < r; ++i) { f(i) = k_it(dlr_it(i), omf, beta); };
  for (int i = 0; i < r; ++i) { g(i) = k_it(dlr_it(i), omg, beta); };

  // Get DLR coefficients of f and g
  auto fc = itops.vals2coefs(f);
  auto gc = itops.vals2coefs(g);

  // Get convolution and time-ordered convolution of f and g directly
  auto h  = itops.convolve(beta, Fermion, fc, gc);
  auto ht = itops.convolve(beta, Fermion, fc, gc, TIME_ORDERED);

  // Get convolution and time-ordered convolution of f and g by first forming
  // matrix of convolution by f and then applying it to g
  auto h2  = itops.convolve(itops.convmat(beta, Fermion, fc), g);
  auto ht2 = itops.convolve(itops.convmat(beta, Fermion, fc, TIME_ORDERED), g);

  // Check that the two methods give the same result
  EXPECT_LT(max_element(abs(h - h2)), 1e-14);
  EXPECT_LT(max_element(abs(ht - ht2)), 1e-14);

  // Check error of convolution and time-ordered convolution

  auto hc  = itops.vals2coefs(h);  // DLR coefficients of h
  auto htc = itops.vals2coefs(ht); // DLR coefficients of ht

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  double gtru = 0, gtst = 0, gttru = 0, gttst = 0;
  double errlinf = 0, errl2 = 0, errtlinf = 0, errtl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru     = (k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta)) / (omg - omf);                                               // Exact result
    gttru    = (k_it(0.0, omf, beta) * k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta) * k_it(0.0, omg, beta)) / (omf - omg); // Exact result
    gtst     = itops.coefs2eval(hc, ttst(i));
    gttst    = itops.coefs2eval(htc, ttst(i));
    errlinf  = std::max(errlinf, abs(gtru - gtst));
    errtlinf = std::max(errtlinf, abs(gttru - gttst));
    errl2 += pow(gtru - gtst, 2);
    errtl2 += pow(gttru - gttst, 2);
  }
  errl2  = sqrt(errl2 / ntst);
  errtl2 = sqrt(errtl2 / ntst);

  EXPECT_LT(errlinf, 25 * eps);
  EXPECT_LT(errtlinf, 25 * eps);
  EXPECT_LT(errl2, 4 * eps);
  EXPECT_LT(errtl2, 3 * eps);
  std::cout << fmt::format("Ordinary convolution: L^inf err = {:e}, L^2 err = {:e}\n", errlinf, errl2);
  std::cout << fmt::format("Time-ordered convolution: L^inf err = {:e}, L^2 err = {:e}\n", errtlinf, errtl2);
}

/**
* @brief Test convolution and time-ordered convolution of two complex-valued
* Green's functions
*
* We use Green's functions f and g given by a scalar multiple of a single exponential, so that the
* result of the convolution is easy to compute analytically.
*/
TEST(imtime_ops, convolve_scalar_cmplx) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-12; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Specify frequency of single exponentials to be used for f and g
  double omf = 0.1234, omg = -0.5678;

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r = itops.rank();
  // [Q] Is this correct or just auto?
  auto const &dlr_it      = itops.get_itnodes();
  auto f                  = nda::array<dcomplex, 1>(r);
  auto g                  = nda::array<dcomplex, 1>(r);
  std::complex<double> c1 = (1.0 + 2.0i) / 3.0, c2 = (2.0 + 1.0i) / 3.0;
  for (int i = 0; i < r; ++i) { f(i) = c1 * k_it(dlr_it(i), omf, beta); };
  for (int i = 0; i < r; ++i) { g(i) = c2 * k_it(dlr_it(i), omg, beta); };

  // Get DLR coefficients of f and g
  auto fc = itops.vals2coefs(f);
  auto gc = itops.vals2coefs(g);

  // Get convolution and time-ordered convolution of f and g directly
  auto h  = itops.convolve(beta, Fermion, fc, gc);
  auto ht = itops.convolve(beta, Fermion, fc, gc, TIME_ORDERED);

  // Get convolution and time-ordered convolution of f and g by first forming
  // matrix of convolution by f and then applying it to g
  auto h2  = itops.convolve(itops.convmat(beta, Fermion, fc), g);
  auto ht2 = itops.convolve(itops.convmat(beta, Fermion, fc, TIME_ORDERED), g);

  // Check that the two methods give the same result
  EXPECT_LT(max_element(abs(h - h2)), 1e-14);
  EXPECT_LT(max_element(abs(ht - ht2)), 1e-14);

  // Check error of convolution and time-ordered convolution

  auto hc  = itops.vals2coefs(h);  // DLR coefficients of h
  auto htc = itops.vals2coefs(ht); // DLR coefficients of ht

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  std::complex<double> gtru = 0, gttru = 0, gtst = 0, gttst = 0;
  double errlinf = 0, errl2 = 0, errtlinf = 0, errtl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru = c1 * c2 * (k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta)) / (omg - omf); // Exact result
    gttru =
       c1 * c2 * (k_it(0.0, omf, beta) * k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta) * k_it(0.0, omg, beta)) / (omf - omg); // Exact result
    gtst     = itops.coefs2eval(hc, ttst(i));
    gttst    = itops.coefs2eval(htc, ttst(i));
    errlinf  = std::max(errlinf, abs(gtru - gtst));
    errtlinf = std::max(errtlinf, abs(gttru - gttst));
    errl2 += pow(abs(gtru - gtst), 2);
    errtl2 += pow(abs(gttru - gttst), 2);
  }
  errl2  = sqrt(errl2 / ntst);
  errtl2 = sqrt(errtl2 / ntst);

  EXPECT_LT(errlinf, 13 * eps);
  EXPECT_LT(errtlinf, 13 * eps);
  EXPECT_LT(errl2, 2 * eps);
  EXPECT_LT(errtl2, 2 * eps);
  std::cout << fmt::format("Ordinary convolution: L^inf err = {:e}, L^2 err = {:e}\n", errlinf, errl2);
  std::cout << fmt::format("Time-ordered convolution: L^inf err = {:e}, L^2 err = {:e}\n", errtlinf, errtl2);
}

/**
* @brief Test convolution and time-ordered convolution of two real matrix-valued
* Green's functions
*
* We use Green's functions f and g given by products of single exponentials and a
* matrix of ones, so that the result of the convolution is easy to compute
* analytically
*/
TEST(imtime_ops, convolve_matrix_real) {

  double lambda = 10;    // DLR cutoff
  double eps    = 1e-12; // DLR tolerance

  double beta = 10;  // Inverse temperature
  int ntst    = 100; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  int norb2 = 2; // Orbital dimensions
  for (int norb1 = 2; norb1 <= 3; ++norb1) {
    std::cout << fmt::format("norb1 = {}, norb2 = {}\n", norb1, norb2);

    // Specify frequency of single exponentials to be used for f and g
    double omf = 0.1234, omg = -0.5678;

    // Get DLR frequencies
    auto dlr_rf = build_dlr_rf(lambda, eps);

    // Get DLR imaginary time object
    auto itops = imtime_ops(lambda, dlr_rf);

    // Sample Green's function G at DLR imaginary time nodes
    int r = itops.rank();
    // [Q] Is this correct or just auto?
    auto const &dlr_it = itops.get_itnodes();
    auto f             = nda::array<double, 3>(r, norb1, norb2);
    auto g             = nda::array<double, 3>(r, norb2, norb2);
    for (int i = 0; i < r; ++i) { f(i, _, _) = k_it(dlr_it(i), omf, beta) / sqrt(1.0 * norb2); };
    for (int i = 0; i < r; ++i) { g(i, _, _) = k_it(dlr_it(i), omg, beta) / sqrt(1.0 * norb2); };

    // Get DLR coefficients of f and g
    auto fc = itops.vals2coefs(f);
    auto gc = itops.vals2coefs(g);

    // Get convolution and time-ordered convolution of f and g directly
    auto h  = itops.convolve(beta, Fermion, fc, gc);
    auto ht = itops.convolve(beta, Fermion, fc, gc, TIME_ORDERED);

    auto h2  = itops.convolve(itops.convmat(beta, Fermion, fc), g);
    auto ht2 = itops.convolve(itops.convmat(beta, Fermion, fc, TIME_ORDERED), g);

    // Check that the two methods give the same result
    EXPECT_LT(max_element(abs(h - h2)), 1e-14);
    EXPECT_LT(max_element(abs(ht - ht2)), 1e-14);

    // Check error of convolution and time-ordered convolution
    auto hc  = itops.vals2coefs(h);  // DLR coefficients of h
    auto htc = itops.vals2coefs(ht); // DLR coefficients of ht

    // Get test points in relative format
    auto ttst = eqptsrel(ntst);

    // Compute error in imaginary time
    auto htru      = nda::array<double, 2>(norb1, norb2);
    auto htst      = nda::array<double, 2>(norb1, norb2);
    auto httru     = nda::array<double, 2>(norb1, norb2);
    auto httst     = nda::array<double, 2>(norb1, norb2);
    double errlinf = 0, errl2 = 0, errtlinf = 0, errtl2 = 0;
    for (int i = 0; i < ntst; ++i) {
      htru     = (k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta)) / (omg - omf);                                               // Exact result
      httru    = (k_it(0.0, omf, beta) * k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta) * k_it(0.0, omg, beta)) / (omf - omg); // Exact result
      htst     = itops.coefs2eval(hc, ttst(i));
      httst    = itops.coefs2eval(htc, ttst(i));
      errlinf  = std::max(errlinf, max_element(abs(htru - htst)));
      errtlinf = std::max(errtlinf, max_element(abs(httru - httst)));
      errl2 += pow(frobenius_norm(htru - htst), 2);
      errtl2 += pow(frobenius_norm(httru - httst), 2);
    }
    errl2  = sqrt(errl2 / ntst);
    errtl2 = sqrt(errtl2 / ntst);

    EXPECT_LT(errlinf, 23 * eps);
    EXPECT_LT(errl2, 7 * eps);
    EXPECT_LT(errtlinf, 23 * eps);
    EXPECT_LT(errtl2, 6 * eps);
    std::cout << fmt::format("Ordinary convolution: L^inf err = {:e}, L^2 err = {:e}\n", errlinf, errl2);
    std::cout << fmt::format("Time-ordered convolution: L^inf err = {:e}, L^2 err = {:e}\n", errtlinf, errtl2);
  }
}

/**
* @brief Test convolution and time-ordered convolution of two complex
* matrix-valued Green's functions
*
* We use Green's functions f and g given by scalar multiples of products of
* single exponentials and a matrix of ones, so that the result of the convolution
* is easy to compute analytically
*/
TEST(imtime_ops, convolve_matrix_cmplx) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-12; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  int norb2 = 2; // Orbital dimensions
  for (int norb1 = 2; norb1 <= 3; ++norb1) {
    std::cout << fmt::format("norb1 = {}, norb2 = {}\n", norb1, norb2);

    // Specify frequency of single exponentials to be used for f and g
    double omf = 0.1234, omg = -0.5678;

    // Get DLR frequencies
    auto dlr_rf = build_dlr_rf(lambda, eps);

    // Get DLR imaginary time object
    auto itops = imtime_ops(lambda, dlr_rf);

    // Sample Green's function G at DLR imaginary time nodes
    int r = itops.rank();
    // [Q] Is this correct or just auto?
    auto const &dlr_it      = itops.get_itnodes();
    auto f                  = nda::array<dcomplex, 3>(r, norb1, norb2);
    auto g                  = nda::array<dcomplex, 3>(r, norb2, norb2);
    std::complex<double> c1 = (1.0 + 2.0i) / 3.0, c2 = (2.0 + 1.0i) / 3.0;
    for (int i = 0; i < r; ++i) { f(i, _, _) = c1 * k_it(dlr_it(i), omf, beta) / sqrt(1.0 * norb2); };
    for (int i = 0; i < r; ++i) { g(i, _, _) = c2 * k_it(dlr_it(i), omg, beta) / sqrt(1.0 * norb2); };

    // Get DLR coefficients of f and g
    auto fc = itops.vals2coefs(f);
    auto gc = itops.vals2coefs(g);

    // Get convolution and time-ordered convolution of f and g directly
    auto h  = itops.convolve(beta, Fermion, fc, gc);
    auto ht = itops.convolve(beta, Fermion, fc, gc, TIME_ORDERED);

    // Get convolution and time-ordered convolution of f and g by first forming
    // matrix of convolution by f and then applying it to g
    auto h2  = itops.convolve(itops.convmat(beta, Fermion, fc), g);
    auto ht2 = itops.convolve(itops.convmat(beta, Fermion, fc, TIME_ORDERED), g);

    // Check that the two methods give the same result
    EXPECT_LT(max_element(abs(h - h2)), 1e-14);
    EXPECT_LT(max_element(abs(ht - ht2)), 1e-14);

    // Check error of convolution and time-ordered convolution

    auto hc  = itops.vals2coefs(h);  // DLR coefficients of h
    auto htc = itops.vals2coefs(ht); // DLR coefficients of ht

    // Get test points in relative format
    auto ttst = eqptsrel(ntst);

    // Compute error in imaginary time
    auto htru      = nda::array<dcomplex, 2>(norb1, norb2);
    auto htst      = nda::array<dcomplex, 2>(norb1, norb2);
    auto httru     = nda::array<dcomplex, 2>(norb1, norb2);
    auto httst     = nda::array<dcomplex, 2>(norb1, norb2);
    double errlinf = 0, errl2 = 0, errtlinf = 0, errtl2 = 0;
    for (int i = 0; i < ntst; ++i) {
      htru = c1 * c2 * (k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta)) / (omg - omf); // Exact result
      httru =
         c1 * c2 * (k_it(0.0, omf, beta) * k_it(ttst(i), omg, beta) - k_it(ttst(i), omf, beta) * k_it(0.0, omg, beta)) / (omf - omg); // Exact result
      htst    = itops.coefs2eval(hc, ttst(i));
      httst   = itops.coefs2eval(htc, ttst(i));
      errlinf = std::max(errlinf, max_element(abs(htru - htst)));
      errl2 += pow(frobenius_norm(htru - htst), 2);
      errtlinf = std::max(errtlinf, max_element(abs(httru - httst)));
      errtl2 += pow(frobenius_norm(httru - httst), 2);
    }

    EXPECT_LT(errlinf, 13 * eps);
    EXPECT_LT(errl2, 2 * eps);
    EXPECT_LT(errtlinf, 13 * eps);
    EXPECT_LT(errtl2, 2 * eps);
    std::cout << fmt::format("Ordinary convolution: L^inf err = {:e}, L^2 err = {:e}\n", errlinf, errl2);
    std::cout << fmt::format("Time-ordered convolution: L^inf err = {:e}, L^2 err = {:e}\n", errtlinf, errtl2);
  }
}

/**
* @brief Test reflection of matrix-valued Green's function
*/
TEST(imtime_ops, refl_matrix) {

  double lambda = 10;    // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 10; // Inverse temperature
  int ntst    = 10; // # imag time test points

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  int norb = 2; // Orbital dimensions

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Sample Green's function G at DLR imaginary time nodes
  int r              = itops.rank();
  auto const &dlr_it = itops.get_itnodes();
  auto g             = nda::array<double, 3>(r, norb, norb);
  double om0 = 0.1234, om1 = -0.5678, om2 = 0.9012, om3 = -0.3456;
  for (int i = 0; i < r; ++i) {
    g(i, 0, 0) = k_it(dlr_it(i), om0, beta);
    g(i, 1, 0) = k_it(dlr_it(i), om1, beta);
    g(i, 0, 1) = k_it(dlr_it(i), om2, beta);
    g(i, 1, 1) = k_it(dlr_it(i), om3, beta);
  }

  auto gr  = itops.reflect(g);     // Reflection
  auto grc = itops.vals2coefs(gr); // DLR coefficients of reflection

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<double>(norb, norb);
  auto gtst      = nda::matrix<double>(norb, norb);
  double errlinf = 0, errl2 = 0, t = 0;
  for (int i = 0; i < ntst; ++i) { // Get evaluation point beta - t in relative format
    if (ttst(i) == 0) {
      t = 1;
    } else {
      t = -ttst(i);
    }
    gtru(0, 0) = k_it(t, om0, beta);
    gtru(1, 0) = k_it(t, om1, beta);
    gtru(0, 1) = k_it(t, om2, beta);
    gtru(1, 1) = k_it(t, om3, beta);

    gtst    = itops.coefs2eval(grc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test symmetrized DLR interpolation and evaluation for fermionic
* matrix-valued Green's function
*/
TEST(imtime_ops, interp_matrix_sym_fer) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps, SYM);
  int r       = dlr_rf.size();

  // Verify DLR rank is even
  EXPECT_EQ(r % 2, 0);

  // Verify symmetry
  EXPECT_EQ(max_element(abs(dlr_rf(range(r / 2)) + dlr_rf(range(r - 1, r / 2 - 1, -1)))), 0);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf, SYM);

  // Sample Green's function G at DLR imaginary time nodes
  auto const &dlr_it = itops.get_itnodes();

  // Verify symmetry
  EXPECT_EQ(max_element(abs(dlr_it(range(r / 2)) + dlr_it(range(r - 1, r / 2 - 1, -1)))), 0);

  auto g = nda::array<double, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<double>(norb, norb);
  auto gtst      = nda::matrix<double>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary frequency
  auto ifops = imfreq_ops(lambda, dlr_rf, Fermion);

  auto gtru_if = nda::matrix<dcomplex>(norb, norb);
  auto gtst_if = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru_if = gfun(norb, beta, n, Fermion);
    gtst_if = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru_if - gtst_if)));
    errl2 += pow(frobenius_norm(gtru_if - gtst_if), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test symmetrized DLR interpolation and evaluation for bosonic
* matrix-valued Green's function
*/
TEST(imtime_ops, interp_matrix_sym_bos) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 1000;  // Inverse temperature
  int ntst    = 10000; // # imag time test points
  int nmaxtst = 10000; // # imag freq test points

  int norb = 2; // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps, SYM);
  int r       = dlr_rf.size();

  // Verify DLR rank is even
  EXPECT_EQ(r % 2, 0);

  // Verify symmetry
  EXPECT_EQ(max_element(abs(dlr_rf(range(r / 2)) + dlr_rf(range(r - 1, r / 2 - 1, -1)))), 0);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf, SYM);

  // Obtain DLR imaginary time nodes
  auto const &dlr_it = itops.get_itnodes();

  // Verify symmetry
  EXPECT_EQ(max_element(abs(dlr_it(range(r / 2)) + dlr_it(range(r - 1, r / 2 - 1, -1)))), 0); // r even

  // Sample Green's function at DLR nodes
  auto g = nda::array<double, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_it(i)); }

  // DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Check that G can be recovered at imaginary time nodes
  EXPECT_LT(max_element(abs(itops.coefs2vals(gc) - g)), 1e-14);

  // Get test points in relative format
  auto ttst = eqptsrel(ntst);

  // Compute error in imaginary time
  auto gtru      = nda::matrix<double>(norb, norb);
  auto gtst      = nda::matrix<double>(norb, norb);
  double errlinf = 0, errl2 = 0;
  for (int i = 0; i < ntst; ++i) {
    gtru    = gfun(norb, beta, ttst(i));
    gtst    = itops.coefs2eval(gc, ttst(i));
    errlinf = std::max(errlinf, max_element(abs(gtru - gtst)));
    errl2 += pow(frobenius_norm(gtru - gtst), 2);
  }
  errl2 = sqrt(errl2 / ntst);

  EXPECT_LT(errlinf, 10 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag time: L^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);

  // Compute error in imaginary frequency
  auto ifops = imfreq_ops(lambda, dlr_rf, Fermion);

  auto gtru_if = nda::matrix<dcomplex>(norb, norb);
  auto gtst_if = nda::matrix<dcomplex>(norb, norb);
  errlinf = 0, errl2 = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru_if = gfun(norb, beta, n, Fermion);
    gtst_if = ifops.coefs2eval(beta, gc, n);
    errlinf = std::max(errlinf, max_element(abs(gtru_if - gtst_if)));
    errl2 += pow(frobenius_norm(gtru_if - gtst_if), 2);
  }
  errl2 = sqrt(errl2) / beta;

  EXPECT_LT(errlinf, 100 * eps);
  EXPECT_LT(errl2, eps);
  std::cout << fmt::format("Imag freq: l^2 err = {:e}, L^inf err = {:e}\n", errl2, errlinf);
}

/**
* @brief Test inner product of two DLR expansions
*/
TEST(imtime_ops, innerprod) {

  double lambda = 200;   // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  double beta = 100; // Inverse temperature
  int norb    = 2;   // Orbital dimensions

  std::cout << fmt::format("eps = {:e}, Lambda = {:e}\n", eps, lambda);

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);
  int r       = dlr_rf.size();

  // Get DLR imaginary time object
  auto itops  = imtime_ops(lambda, dlr_rf);
  auto dlr_it = itops.get_itnodes();

  // We take the inner product of two matrix-valued functions, each of the form
  //
  // c_{1,2} * exp(-tau * om_{1,2})/(1+exp(-beta * om_{1,2}))
  //
  // where c_1, c_2 are 2x2 complex-valued matrices, and om_1, om_2 are real
  // frequencies. The inner product is given by
  //
  // (sum_ij conj(c1_ij) * c2_ij) * (1+exp(-beta * om_1))^{-1} * (1+exp(-beta *
  // om_2))^{-1}  * (1 - exp(-beta * (om_1 + om_2))) / (om_1 + om_2)
  //
  // If om_1 + om_2 is positive and not too small, this formula is numerically
  // stable, and can be used for an analytical reference.

  auto c1    = nda::array<dcomplex, 2>(norb, norb);
  auto c2    = nda::array<dcomplex, 2>(norb, norb);
  c1(0, 0)   = 0.12 + 0.21i;
  c1(0, 1)   = 0.34 + 0.43i;
  c1(1, 1)   = 0.56 + 0.65i;
  c1(1, 0)   = conj(c1(0, 1));
  c2(0, 0)   = -0.22 + 0.11i;
  c2(0, 1)   = -0.44 + 0.33i;
  c2(1, 1)   = 0.66 - 0.55i;
  c2(1, 0)   = conj(c2(0, 1));
  double om1 = 0.0789;
  double om2 = 0.456;

  std::complex<double> iptrue = conj(c1(0, 0)) * c2(0, 0) + conj(c1(0, 1)) * c2(0, 1) + conj(c1(1, 0)) * c2(1, 0) + conj(c1(1, 1)) * c2(1, 1);
  iptrue *= (1.0 - exp(-beta * (om1 + om2))) / (beta * (om1 + om2) * (1.0 + exp(-beta * om1)) * (1.0 + exp(-beta * om2)));

  // Sample Green's functions at DLR nodes
  auto f = nda::array<dcomplex, 3>(r, norb, norb);
  auto g = nda::array<dcomplex, 3>(r, norb, norb);
  for (int i = 0; i < r; ++i) {
    f(i, _, _) = k_it(dlr_it(i), om1, beta) * c1;
    g(i, _, _) = k_it(dlr_it(i), om2, beta) * c2;
  }

  // DLR coefficients of f and g
  auto fc = itops.vals2coefs(f);
  auto gc = itops.vals2coefs(g);

  std::complex<double> iptest = itops.innerprod(fc, gc);

  double err = abs(iptest - iptrue);
  std::cout << fmt::format("True inner product: {:.16e}+1i*{:.16e}\n", iptrue.real(), iptrue.imag());
  std::cout << fmt::format("Test inner product: {:.16e}+1i*{:.16e}\n", iptest.real(), iptest.imag());
  std::cout << fmt::format("Inner product error: {:e}\n", err);

  EXPECT_LT(err, eps);
}

/**
* @brief Test transpose of DLR values -> coefficients transformation
*/
TEST(dlr_imtime, vals2coefs_transpose) {

  double lambda = 100;   // DLR cutoff
  double eps    = 1e-12; // DLR tolerance
  double beta   = 100;   // Inverse temperature

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);
  int r       = dlr_rf.size();

  // Get DLR coefficients -> values matrix
  auto itops         = imtime_ops(lambda, dlr_rf);
  auto const &dlr_it = itops.get_itnodes();

  // Define Green's function G(t) = k(t,om0) for fixed om0
  double om0 = 0.4321;
  auto g     = nda::vector<double>(r);
  for (int k = 0; k < r; ++k) { g(k) = k_it(dlr_it(k), beta * om0); }
  auto gc = itops.vals2coefs(g);

  // Row vector of integration taking DLR coefficients
  auto intvec = nda::vector<double>(r);
  for (int k = 0; k < r; ++k) { intvec(k) = -beta * tanh(dlr_rf(k) / 2) / dlr_rf(k); }

  // Compute integral of G via integration vector and compare with true integral
  double inttst = nda::blas::dot(intvec, gc);
  double inttru = -tanh(beta * om0 / 2) / om0;

  EXPECT_LT(abs(inttst - inttru), 5 * eps);

  // Build row vector of integration taking values at DLR nodes; since
  //
  // intvec^T * gc = intvec^T * it2cf * g = (it2cf^T * intvec)^T * g,
  //
  // this is given by applying the transpose of the DLR values -> coefficients
  // transformation to the integration vector
  auto intvec_vals = itops.vals2coefs(intvec, true);

  // Compute integral of G using this vector and compare with true integral
  double inttst2 = nda::blas::dot(intvec_vals, g);

  EXPECT_LT(abs(inttst2 - inttru), 5 * eps);
}

TEST(dlr_imtime, h5_rw) {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  auto filename = "data_imtime_ops_h5_rw.h5";
  auto name     = "itops";

  {
    h5::file file(filename, 'w');
    h5::write(file, name, itops);
  }

  imtime_ops itops_ref;
  {
    h5::file file(filename, 'r');
    h5::read(file, name, itops_ref);
  }

  // Check equal
  EXPECT_EQ(itops.lambda(), itops_ref.lambda());
  EXPECT_EQ(itops.rank(), itops_ref.rank());
  EXPECT_EQ_ARRAY(itops.get_rfnodes(), itops_ref.get_rfnodes());
  EXPECT_EQ_ARRAY(itops.get_itnodes(), itops_ref.get_itnodes());
  EXPECT_EQ_ARRAY(itops.get_cf2it(), itops_ref.get_cf2it());
  EXPECT_EQ_ARRAY(itops.get_it2cf_lu(), itops_ref.get_it2cf_lu());
  EXPECT_EQ_ARRAY(itops.get_it2cf_zlu(), itops_ref.get_it2cf_zlu());
  EXPECT_EQ_ARRAY(itops.get_it2cf_piv(), itops_ref.get_it2cf_piv());
}
