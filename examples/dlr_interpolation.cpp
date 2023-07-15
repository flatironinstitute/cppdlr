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
// Authors: Jason Kaye

/** 
* @file dlr_interpolation.cpp
*
* @brief Example of DLR interpolation from imaginary time data
* The Green's function is taken to be a sum of a few exponentials. The DLR
* expansion of the Green's function is formed from its samples at the DLR
* imaginary time nodes. We test that the DLR expansion is correct on dense grids
* in imaginary time and imaginary frequency.
*/

#include <nda/nda.hpp>
#include <cppdlr/cppdlr.hpp>

using namespace cppdlr;

/**
* @brief Fermionic 2x2 matrix-valued Green's function which is a sum of a few exponentials
* in each entry
*
* The spectral width (maximum non-zero frequency of spectral function) is less
* than 1.
*
* @param beta Inverse temperature
* @param t    Imaginary time evaluation point (in relative time format)
*
* @return Green's function evaluated at t
*/
nda::matrix<double> gfun(double beta, double t) {

  auto g = nda::matrix<double>(2, 2); // Green's function represented by nda matrix

  // The kernel K(t, om) = -exp(-om * t) / (1 + exp(-beta * om)) is implemented
  // here by the function k_it, which uses dimenionless variables
  // (corresponding to beta = 1). Thus, to evaluate K(t, om), we scale the
  // frequency argument of k_it by hand. See the documentation of k_it for
  // more details.

  g(0, 0) = 0.9 * k_it(t, beta * 0.8) + 0.1 * k_it(t, beta * 0.3);
  g(0, 1) = 0.3 * k_it(t, beta * 0.4) + 0.3 * k_it(t, beta * 0.7);
  g(1, 0) = g(0, 1);
  g(1, 1) = 0.7 * k_it(t, beta * 0.5);

  return g;
}

/**
* @brief Fermionic 2x2 matrix-valued Green's function which is a sum of a few poles
* in each entry (the fermionic Fourier transform of the Green's function defined
* in imaginary time above)
*
* The spectral width (maximum non-zero frequency of spectral function) is less
* than 1.
*
* @param beta Inverse temperature
* @param n    Index of fermionic Matsubara frequency evaluation point
*
* @return Green's function evaluated at Matsubara frequency i nu_n
*/
nda::matrix<dcomplex> gfun(double beta, int n) {

  auto g = nda::matrix<dcomplex>(2, 2); // Green's function represented by nda matrix

  // The fermionic kernel K(i nu_n, om) = 1 / (i nu_n - om), with i nu_n =
  // (2n+1) * i * pi / beta, is implemented here by the function k_if,  which
  // uses dimenionless variables (corresponding to beta = 1). Thus, to
  // evaluate K(i nu_n, om), we scale the frequency argument of k_if by hand,
  // and multiply the result by beta. Please see the documentation of k_if for
  // more details.

  g(0, 0) = 0.9 * beta * k_if(n, beta * 0.8, Fermion) + 0.1 * beta * k_if(n, beta * 0.3, Fermion);
  g(0, 1) = 0.3 * beta * k_if(n, beta * 0.4, Fermion) + 0.3 * beta * k_if(n, beta * 0.7, Fermion);
  g(1, 0) = g(0, 1);
  g(1, 1) = 0.7 * beta * k_if(n, beta * 0.5, Fermion);

  return g;
}

/**
* @brief Example of DLR interpolation from imaginary time values
*/
int main() {

  double beta = 100; // Inverse temperature
  int norb    = 2;   // Orbital dimensions

  // DLR cutoff Lambda = beta * omega_max. Here, we know omega_max (spectral
  // width of Green's function) is less than 1, so Lambda = beta is sufficient.
  double lambda = 100;

  // DLR error tolerance
  double eps = 1e-10;

  // # test points
  int ntst_t  = 1000; // In imaginary time
  int nmax_om = 1000; // In imaginary frequency

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);
  int r       = dlr_rf.size(); // # DLR basis functions

  std::cout << "DLR cutoff Lambda = " << lambda << std::endl;
  std::cout << "DLR tolerance epsilon = " << eps << std::endl;
  std::cout << "# DLR basis functions = " << r << std::endl;

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // Evaluate Green's function at DLR imaginary time interpolation nodes
  auto dlr_it = itops.get_itnodes();                  // r DLR imaginary time nodes
  auto g      = nda::array<double, 3>(r, norb, norb); // Green's function represented by ndarray
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(beta, dlr_it(i)); }

  // Obtain DLR coefficients of G
  auto gc = itops.vals2coefs(g);

  // Get dense equispaced grid of imaginary time test points in relative
  // format
  auto ttst = eqptsrel(ntst_t);

  // Compute maximum error of Green's function on test points in imaginary time
  auto gtru_t = nda::matrix<double>(norb, norb);
  auto gtst_t = nda::matrix<double>(norb, norb);
  double err  = 0;
  for (int i = 0; i < ntst_t; ++i) {
    gtru_t = gfun(beta, ttst(i));                              // Evaluate true Green's function
    gtst_t = itops.coefs2eval(gc, ttst(i));                    // Evaluate DLR expansion
    err    = std::max(err, max_element(abs(gtru_t - gtst_t))); // Measure error
  }
  std::cout << "Maximum error of DLR expansion in imaginary time: " << err << std::endl;

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, Fermion);

  // Compute maximum error of Green's function on test points in Matsubara frequency
  auto gtru_om = nda::matrix<dcomplex>(norb, norb);
  auto gtst_om = nda::matrix<dcomplex>(norb, norb);
  err          = 0;
  for (int n = -nmax_om / 2; n < nmax_om / 2; ++n) {
    gtru_om = gfun(beta, n);                                      // Evaluate true Green's function
    gtst_om = ifops.coefs2eval(beta, gc, n);                      // Evaluate DLR expansion
    err     = std::max(err, max_element(abs(gtru_om - gtst_om))); // Measure error
  }
  std::cout << "Maximum error of DLR expansion in imaginary frequency: " << err << std::endl;
}