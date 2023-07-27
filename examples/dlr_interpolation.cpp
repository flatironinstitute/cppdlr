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

#include <cppdlr/cppdlr.hpp>
#include <filesystem> // For outputting data
#include <fstream>    // For outputting data

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
  // here by the function k_it.

  g(0, 0) = 0.9 * k_it(t, 0.8, beta) + 0.7 * k_it(t, -0.3, beta);
  g(0, 1) = 0.3 * k_it(t, 0.4, beta) + 0.3 * k_it(t, 0.7, beta);
  g(1, 0) = g(0, 1);
  g(1, 1) = 0.7 * k_it(t, -0.5, beta);

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
  // (2n+1) * i * pi / beta, is implemented here by the function k_if.

  g(0, 0) = 0.9 * k_if(n, 0.8, Fermion, beta) + 0.7 * k_if(n, -0.3, Fermion, beta);
  g(0, 1) = 0.3 * k_if(n, 0.4, Fermion, beta) + 0.3 * k_if(n, 0.7, Fermion, beta);
  g(1, 0) = g(0, 1);
  g(1, 1) = 0.7 * k_if(n, -0.5, Fermion, beta);

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

  // Output data?
  bool output = true;

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
  auto gtru_t = nda::array<double, 3>(ntst_t, norb, norb);
  auto gtst_t = nda::array<double, 3>(ntst_t, norb, norb);
  for (int i = 0; i < ntst_t; ++i) {
    gtru_t(i, _, _) = gfun(beta, ttst(i));           // Evaluate true Green's function
    gtst_t(i, _, _) = itops.coefs2eval(gc, ttst(i)); // Evaluate DLR expansion
  }
  std::cout << "Maximum error of DLR expansion in imaginary time: " << max_element(abs(gtru_t - gtst_t)) << std::endl;

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, Fermion);

  // Compute maximum error of Green's function on test points in Matsubara frequency
  auto gtru_om = nda::array<dcomplex, 3>(nmax_om, norb, norb);
  auto gtst_om = nda::array<dcomplex, 3>(nmax_om, norb, norb);
  for (int n = -nmax_om / 2; n < nmax_om / 2; ++n) {
    gtru_om(n + nmax_om / 2, _, _) = gfun(beta, n);                 // Evaluate true Green's function
    gtst_om(n + nmax_om / 2, _, _) = ifops.coefs2eval(beta, gc, n); // Evaluate DLR expansion
  }
  std::cout << "Maximum error of DLR expansion in imaginary frequency: " << max_element(abs(gtru_om - gtst_om)) << std::endl;

  // ---------------------------------------------------------------------------

  // Output data
  if (output) {

    std::cout << "Writing results to `data` directory..." << std::endl;

    std::filesystem::create_directory("data"); // Create directory for data

    // DLR frequencies
    std::ofstream dlr_rf_file("data/dlr_rf");
    dlr_rf_file << std::setprecision(16);
    for (int i = 0; i < r; ++i) { dlr_rf_file << dlr_rf(i) << std::endl; }
    dlr_rf_file.close();

    // DLR imaginary time nodes and Green's function values at DLR nodes
    auto dlr_it_abs = rel2abs(dlr_it); // Convert DLR imaginary time nodes from relative to absolute time format
    std::ofstream dlr_it_file("data/dlr_it");
    std::ofstream g_dlr_it_file("data/g_dlr_it");
    dlr_it_file << std::setprecision(16);
    g_dlr_it_file << std::setprecision(16);
    for (int i = 0; i < r; ++i) {
      dlr_it_file << dlr_it_abs(i) << std::endl;
      for (int j = 0; j < norb; ++j) {
        for (int k = 0; k < norb; ++k) { g_dlr_it_file << g(i, j, k) << std::endl; }
      }
    }
    dlr_it_file.close();
    g_dlr_it_file.close();

    // Imaginary time data
    auto ttst_abs = rel2abs(ttst); // Convert test points from relative to absolute time format
    std::ofstream ttst_file("data/ttst");
    std::ofstream gtru_t_file("data/gtru_t");
    std::ofstream gtst_t_file("data/gtst_t");
    ttst_file << std::setprecision(16);
    gtru_t_file << std::setprecision(16);
    gtst_t_file << std::setprecision(16);
    for (int i = 0; i < ntst_t; ++i) {
      ttst_file << ttst_abs(i) << std::endl;
      for (int j = 0; j < norb; ++j) {
        for (int k = 0; k < norb; ++k) {
          gtru_t_file << gtru_t(i, j, k) << std::endl;
          gtst_t_file << gtst_t(i, j, k) << std::endl;
        }
      }
    }
    ttst_file.close();
    gtru_t_file.close();
    gtst_t_file.close();

    // Matsubara frequency data
    std::ofstream n_file("data/n");
    std::ofstream gtru_om_file("data/gtru_om");
    std::ofstream gtst_om_file("data/gtst_om");
    gtru_om_file << std::setprecision(16);
    gtst_om_file << std::setprecision(16);
    for (int n = -nmax_om / 2; n < nmax_om / 2; ++n) {
      n_file << n << std::endl;
      for (int j = 0; j < norb; ++j) {
        for (int k = 0; k < norb; ++k) {
          gtru_om_file << real(gtru_om(n + nmax_om / 2, j, k)) << " " << imag(gtru_om(n + nmax_om / 2, j, k)) << std::endl;
          gtst_om_file << real(gtst_om(n + nmax_om / 2, j, k)) << " " << imag(gtst_om(n + nmax_om / 2, j, k)) << std::endl;
        }
      }
    }
    n_file.close();
    gtru_om_file.close();
    gtst_om_file.close();
  }
}