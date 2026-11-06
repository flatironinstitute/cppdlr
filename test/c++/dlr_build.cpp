// Copyright (c) 2022-2023 Simons Foundation
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

#include <gtest/gtest.h>
#include <cppdlr/dlr_build.hpp>
#include <cppdlr/dlr_kernels.hpp>
#include <fmt/format.h>
#include <cmath>
#include <complex>
#include <numbers>

using namespace cppdlr;

TEST(dlr_build, fineparams) {

  fineparams fine1(100.0);

  // Check expected default values
  EXPECT_EQ(fine1.p, 24);
  EXPECT_EQ(fine1.npom, 7);
  EXPECT_EQ(fine1.npt, 5);
  EXPECT_EQ(fine1.nt, 2 * fine1.p * fine1.npt);
  EXPECT_EQ(fine1.nom, 2 * fine1.p * fine1.npom);

  // Check non-default value of p
  fineparams fine2(100.0, 7);
  EXPECT_EQ(fine2.p, 7);

  // Check error triggered for non-permissible parameters.
  EXPECT_THROW({ fineparams fine3(0.0); }, std::runtime_error);
  EXPECT_THROW({ fineparams fine3(1.0, 0); }, std::runtime_error);
}

/** 
* @brief Test that matrix of analytic continuation kernel in imaginary time is
* correct to nearly machine precision
*/
TEST(dlr_build, get_kfine) {

  fineparams fine(100.0);

  auto [t, w] = build_it_fine(fine);
  auto om     = build_rf_fine(fine);

  auto kmat          = build_k_it(t, om);
  auto [errt, errom] = geterr_k_it(fine, t, om, kmat);

  EXPECT_LT(errt, 1e-14);
  EXPECT_LT(errom, 1e-14);

  std::cout << fmt::format("Max imag time err = {:e}, Max freq err = {:e}\n", errt, errom);
}

/**
* @brief Test that the symmetric DLR construction includes the self-symmetric fixed
* points, omega=0 and tau=beta/2, and leaves the non-symmetric grids unchanged.
*/
TEST(dlr_build, symmetric_fixed_points) {

  double eps = 1e-10;

  for (double lambda : {10.0, 100.0, 1000.0}) {
    auto dlr_rf = build_dlr_rf(lambda, eps, SYM);
    long r      = dlr_rf.size();

    // Symmetric rank is odd, with omega=0 as the central self-paired node
    EXPECT_EQ(r % 2, 1);
    EXPECT_DOUBLE_EQ(dlr_rf((r - 1) / 2), 0.0);

    // Grid is mirror-symmetric about omega=0: omega_j = -omega_{r-1-j}, exactly, since
    // the negative half of the fine grid is built by negating the positive half
    for (int j = 0; j < r / 2; ++j) { EXPECT_EQ(dlr_rf(j), -dlr_rf(r - 1 - j)); }
  }

  // Symmetric fine grids: omega=0 and tau=beta/2 (relative t=0.5) are present as
  // central nodes, with one extra point relative to the even default grid.
  fineparams fine(100.0);

  auto om_sym = build_rf_fine(fine, SYM);
  EXPECT_EQ(om_sym.size(), fine.nom + 1);
  EXPECT_DOUBLE_EQ(om_sym(fine.nom / 2), 0.0);

  auto [t_sym, w_sym] = build_it_fine(fine, SYM);
  EXPECT_EQ(t_sym.size(), fine.nt + 1);
  EXPECT_DOUBLE_EQ(t_sym(fine.nt / 2), 0.5);
  EXPECT_DOUBLE_EQ(w_sym(fine.nt / 2), 0.0);

  // Non-symmetric grids are unchanged (even count, no central fixed-point node)
  auto [t_non, w_non] = build_it_fine(fine);
  EXPECT_EQ(build_rf_fine(fine).size(), fine.nom);
  EXPECT_EQ(t_non.size(), fine.nt);
  EXPECT_EQ(w_non.size(), fine.nt);
}

/**
* @brief Test structure of dyadic fine imaginary frequency grid
*/
TEST(dlr_build, build_if_fine) {

  for (auto statistic : {Fermion, Boson}) {
    for (int nmax : {1, 20, 1000, 1000000}) {
      for (int npt : {1, 7, 128}) {

        auto n = build_if_fine(nmax, statistic, npt);
        int m  = n.size();

        // Covers cutoff range, sorted, symmetric in the sense of pivrgs_sym
        int shift = (statistic == Fermion ? 1 : 0);
        EXPECT_EQ(n(0), -nmax);
        EXPECT_EQ(n(m - 1), nmax - shift);
        EXPECT_EQ(m % 2, 1 - shift);
        for (int i = 0; i + 1 < m; ++i) { EXPECT_LT(n(i), n(i + 1)); }
        for (int i = 0; i < m; ++i) { EXPECT_EQ(n(i) + n(m - 1 - i), -shift); }

        // Spacing on the half with n >= 1 - shift doubles every npt points, except for
        // the last step to the cutoff
        int i0 = (m + 1) / 2;
        for (int i = i0; i + 2 < m; ++i) { EXPECT_EQ(n(i + 1) - n(i), 1 << ((i - i0) / npt)); }
        if (m - i0 >= 2) { EXPECT_LE(n(m - 1) - n(m - 2), 1 << ((m - 2 - i0) / npt)); }

        // Full grid if panels cover cutoff range
        if (npt >= nmax) { EXPECT_EQ(m, 2 * nmax + 1 - shift); }
      }
    }
  }

  EXPECT_THROW(build_if_fine(0, Fermion, 128), std::runtime_error);
  EXPECT_THROW(build_if_fine(100, Fermion, 0), std::runtime_error);
}

/**
* @brief Test that kernel on given imaginary frequency indices matches kernel
* on full grid
*/
TEST(dlr_build, build_k_if_grid) {

  int nmax = 300;
  auto om  = nda::vector<double>{-250.0, -1.0, 0.0, 0.5, 300.0};

  for (auto statistic : {Fermion, Boson}) {
    auto n     = build_if_fine(nmax, statistic, 16);
    auto kfull = build_k_if(nmax, om, statistic);
    auto kgrid = build_k_if(n, om, statistic);
    for (int i = 0; i < n.size(); ++i) { EXPECT_EQ(nda::max_element(nda::abs(kgrid(i, nda::range::all) - kfull(n(i) + nmax, nda::range::all))), 0); }
  }
}

/**
* @brief Test analytic continuation kernel in imaginary frequency for indices
* where 2n+1 overflows int
*/
TEST(dlr_build, k_if_large_n) {

  using std::numbers::pi;
  for (int n : {(1 << 30) + 5, -(1 << 30) - 7, 2147483000, -2147483000}) {
    double om = 3.0;
    auto kf   = 1.0 / std::complex<double>(-om, (2.0 * n + 1) * pi);
    auto kb   = std::tanh(0.5 * om) / std::complex<double>(-om, 2.0 * n * pi);
    EXPECT_LE(std::abs(k_if(n, om, Fermion) - kf), 1e-15 * std::abs(kf));
    EXPECT_LE(std::abs(k_if(n, om, Boson) - kb), 1e-15 * std::abs(kb));
  }
}
