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
#include <cppdlr/utils.hpp>
#include <cmath>
#include <vector>
#include <nda/blas.hpp>

using namespace cppdlr;
using namespace nda;

/**
 * Test pivoted reorthogonalized Gram-Schmidt function pivrgs for real-valued
 * matrix
 */
TEST(pivrgs, pivrgs_real) {

  // Matrix size and rank cutoff tolerance
  int m      = 50;
  int n      = 40;
  double eps = 1e-6;

  // Generate random numerically low rank mxn matrix: rapidly decaying singular
  // values 2^{-k}

  // Get random matrices
  auto a1 = nda::matrix<double>::rand(m, m);
  auto a2 = nda::matrix<double>::rand(n, n);

  // Orthonormalize to get random orthogonal matrices U, V. Note: this also
  // tests whether behavior of pivrgs is correct for numerically full rank
  // matrix.
  auto [u, norms1, piv1] = pivrgs(a1, 1e-100);
  auto [v, norms2, piv2] = pivrgs(a2, 1e-100);

  EXPECT_EQ(u.shape(), (std::array<long, 2>{m, m}));
  EXPECT_EQ(v.shape(), (std::array<long, 2>{n, n}));

  // Check u and v are orthogonal as well
  EXPECT_LE(frobenius_norm(eye<double>(m) - transpose(u) * u), 1e-14);
  EXPECT_LE(frobenius_norm(eye<double>(n) - transpose(v) * v), 1e-14);

  // Multiply smaller matrix by singular values S
  for (int i = 0; i < n; ++i) { v(i, _) *= pow(2.0, -i); }

  // Matrix A = U*S*V has given singular values
  auto a = u(_, range(n)) * v;

  // Gram-Schmidt to obtain orthonormal basis Q of row space of A
  auto [q, norms, piv] = pivrgs(a, eps);
  int r                = norms.size(); // Estimated rank

  // Verify rank is almost correct
  EXPECT_LE(r, ceil(log2(1.0 / eps)) + 3);
  EXPECT_GE(r, ceil(log2(1.0 / eps)) - 3);

  // Verify rows of Q are orthonormal to near double precision
  EXPECT_LE(frobenius_norm(eye<double>(r) - q * transpose(q)), 1e-14);

  // Verify row space of A contained in that of Q; make sure projection onto
  // row(Q) of random linear combination of rows of A yields identity to
  // within roughly target accuracy.
  auto x = nda::vector<double>::rand(m);
  x      = 2 * x - 1;
  x /= sqrt(blas::dot(x, x));
  auto b = transpose(a) * x;

  auto tmp = make_regular(b - transpose(q) * (q * b));

  EXPECT_LT(sqrt(blas::dot(tmp, tmp)), 10 * eps);

  // More stringent test: make sure projection of A onto row space of Q is A
  // to within roughly target accuracy.
  EXPECT_LT(frobenius_norm(a - (a * transpose(q)) * q), 10 * eps);

  // Verify that calling pivoted GS on matrix of pivot rows yields an
  // identical Q matrix, and the correct order of pivots; this tests that the
  // function returns the correct pivots
  auto athin = nda::matrix<double>(r, n);
  for (int i = 0; i < r; ++i) { athin(i, _) = a(piv(i), _); }

  auto [qthin, normthin, pivthin] = pivrgs(athin, eps);

  EXPECT_EQ(pivthin.size(), r);
  EXPECT_EQ(pivthin, arange(r));
  EXPECT_LE(frobenius_norm(q - qthin), 1e-14);
}

/**
 * Test that the eps-cutoff symmetrized Gram-Schmidt pivrgs_sym handles an odd
 * number of rows: no throw, the central row (m-1)/2 is selected as a
 * self-paired pivot giving an odd rank, the remaining pivots come in mirror
 * pairs, and the returned rows are orthonormal.
 */
TEST(pivrgs, pivrgs_sym_odd) {

  int m      = 51; // Odd number of rows
  int n      = 40;
  double eps = 1e-6;

  // Random matrix with rapidly decaying column norms, so numerically low-rank
  auto a = nda::matrix<double>::rand(m, n);
  for (int j = 0; j < n; ++j) { a(_, j) *= pow(2.0, -j); }

  // Symmetrized pivoted Gram-Schmidt with eps cutoff (must not throw for odd m)
  auto [q, norms, piv] = pivrgs_sym(a, eps);
  int r                = norms.size();

  // Rank is odd, and below the number of columns
  EXPECT_EQ(r % 2, 1);
  EXPECT_LT(r, n);

  // The central row (m-1)/2, its own mirror, is selected as the first pivot
  EXPECT_EQ(piv(0), (m - 1) / 2);

  // The remaining pivots come in mirror pairs: i selected => m-1-i selected
  auto selected = std::vector<bool>(m, false);
  for (auto p : piv) { selected[p] = true; }
  for (auto p : piv) { EXPECT_TRUE(selected[m - 1 - p]); }

  // Returned rows are orthonormal to near double precision
  EXPECT_LE(frobenius_norm(eye<double>(r) - q * transpose(q)), 1e-13);
}

/**
 * Test that the eps-cutoff pivrgs_sym rejects an odd-row matrix whose middle row
 * is negligible. That row is a forced pivot rather than one chosen by norm, so
 * there is no symmetrized basis of odd rank, however large the other rows are.
 */
TEST(pivrgs, pivrgs_sym_negligible_middle_row) {

  int m      = 51; // Odd number of rows
  int n      = 40;
  double eps = 1e-6;

  auto a = nda::matrix<double>::rand(m, n);
  for (int j = 0; j < n; ++j) { a(_, j) *= pow(2.0, -j); }
  EXPECT_NO_THROW(pivrgs_sym(a, eps));

  a((m - 1) / 2, _) = 0;
  EXPECT_THROW(pivrgs_sym(a, eps), std::runtime_error);
}

/**
 * Test pivoted reorthogonalized Gram-Schmidt function pivrgs for complex-valued
 * matrix
 */
TEST(pivrgs, pivrgs_cplx) {

  // Matrix size and rank cutoff tolerance
  int m      = 50;
  int n      = 40;
  double eps = 1e-6;

  // Generate random numerically low rank mxn matrix: rapidly decaying singular
  // values 2^{-k}

  // Get random matrices
  auto a1 = nda::matrix<dcomplex>::rand(m, m);
  auto a2 = nda::matrix<dcomplex>::rand(n, n);

  // Orthonormalize to get random orthogonal matrices U, V. Note: this also
  // tests whether behavior of pivrgs is correct for numerically full rank
  // matrix.
  auto [u, norms1, piv1] = pivrgs(a1, 1e-100);
  auto [v, norms2, piv2] = pivrgs(a2, 1e-100);

  EXPECT_EQ(u.shape(), (std::array<long, 2>{m, m}));
  EXPECT_EQ(v.shape(), (std::array<long, 2>{n, n}));

  // Check u and v are orthogonal as well
  EXPECT_LE(frobenius_norm(eye<double>(m) - conj(transpose(u)) * u), 1e-14);
  EXPECT_LE(frobenius_norm(eye<double>(n) - conj(transpose(v)) * v), 1e-14);

  // Multiply smaller matrix by singular values S
  for (int i = 0; i < n; ++i) { v(i, _) *= pow(2.0, -i); }

  // Matrix A = U*S*V has given singular values
  auto a = u(_, range(n)) * v;

  // Gram-Schmidt to obtain orthonormal basis Q of row space of A
  auto [q, norms, piv] = pivrgs(a, eps);
  int r                = norms.size(); // Estimated rank

  // Verify rank is almost correct
  EXPECT_LE(r, ceil(log2(1.0 / eps)) + 3);
  EXPECT_GE(r, ceil(log2(1.0 / eps)) - 3);

  // Verify rows of Q are orthonormal to near double precision
  EXPECT_LE(frobenius_norm(eye<double>(r) - q * conj(transpose(q))), 1e-14);

  // Verify row space of A contained in that of Q; make sure projection onto
  // row(Q) of random linear combination of rows of A yields identity to
  // within roughly target accuracy.
  auto x = nda::vector<dcomplex>::rand(m);
  x      = 2 * x - 1;
  x /= sqrt(real(blas::dotc(x, x)));
  auto b = conj(transpose(a)) * x;

  auto tmp = make_regular(b - conj(transpose(q)) * (q * b));

  EXPECT_LT(sqrt(real(blas::dotc(tmp, tmp))), 10 * eps);

  // More stringent test: make sure projection of A onto row space of Q is A
  // to within roughly target accuracy.
  EXPECT_LT(frobenius_norm(a - (a * conj(transpose(q))) * q), 10 * eps);

  // Verify that calling pivoted GS on matrix of pivot rows yields an
  // identical Q matrix, and the correct order of pivots; this tests that the
  // function returns the correct pivots
  auto athin = nda::matrix<dcomplex>(r, n);
  for (int i = 0; i < r; ++i) { athin(i, _) = a(piv(i), _); }

  auto [qthin, normthin, pivthin] = pivrgs(athin, eps);

  EXPECT_EQ(pivthin.size(), r);
  EXPECT_EQ(pivthin, arange(r));
  EXPECT_LE(frobenius_norm(q - qthin), 1e-14);
}
