#include <gtest/gtest.h>
#include <cppdlr/utils.hpp>
#include <cmath>
#include <nda/blas.hpp>

using namespace cppdlr;

// Test pivoted reorthogonalized Gram-Schmidt

TEST(pivrgs, pivrgs) {

  auto _ = nda::range::all;

  int m      = 40;
  int n      = 50;
  double eps = 1e-6;

  // Generate random numerically low rank mxn matrix: rapidly decaying singular
  // values 2^{-k}

  // Get random matrices

  auto a1 = nda::matrix<double>::rand(m, m);
  auto a2 = nda::matrix<double>::rand(n, n);

  // Orthonormalize to get random orthogonal matrices U, V. Note: this also
  // tests whether behavior of pivrgs is correct for numerically full rank
  // matrix.

  auto [u, norms1] = pivrgs(a1, 1e-100);
  auto [v, norms2] = pivrgs(a2, 1e-100);

  // [Q] How to do this more concisely?

  EXPECT_EQ(u.shape(0), m);
  EXPECT_EQ(u.shape(1), m);
  EXPECT_EQ(v.shape(0), n);
  EXPECT_EQ(v.shape(1), n);

  // Check u and v are orthogonal as well

  EXPECT_LE(frobenius_norm(nda::eye<double>(m) - transpose(u) * u), 1e-14);
  EXPECT_LE(frobenius_norm(nda::eye<double>(n) - transpose(v) * v), 1e-14);

  // Multiply first matrix by singular values S

  for (int i = 0; i < m; ++i) { u(_, i) *= pow(2.0, -i); }

  // Matrix A = U*S*V has given singular values

  auto a = u * v(nda::range(0, m), _);

  // Gram-Schmidt to obtain orthonormal basis Q of column space of A

  auto [q, norms] = pivrgs(a, eps);
  int r           = norms.size(); // Estimated rank

  // Verify rank is almost correct

  EXPECT_LE(r, ceil(log2(1.0 / eps)) + 3);
  EXPECT_GE(r, ceil(log2(1.0 / eps)) - 3);

  // Verify columns of Q are orthonormal to near double precision

  EXPECT_LE(frobenius_norm(nda::eye<double>(r) - transpose(q) * q), 1e-14);

  // Verify column space of A contained in that of Q; make sure projection
  // onto col(Q) of random linear combination of columns of A yields identity to
  // within target accuracy

  auto x = nda::vector<double>::rand(n);
  x      = 2 * x - 1;
  x /= sqrt(nda::blas::dot(x, x));
  auto b = a * x;

  auto tmp = b-q*(transpose(q)*b);

  EXPECT_LT(sqrt(nda::blas::dot(tmp, tmp)), 10 * eps);

  // EXPECT_LT(frobenius_norm(a - q * transpose(q) * a), 10 * eps); // More stringent test...
  //
}
