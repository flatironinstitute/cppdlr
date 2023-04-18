#include <gtest/gtest.h>
#include <cppdlr/utils.hpp>

using namespace cppdlr;

/**
 * Test arraymult function: double matrix times complex array
 */
TEST(arraymult, matrix_array) {

  int m = 3;
  int n = 4;
  int p = 5;
  int q = 6;

  auto a = nda::matrix<double>(m, n);
  auto b = nda::array<dcomplex, 3>(n, p, q);

  // Populate with random entries
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) { a(i, j) = sin(10000.0 * (i + j)); }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < q; ++k) { b(i, j, k) = dcomplex(sin(10000.0 * (i + j + k)), cos(10000.0 * (i + j + k))); }
    }
  }

  // Get correct answer
  auto ctrue = nda::array<dcomplex, 3>(m, p, q);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < q; ++k) {
        ctrue(i, j, k) = 0.0;
        for (int l = 0; l < n; ++l) { ctrue(i, j, k) += a(i, l) * b(l, j, k); }
      }
    }
  }

  // Get arraymult result and compare
  auto c = arraymult(a, b);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < q; ++k) { EXPECT_LT(abs(ctrue(i, j, k) - c(i, j, k)), 1e-15); }
    }
  }
}

/**
 * Test arraymult function: complex array times complex array
 */
TEST(arraymult, array_array) {

  int m = 3;
  int n = 4;
  int p = 5;
  int q = 6;

  auto a = nda::array<dcomplex, 3>(m, n, p);
  auto b = nda::array<dcomplex, 2>(p, q);

  // Populate with random entries
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < p; ++k) { a(i, j, k) = dcomplex(sin(10000.0 * (i + j + k)), cos(10000.0 * (i + j + k))); }
    }
  }

  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < q; ++j) { b(i, j) = dcomplex(sin(10000.0 * (i + j)), cos(10000.0 * (i + j))); }
  }

  // Get correct answer
  auto ctrue = nda::array<dcomplex, 3>(m, n, q);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < q; ++k) {
        ctrue(i, j, k) = 0.0;
        for (int l = 0; l < p; ++l) { ctrue(i, j, k) += a(i, j, l) * b(l, k); }
      }
    }
  }

  // Get arraymult result and compare
  auto c = arraymult(a, b);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < q; ++k) { EXPECT_LT(abs(ctrue(i, j, k) - c(i, j, k)), 1e-15); }
    }
  }
}

/**
 * Test arraymult function: complex matrix times double vector
 */
TEST(arraymult, matrix_vector) {

  int m = 3;
  int n = 4;

  auto a = nda::matrix<dcomplex>(m, n);
  auto b = nda::vector<double>(n);

  // Populate with random entries
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) { a(i, j) = dcomplex(sin(10000.0 * (i + j)), cos(10000.0 * (i + j))); }
  }

  for (int i = 0; i < n; ++i) { b(i) = sin(10000.0 * i); }

  // Get correct answer
  auto ctrue = nda::array<dcomplex, 1>(m);
  for (int i = 0; i < m; ++i) {
    ctrue(i) = 0.0;
    for (int j = 0; j < n; ++j) { ctrue(i) += a(i, j) * b(j); }
  }

  // Get arraymult result and compare
  auto c = arraymult(a, b);

  for (int i = 0; i < m; ++i) { EXPECT_LT(abs(ctrue(i) - c(i)), 1e-15); }
}

/**
 * Test arraymult function: double vector times double array
 */
TEST(arraymult, vector_array) {

  int m = 3;
  int n = 4;
  int p = 5;

  auto a = nda::vector<double>(m);
  auto b = nda::array<double, 3>(m, n, p);

  // Populate with random entries
  for (int i = 0; i < m; ++i) { a(i) = sin(10000.0 * i); }

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < p; ++k) { b(i, j, k) = sin(10000.0 * (i + j + k)); }
    }
  }

  // Get correct answer
  auto ctrue = nda::array<double, 2>(n, p);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      ctrue(i, j) = 0.0;
      for (int k = 0; k < m; ++k) { ctrue(i, j) += a(k) * b(k, i, j); }
    }
  }

  // Get arraymult result and compare
  auto c = arraymult(a, b);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) { EXPECT_LT(abs(ctrue(i, j) - c(i, j)), 1e-15); }
  }
}