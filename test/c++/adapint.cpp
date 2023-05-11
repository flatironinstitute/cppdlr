// Copyright (c) 2023 Simons Foundation
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
// Authors: Nils Wentzell, jasonkaye

/**
* @file adapint.cpp 
*
*@brief Test adaptive integration. 
*/

#include <numbers>

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>

#include <cppdlr/utils.hpp>
#include <numbers>

using namespace cppdlr;
using namespace nda;

/**
* @brief Test Gaussian quadrature nodes and weights.
*/
TEST(adapinterp, gaussquad) {

  auto [xgl1, wgl1] = gaussquad(1);
  EXPECT_NEAR(xgl1(0), 0.0, 1.0e-15);

  EXPECT_NEAR(wgl1(0), 2.0, 1.0e-15);

  auto [xgl2, wgl2] = gaussquad(2);
  EXPECT_NEAR(xgl2(0), -1.0 / sqrt(3.0), 1.0e-15);
  EXPECT_NEAR(xgl2(1), 1.0 / sqrt(3.0), 1.0e-15);

  EXPECT_NEAR(wgl2(0), 1.0, 1.0e-15);
  EXPECT_NEAR(wgl2(1), 1.0, 1.0e-15);

  auto [xgl3, wgl3] = gaussquad(3);

  EXPECT_NEAR(xgl3(0), -sqrt(3.0/5.0), 1.0e-15);
  EXPECT_NEAR(xgl3(1), 0.0, 1.0e-15);
  EXPECT_NEAR(xgl3(2), sqrt(3.0/5.0), 1.0e-15);

  EXPECT_NEAR(wgl3(0), 5.0/9.0, 1.0e-15);
  EXPECT_NEAR(wgl3(1), 8.0/9.0, 1.0e-15);
  EXPECT_NEAR(wgl3(2), 5.0/9.0, 1.0e-15);

  auto [xgl4, wgl4] = gaussquad(4);

  EXPECT_NEAR(xgl4(0), -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)), 1.0e-15);
  EXPECT_NEAR(xgl4(1), -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)), 1.0e-15);
  EXPECT_NEAR(xgl4(2), sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)), 1.0e-15);
  EXPECT_NEAR(xgl4(3), sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)), 1.0e-15);

  EXPECT_NEAR(wgl4(0), (18.0 - sqrt(30.0))/36.0, 1.0e-15);
  EXPECT_NEAR(wgl4(1), (18.0 + sqrt(30.0))/36.0, 1.0e-15);
  EXPECT_NEAR(wgl4(2), (18.0 + sqrt(30.0))/36.0, 1.0e-15);
  EXPECT_NEAR(wgl4(3), (18.0 - sqrt(30.0))/36.0, 1.0e-15);

  auto [xgl5, wgl5] = gaussquad(5);

  EXPECT_NEAR(xgl5(0), -1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0)), 1.0e-15);
  EXPECT_NEAR(xgl5(1), -1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)), 1.0e-15);
  EXPECT_NEAR(xgl5(2), 0.0, 1.0e-15);
  EXPECT_NEAR(xgl5(3), 1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0)), 1.0e-15);
  EXPECT_NEAR(xgl5(4), 1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0)), 1.0e-15);

  EXPECT_NEAR(wgl5(0), (322.0 - 13.0 * sqrt(70.0))/900.0, 1.0e-15);
  EXPECT_NEAR(wgl5(1), (322.0 + 13.0 * sqrt(70.0))/900.0, 1.0e-15);
  EXPECT_NEAR(wgl5(2), 128.0/225.0, 1.0e-15);
  EXPECT_NEAR(wgl5(3), (322.0 + 13.0 * sqrt(70.0))/900.0, 1.0e-15);
  EXPECT_NEAR(wgl5(4), (322.0 - 13.0 * sqrt(70.0))/900.0, 1.0e-15);


}

/**
* @brief Test adaptive Gauss quadrature for semi-circular function.
*/
TEST(adapinterp, semicirc) {

  int n      = 16;      // Quadrature order
  double tol = 1.0e-14; // Absolute tolerance

  // Define integrand
  double a = 0.0, b = 2.0;
  auto f = [a, b](nda::array<double, 1> x) -> nda::array<double, 1> {
    double c = (a + b) / 2;
    double r = (b - a) / 2;
    return sqrt(r * r - pow(x - c, 2));
  };

  // Gauss-Legendre nodes and weights
  auto [xgl, wgl] = gaussquad(n);

  // Compute integral
  double intgrl = adapgl<double>(f, a, b, tol, xgl, wgl);

  // Compare with exact result
  EXPECT_NEAR(intgrl, std::numbers::pi / 2, tol);
}
