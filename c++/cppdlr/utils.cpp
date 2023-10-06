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

#include <numbers>
#include "utils.hpp"
#include <nda/blas.hpp>

using namespace nda;
using std::numbers::pi;

namespace cppdlr {

  barycheb::barycheb(int n) : x(n), w(n) {
    for (int i = 0; i < n; i++) {
      auto c       = (2 * i + 1) * pi / (2 * n);
      x(n - i - 1) = std::cos(c);
      w(n - i - 1) = (1 - 2 * (i % 2)) * std::sin(c);
    }
  }

  nda::vector<double> const &barycheb::getnodes() { return x; }

  double barycheb::interp(double xeval, nda::vector_const_view<double> f) { return baryinterp(x, w, f, xeval); }

  baryleg::baryleg(int n) : x(n), w(n) {
    auto [xgl, wgl] = gaussquad(n); // Get Gauss-Legendre nodes and weights
    for (int j = 0; j < n; j++) {
      x(j) = xgl(j);
      w(j) = pow(-1, j) * sqrt((1 - pow(x(j), 2)) * wgl(j)); // Barycentric weights
    }
  }

  nda::vector<double> const &baryleg::getnodes() { return x; }

  double baryleg::interp(double xeval, nda::vector_const_view<double> f) { return baryinterp(x, w, f, xeval); }

  double baryinterp(nda::vector_const_view<double> x, nda::vector_const_view<double> w, nda::vector_const_view<double> f, double xeval) {

    int n = x.size();

    for (int i = 0; i < n; ++i) {
      if (xeval == x(i)) { return f(i); }
    }

    double num = 0, den = 0, dif = 0, q = 0;

    for (int i = 0; i < n; ++i) {
      dif = xeval - x(i);
      q   = w(i) / dif;
      num = num + q * f(i);
      den = den + q;
    }

    return num / den;
  }

  nda::vector<double> eqptsrel(int n) {

    auto t = nda::vector<double>(n);

    for (int i = 0; i < n - 1; ++i) {
      if (i <= (n - 1) / 2) {
        t(i) = i * 1.0 / (n - 1);
      } else {
        t(i) = -(n - 1 - i) * 1.0 / (n - 1);
      }
    }
    t(n - 1) = 1;

    return t;
  }

  nda::vector<double> rel2abs(nda::vector_const_view<double> t) {

    auto t_abs = nda::vector<double>(t.size());
    for (int i = 0; i < t.size(); ++i) {
      if (t(i) < 0) {
        t_abs(i) = t(i) + 1.0;
      } else {
        t_abs(i) = t(i);
      }
    }

    return t_abs;
  }

  double rel2abs(double t) {
    if (t < 0) {
      return t + 1.0;
    } else {
      return t;
    }
  }

  nda::vector<double> abs2rel(nda::vector_const_view<double> t_abs) {

    auto t = nda::vector<double>(t_abs.size());
    for (int i = 0; i < t_abs.size(); ++i) {
      if (t_abs(i) > 0.5 && t_abs(i) < 1.0) {
        t(i) = t_abs(i) - 1.0;
      } else {
        t(i) = t_abs(i);
      }
    }

    return t;
  }

  double abs2rel(double t_abs) {
    if (t_abs > 0.5 && t_abs < 1.0) {
      return t_abs - 1.0;
    } else {
      return t_abs;
    }
  }

  std::tuple<nda::vector<double>, nda::vector<double>> gaussquad(int n) {

    auto xgl = nda::vector<double>(n); // Gauss-Legendre nodes
    auto wgl = nda::vector<double>(n); // Gauss-Legendre weights
    double x = 0, dx = 0;
    int convcount = 0;

    // Get Gauss-Legendre nodes
    xgl(n / 2) = 0;                   // If odd number of nodes, middle node is 0
    for (int i = 0; i < n / 2; i++) { // Loop through nodes
      convcount = 0;
      x         = cos((2 * i + 1) * std::numbers::pi / (2 * n)); // Initial guess: Chebyshev node
      while (true) {                                             // Newton iteration
        auto [p, dp] = leg_eval(n, x);
        dx           = -p / dp;
        x += dx; // Newton step
        if (std::abs(dx) < 1e-14) { convcount++; }
        if (convcount == 3) { break; } // If convergence tol hit 3 times, stop
      }
      xgl(i)         = -x;
      xgl(n - i - 1) = x; // Symmetric nodes
    }

    // Get Gauss-Legendre weights from formula w_i = -2 / ((n+1)*P_n'(x_i)*P_{n+1}(x_i)) (Atkinson '89, pg. 276)
    for (int i = 0; i < n / 2 + 1; i++) {
      auto [junk1, dp] = leg_eval(n, xgl(i));
      auto [p, junk2]  = leg_eval(n + 1, xgl(i)); // This is a bit inefficient, but who cares...
      wgl(i)           = -2 / ((n + 1) * dp * p);
      wgl(n - i - 1)   = wgl(i);
    }

    return {xgl, wgl};
  }

  std::tuple<double, double> leg_eval(int n, double x) {

    if (n == 0) { return {1.0, 0.0}; }
    if (n == 1) { return {x, 1.0}; }

    // Three-term recurrence and formula for derivative
    double p0 = 0.0, p1 = 1.0, p2 = x;
    for (int i = 1; i < n; i++) {
      p0 = p1;
      p1 = p2;
      p2 = ((2 * i + 1) * x * p1 - i * p0) / (i + 1);
    }
    return {p2, n * (x * p2 - p1) / (x * x - 1)};
  }

} // namespace cppdlr
