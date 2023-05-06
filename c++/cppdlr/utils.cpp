#include <numbers>
#include "utils.hpp"
#include <nda/blas.hpp>

using namespace nda;
using std::numbers::pi;

static constexpr auto _ = range::all;

namespace cppdlr {

  barycheb::barycheb(int n) : xc(n), wc(n), n(n) {
    for (int i = 0; i < n; i++) {
      auto c        = (2 * i + 1) * pi / (2 * n);
      xc(n - i - 1) = std::cos(c);
      wc(n - i - 1) = (1 - 2 * (i % 2)) * std::sin(c);
    }
  }

  nda::vector<double> const &barycheb::getnodes() { return xc; }

  double barycheb::interp(double x, nda::vector_const_view<double> f) {

    for (int i = 0; i < n; ++i) {
      if (x == xc(i)) { return f(i); }
    }

    double num = 0, den = 0, dif = 0, q = 0;

    for (int i = 0; i < n; ++i) {
      dif = x - xc(i);
      q   = wc(i) / dif;
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
