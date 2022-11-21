#include <numbers>
#include "utils.hpp"

using namespace nda;
using std::numbers::pi;

namespace cppdlr {

  barycheb::barycheb(int n) : xc(n), wc(n) {

    // vector<int> ns(n);
    // for (int i = 0; i < n; i++) { ns(i) = i + 1; }

    // auto c = (2 * ns - 1) * pi / (2 * n);
    // xc     = std::reverse(cos(c));
    // wc     = std::reverse(1-2*((ns-1)%2))*sin(c);

    for (int i = 0; i < n; i++) {
      auto c        = (2 * i + 1) * pi / (2 * n);
      xc(n - i - 1) = std::cos(c);
      wc(n - i - 1) = (1 - 2 * (i % 2)) * std::sin(c);
    }
  }

  nda::vector_view<double> barycheb::getnodes() { return xc; }

  double barycheb::interp(double &x, nda::vector<double> &f) {

    int n = xc.size(); // [Q] Do it this way, or make n a private variable of the class?

    for (int i = 0; i < n; i++) {
      if (x == xc(i)) { return f(i); }
    }

    double num = 0.0, den = 0.0;
    double dif, q;

    for (int i = 0; i < n; i++) {
      dif = x - xc(i);
      q   = wc(i) / dif;
      num = num + q * f(i);
      den = den + q;
    }

    return num / den;
  }

} // namespace cppdlr