#include <gtest/gtest.h>
#include <cppdlr/utils.hpp>

using namespace cppdlr;

TEST(barycheb, interp) {

  int n = 16;

  barycheb mybarycheb(n);

  auto& xc = mybarycheb.getnodes();

  nda::vector<double> f(n);

  for (int i = 0; i < n; i++) { f(i) = std::cos(xc(i)); }

  double x = 0.378492;

  auto val1 = mybarycheb.interp(x,f);
  auto val2 = std::cos(x);

  EXPECT_NEAR(val1, val2, 1.0e-14);
}
