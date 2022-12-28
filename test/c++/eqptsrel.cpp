#include <gtest/gtest.h>
#include <cppdlr/utils.hpp>

using namespace cppdlr;

TEST(eqptsrel, eqptsrel_test) {

  int n = 4;
  auto t     = eqptsrel(n);
  auto ttrue = nda::vector<double>(n);
  ttrue(0)   = 0;
  ttrue(1)   = 1.0 / 3;
  ttrue(2)   = -1.0 / 3;
  ttrue(3)   = 1;

  EXPECT_EQ(t, ttrue);

  n        = 5;
  t        = eqptsrel(n);
  ttrue    = nda::vector<double>(n);
  ttrue(0) = 0;
  ttrue(1) = 1.0 / 4;
  ttrue(2) = 1.0 / 2;
  ttrue(3) = -1.0 / 4;
  ttrue(4) = 1;

  EXPECT_EQ(t, ttrue);
}
