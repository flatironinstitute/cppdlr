#include <gtest/gtest.h>
#include <cppdlr/dlr_build.hpp>
#include <cppdlr/dlr_kernels.hpp>

using namespace cppdlr;

TEST(dlr_build, fineparams) {

  fineparams fine1(100.0);

  // Check expected default values
  EXPECT_EQ(fine1.p, 24);
  EXPECT_EQ(fine1.npom, 7);
  EXPECT_EQ(fine1.npt, 5);
  EXPECT_EQ(fine1.nt , 2*fine1.p*fine1.npt);
  EXPECT_EQ(fine1.nom , 2*fine1.p*fine1.npom);

  // Check non-default value of p
  fineparams fine2(100.0,7);
  EXPECT_EQ(fine2.p, 7);

  // Check error triggered for non-permissible parameters.
  EXPECT_THROW({fineparams fine3(0.0);}, std::runtime_error);
  EXPECT_THROW({fineparams fine3(1.0, 0);}, std::runtime_error);
}

TEST(dlr_build, get_kfine) {

  fineparams fine(100.0);

  auto t = get_tfine(fine);
  auto om = get_omfine(fine);

  auto kmat = get_kfine(t,om);
  auto [errt,errom] = get_kfineerr(fine,t,om,kmat);

  EXPECT_LT(errt,   1e-14);
  EXPECT_LT(errom,  1e-14);

}
