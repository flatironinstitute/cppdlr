#include <gtest/gtest.h>
#include <cppdlr/dlr_build.hpp>
#include <cppdlr/dlr_kernels.hpp>

using namespace cppdlr;

TEST(dlr_build, fineparams) {

  fineparams fine1(100.0,7);
  EXPECT_EQ(fine1.p, 7);

  fineparams fine(100.0);

  EXPECT_EQ(fine.p, 24);
  EXPECT_EQ(fine.npt, 5);
  EXPECT_EQ(fine.npo, 7);
  EXPECT_EQ(fine.nt , 2*fine.p*fine.npt);
  EXPECT_EQ(fine.no , 2*fine.p*fine.npo);
}

TEST(dlr_build, get_kfine) {

  fineparams fine(100.0);

  auto t = get_tfine(fine);
  auto om = get_omfine(fine);

  auto kmat = get_kfine(fine,t,om);
  auto [errt,errom] = get_kfineerr(fine,t,om,kmat);

  EXPECT_LT(errt,   1e-14);
  EXPECT_LT(errom,  1e-14);

}
