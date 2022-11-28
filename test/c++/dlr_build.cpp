#include <gtest/gtest.h>
#include <cppdlr/dlr_build.hpp>
#include <cppdlr/dlr_kernels.hpp>

using namespace cppdlr;

TEST(dlr_build, fineparams) {

  fineparams fine(100.0);

  EXPECT_EQ(fine.p, 24);
  EXPECT_EQ(fine.npt, 5);
  EXPECT_EQ(fine.npo, 7);
  EXPECT_EQ(fine.nt , 2*fine.p*fine.npt);
  EXPECT_EQ(fine.no , 2*fine.p*fine.npo);
}

TEST(dlr_build, get_finegrids) {

  fineparams fine(10.0);

  auto [t,om] = get_finegrids(fine);

  std::cout << t << std::endl;
  std::cout << om << std::endl;

  auto kmat = get_kfine(fine,t,om);

}
