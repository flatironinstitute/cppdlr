#include <gtest/gtest.h>
#include <cppdlr/dlr_basis.hpp>

using namespace cppdlr;

TEST(dlr_basis, construct) {

    double lambda = 1000.0;
    double eps = 1e-14;
    auto om = dlr_freq(lambda,eps);

    int r = om.size();

    std::cout << r << std::endl;

  }