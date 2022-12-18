#include <gtest/gtest.h>
#include <cppdlr/dlr_basis.hpp>

using namespace cppdlr;

TEST(dlr_basis, construct) {

    double lambda = 1000.0;
    double eps = 1e-14;
    auto dlrb = dlr_basis(lambda,eps);
    //dlr_basis dlrb(lambda,eps);

    auto om = dlrb.get_rf();
    int r = om.size();

    std::cout << r << std::endl;

  }