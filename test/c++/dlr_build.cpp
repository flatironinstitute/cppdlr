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

  auto t = build_it_fine(fine);
  auto om = build_rf_fine(fine);

  auto kmat = build_k_it(t,om);
  auto [errt,errom] = geterr_k_it(fine,t,om,kmat);

  EXPECT_LT(errt,   1e-14);
  EXPECT_LT(errom,  1e-14);

}
