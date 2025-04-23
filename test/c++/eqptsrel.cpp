// Copyright (c) 2022 Simons Foundation
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
#include <cppdlr/utils.hpp>

using namespace cppdlr;

TEST(eqptsrel, eqptsrel_test) {

  int n      = 4;
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
