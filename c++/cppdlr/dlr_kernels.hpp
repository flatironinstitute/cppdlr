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
// Authors: jasonkaye

#pragma once
#include <complex>

namespace cppdlr {

  /** Evaluate Lehmann kernel in imaginary time with relative format for time variable */

  /**
  * @brief Evaluate analytic continuation kernel in imaginary time (relative time format)
  *
  * @param[in] t Imaginary time value (relative format)
  * @param[in] om Real frequency value
  *
  * @return Value K(t,om) of analytic continuation kernel
  *
  * \note We use dimensionless variables, i.e. beta = 1, so K(t,om) =
  * exp(-tau*om)/(1+exp(-om))
  */
  double k_it(double t, double om);

  /**
  * @brief Evaluate analytic continuation kernel in imaginary time (absolute time format)
  *
  * @param[in] t Imaginary time value (absolute format)
  * @param[in] om Real frequency value
  *
  * @return Value K(t,om) of analytic continuation kernel
  *
  * \note We use dimensionless variables, i.e. beta = 1, so K(t,om) =
  * exp(-tau*om)/(1+exp(-om))
  */
  double k_it_abs(double t, double om);

  /** Evaluate Lehmann kernel in imaginary frequency */

  /**
  * @brief Evaluate analytic continuation kernel in imaginary frequency
  *
  * @param[in] n Imaginary frequency index
  * @param[in] om Real frequency value
  *
  * @return Value K(n,om) of analytic continuation kernel
  *
  * \note We use dimensionless variables, i.e. beta = 1, and do not specify
  * fermionic or bosonic imaginary frequencies here, so K(n,om) = 1/(n*pi*i +
  * om)
  */
  std::complex<double> k_if(int n, double om);

} // namespace cppdlr
