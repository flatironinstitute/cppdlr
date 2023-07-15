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

#pragma once
#include <complex>

namespace cppdlr {

  /**
   * The Particle Statistic: Boson or Fermion
   */

  enum statistic_t { Boson = 0, Fermion = 1 };

  /**
  * @brief Evaluate analytic continuation kernel in imaginary time (relative time format)
  *
  * @param[in] t Imaginary time value (relative format)
  * @param[in] om Real frequency value
  *
  * @return Value K(t,om) of analytic continuation kernel
  *
  * \note We use dimensionless variables, i.e. beta = 1, so k_it(t,om) = -exp(-t
  * * om)/(1 + exp(-om)). Thus, to obtain the analytic continuation kernel
  * K(t,om) = -exp(-t*om)/(1+exp(-beta * om)) for a given inverse temperature
  * beta, multiply the frequency argument of k_it by beta. In other words, 
  * K(t,om) = k_it(t, beta * om).
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
  * \note We use dimensionless variables, i.e. beta = 1, so k_it(t,om) = -exp(-t
  * * om)/(1 + exp(-om)). Thus, to obtain the analytic continuation kernel
  * K(t,om) = -exp(-t*om)/(1+exp(-beta * om)) for a given inverse temperature
  * beta, multiply the frequency argument of k_it by beta. In other words, 
  * K(t,om) = k_it(t, beta * om).
  */
  double k_it_abs(double t, double om);

  /**
  * @brief Evaluate analytic continuation kernel in imaginary frequency
  *
  * @param[in] n Imaginary frequency index
  * @param[in] om Real frequency value
  * @param[in] statistic Particle Statistic: Boson or Fermion
  *
  * @return Value K(n,om) of analytic continuation kernel
  *
  * \note We use dimensionless variables, i.e. beta = 1, so k_if(n,om,Fermion) =
  * 1/((2*n+1)*pi*i - om), and k_if(n,om,Boson) = 1/(2*n*pi*i - om). Thus, to
  * obtain the analytic continuation kernel K(i nu_n, om) = 1 / (i nu_n - om)
  * with i nu_n = (2n+1) * i * pi / beta (fermionic case) or i nu_n = 2n * i *
  * pi / beta (bosonic case) for a given inverse temperature beta, multiply the
  * frequency argument of k_if by beta, and multiply the result by beta. In
  * other words, K(i nu_n, om) = beta * k_if(n, beta * om, statistic). One can
  * see this is correct by multiplying the numerator and denominator of K(i
  * nu_n, om) by beta.
  */
  std::complex<double> k_if(int n, double om, statistic_t statistic);

} // namespace cppdlr
