// Copyright (c) 2022-2023 Simons Foundation
// Copyright (c) 2023 Hugo Strand
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
// Authors: Hugo U. R. Strand, Jason Kaye

#include "dlr_imtime.hpp"
#include "cppdlr/dlr_kernels.hpp"
#include "dlr_build.hpp"
#include "utils.hpp"

using namespace nda;

namespace cppdlr {

  imtime_ops::imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf, statistic_t statistic, bool symmetrize)
     : lambda_(lambda), r(dlr_rf.size()), dlr_rf(dlr_rf) {

    dlr_it    = nda::vector<double>(r);
    cf2it     = nda::matrix<double>(r, r);
    it2cf.lu  = nda::matrix<double>(r, r);
    it2cf.piv = nda::vector<int>(r);

    // Get discretization of analytic continuation kernel on fine grid in
    // imaginary time, at DLR frequencies
    auto fine = fineparams(lambda);
    auto t    = build_it_fine(fine);
    auto kmat = build_k_it(t, dlr_rf);

    if (!(symmetrize && statistic == Boson)) { // Treat symmetrized bosonic case separately

      // Pivoted Gram-Schmidt to obtain DLR imaginary time nodes
      auto [q, norms, piv] = (symmetrize ? pivrgs_sym(kmat, 1e-100) : pivrgs(kmat, 1e-100));
      std::sort(piv.begin(), piv.end()); // Sort pivots in ascending order
      for (int i = 0; i < r; ++i) { dlr_it(i) = t(piv(i)); }

      // Obtain coefficients to imaginary time values transformation matrix
      for (int i = 0; i < r; ++i) {
        for (int j = 0; j < r; ++j) { cf2it(i, j) = kmat(piv(i), j); }
      }

    } else { // Symmetrized bosonic case: enforce tau = 1/2 as DLR imaginary time node

      auto kvec12 = build_k_it(0.5, dlr_rf); // K at tau = 1/2: K(1/2, om)

      // Pivoted Gram-Schmidt on rows of K matrix, augmented by vector K(tau, om),
      // to obtain DLR imaginary time nodes
      auto [q, norms, piv] = pivrgs_sym(kmat, kvec12, 1e-100);

      std::sort(piv.begin(), piv.end());      // Sort pivots in ascending order
      for (int i = 0; i < (r - 1) / 2; ++i) { // Fill in time nodes on [0,1/2)
        dlr_it(i) = t(piv(i + 1) - 1);        // Shift by 1 to obtain pivots of original matrix, not augmented matrix
      }
      dlr_it((r - 1) / 2) = 0.5;                  // tau = 1/2 is always chosen
      for (int i = (r - 1) / 2 + 1; i < r; ++i) { // Fill in time nodes on (1/2,1]
        dlr_it(i) = t(piv(i) - 1);
      }

      // Obtain coefficients to imaginary time values transformation matrix
      for (int i = 0; i < (r - 1) / 2; ++i) { // Entries corresponding to time nodes on [0,1/2)
        for (int j = 0; j < r; ++j) { cf2it(i, j) = kmat(piv(i + 1) - 1, j); }
      }
      for (int j = 0; j < r; ++j) { cf2it((r - 1) / 2, j) = kvec12(j); } // Entries corresponding to tau = 1/2
      for (int i = (r - 1) / 2 + 1; i < r; ++i) {                        // Entries corresponding to time nodes on (1/2,1]
        for (int j = 0; j < r; ++j) { cf2it(i, j) = kmat(piv(i) - 1, j); }
      }
    }

    // Prepare imaginary time values to coefficients transformation by computing
    // LU factors of coefficient to imaginary time matrix
    it2cf.lu = cf2it;
    lapack::getrf(it2cf.lu, it2cf.piv);
  }

  imtime_ops::imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf) : imtime_ops(lambda, dlr_rf, Fermion, NONSYM) {}

  nda::vector<double> imtime_ops::build_evalvec(double t) const {

    // TODO: can be further optimized to reduce # exponential evals.
    auto kvec = nda::vector<double>(r);
    if (t >= 0) {
      for (int l = 0; l < r; ++l) { kvec(l) = k_it_abs(t, dlr_rf(l)); }
    } else {
      for (int l = 0; l < r; ++l) { kvec(l) = k_it_abs(-t, -dlr_rf(l)); }
    }

    return kvec;
  }

} // namespace cppdlr
