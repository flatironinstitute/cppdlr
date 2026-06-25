// Copyright (c) 2023 Simons Foundation
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
// Authors: Nils Wentzell, Jason Kaye

#include "dlr_imfreq.hpp"
#include "dlr_imtime.hpp"
#include "cppdlr/dlr_kernels.hpp"

using namespace nda;

namespace cppdlr {

  imfreq_ops::imfreq_ops(double lambda, nda::vector_const_view<double> dlr_rf, statistic_t statistic, bool symmetrize)
     : lambda_(lambda), statistic(statistic), r(dlr_rf.size()), dlr_rf(dlr_rf) {

    // Get # DLR imaginary frequency nodes; for symmetrized bosonic case, this
    // is DLR rank + 1, otherwise it is DLR rank
    niom   = (statistic == Boson && symmetrize) ? r + 1 : r;
    dlr_if = nda::vector<int>(niom);
    cf2if  = nda::matrix<dcomplex>(niom, r);

    // Get analytic continuation kernel at DLR frequencies, up to imaginary
    // frequency cutoff
    auto nmax = fineparams(lambda).nmax;
    auto kmat = build_k_if(nmax, dlr_rf, statistic);

    // Pivoted Gram-Schmidt to obtain DLR imaginary frequency nodes
    auto [q, norms, piv] = (symmetrize ? pivrgs_sym(kmat, niom) : pivrgs(kmat, 1e-100));
    std::sort(piv.begin(), piv.end()); // Sort pivots in ascending order

    // Sign-canonicalize for reproducible grid selection: the Matsubara node set is degenerate with
    // its mirror under nu_n -> -nu_n (row p -> max_piv - p). See flip_to_mirror.
    int max_piv = static_cast<int>(kmat.extent(0)) - 1; // mirror row of p under nu_n -> -nu_n
    if (flip_to_mirror(piv, max_piv)) {
      for (auto &p : piv) p = max_piv - p; // ascending -> descending,
      std::ranges::reverse(piv);           // so reverse restores ascending order
    }

    for (int i = 0; i < niom; ++i) { dlr_if(i) = piv(i) - nmax; }

    // Obtain coefficients to imaginary frequency values transformation matrix
    for (int i = 0; i < niom; ++i) {
      for (int j = 0; j < r; ++j) { cf2if(i, j) = kmat(piv(i), j); }
    }

    if (!(symmetrize && statistic == Boson)) {
      // Prepare imaginary time values to coefficients transformation by computing
      // LU factors of coefficient to imaginary time matrix
      if2cf.lu  = nda::matrix<dcomplex>(cf2if);
      if2cf.piv = nda::vector<int>(r);
      lapack::getrf(if2cf.lu, if2cf.piv);
    }
  }

  nda::vector<dcomplex> imfreq_ops::build_evalvec(double beta, int n) const {

    auto kvec = nda::vector<dcomplex>(r);
    for (int l = 0; l < r; ++l) { kvec(l) = beta * k_if(n, dlr_rf(l), statistic); }

    return kvec;
  }

  nda::vector<dcomplex> imfreq_ops::build_evalvec(int n) const {

    auto kvec = nda::vector<dcomplex>(r);
    for (int l = 0; l < r; ++l) { kvec(l) = k_if(n, dlr_rf(l), statistic); }

    return kvec;
  }

  nda::matrix<dcomplex> build_if2it(imfreq_ops const &ifops, imtime_ops const &itops) {

    int r    = itops.rank();
    int niom = ifops.get_ifnodes().size();

    if (niom == r) {
      // Square case: use LU solve
      auto if2it = nda::matrix<dcomplex>(itops.get_cf2it());

      auto lu  = nda::matrix<dcomplex>(ifops.get_if2cf_lu());
      auto piv = ifops.get_if2cf_piv();

      // Solve cf2if^T C^T = cf2it^T to obtain C = cf2it * cf2if^{-1}
      lapack::getrs(nda::transpose(lu), nda::transpose(if2it), piv);

      return if2it;
    } else {
      // Non-square case (niom > r): use least squares solve
      // Solve cf2if^T C^T = cf2it^T to obtain C = cf2it * cf2if^+
      auto A                            = nda::matrix<dcomplex, nda::F_layout>(nda::transpose(ifops.get_cf2if())); // r x niom
      auto B                            = nda::matrix<dcomplex, nda::F_layout>(niom, r);                           // niom x r
      B                                 = 0;
      B(nda::range(r), nda::range::all) = nda::transpose(itops.get_cf2it());

      auto s       = nda::vector<double>(r);
      double rcond = 0;
      int rank     = 0;
      lapack::gelss(A, B, s, rcond, rank);

      return nda::transpose(B);
    }
  }

} // namespace cppdlr
