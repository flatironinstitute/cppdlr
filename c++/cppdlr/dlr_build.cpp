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
// Authors: Nils Wentzell, Jason Kaye

#include "dlr_build.hpp"
#include "utils.hpp"
#include "dlr_kernels.hpp"
#include <numbers>

using namespace std;
using namespace nda;

using std::numbers::pi;

namespace cppdlr {

  fineparams::fineparams(double lambda, int p)
     : // TODO: make alt constructor in which iommax is a parameter

       lambda(lambda),
       p(p),
       nmax(max((int)ceil(lambda), 20)), // TODO: is this a good choice?
       npom(static_cast<int>(max(ceil(log(lambda) / log(2.0)), 1.0))),
       npt(static_cast<int>(max(ceil(log(lambda) / log(2.0)) - 2, 1.0))),
       nom(2 * p * npom),
       nt(2 * p * npt) {

    if (lambda <= 0) throw std::runtime_error("Choose lambda > 0.");
    if (p <= 0) throw std::runtime_error("Choose p > 0.");
  }

  nda::vector<double> build_rf_fine(fineparams &fine) {

    int p    = fine.p;
    int npom = fine.npom;

    auto bc = barycheb(p);             // Get barycheb object for Chebyshev nodes
    auto xc = (bc.getnodes() + 1) / 2; // Cheb nodes on [0,1]

    // Real frequency grid points

    auto om = nda::vector<double>(fine.nom);

    double a = 0, b = 0;

    // Points on (0,lambda)

    a = 0.0;
    for (int i = 0; i < npom; ++i) {
      b                                             = fine.lambda / pow(2.0, npom - i - 1);
      om(range((npom + i) * p, (npom + i + 1) * p)) = a + (b - a) * xc;
      a                                             = b;
    }

    // Points on (-lambda,0)

    om(range(0, npom * p)) = -om(range(2 * npom * p - 1, npom * p - 1, -1));

    return om;
  }

  nda::vector<double> build_it_fine(fineparams &fine) {

    int p   = fine.p;
    int npt = fine.npt;

    auto bc = barycheb(p);             // Get barycheb object for Chebyshev nodes
    auto xc = (bc.getnodes() + 1) / 2; // Cheb nodes on [0,1]

    // Imaginary time grid points

    auto t = nda::vector<double>(fine.nt);

    double a = 0, b = 0;

    // Points on (0,1/2)

    for (int i = 0; i < npt; ++i) {
      b                            = 1.0 / pow(2.0, npt - i);
      t(range(i * p, (i + 1) * p)) = a + (b - a) * xc;
      a                            = b;
    }

    // Points on (1/2,1) in relative format

    t(range(npt * p, 2 * npt * p)) = -t(range(npt * p - 1, -1, -1));

    return t;
  }

  std::tuple<nda::vector<double>, nda::vector<double>> build_it_fine_wgt(fineparams &fine) {

    int p   = fine.p;
    int npt = fine.npt;

    auto [xgl, wgl] = gaussquad(p);  // Gauss-Legendre nodes and weights on [-1,1]
    xgl             = (xgl + 1) / 2; // Transform to [0,1]
    //auto bc = barycheb(p);             // Get barycheb object for Chebyshev nodes
    //auto xc = (bc.getnodes() + 1) / 2; // Cheb nodes on [0,1]

    // Imaginary time grid points

    auto t = nda::vector<double>(fine.nt);
    auto w = nda::vector<double>(fine.nt);

    double a = 0, b = 0;

    // Nodes and weights on (0,1/2)

    for (int i = 0; i < npt; ++i) {
      b                            = 1.0 / pow(2.0, npt - i);
      t(range(i * p, (i + 1) * p)) = a + (b - a) * xgl;
      w(range(i * p, (i + 1) * p)) = sqrt(((b - a) / 2) * wgl);
      a                            = b;
    }

    // Points on (1/2,1) in relative format

    t(range(npt * p, 2 * npt * p)) = -t(range(npt * p - 1, -1, -1));
    w(range(npt * p, 2 * npt * p)) = w(range(npt * p - 1, -1, -1));

    return {t, w};
  }

  nda::matrix<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> om) {

    int nt  = t.size();
    int nom = om.size();

    auto kmat = nda::matrix<double>(nt, nom);

    for (int i = 0; i < nt; ++i) {
      for (int j = 0; j < nom; ++j) { kmat(i, j) = k_it(t(i), om(j)); }
    }

    return kmat;
  }

  nda::matrix<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> w, nda::vector_const_view<double> om) {

    int nt  = t.size();
    int nom = om.size();

    auto kmat = nda::matrix<double>(nt, nom);

    for (int i = 0; i < nt; ++i) {
      for (int j = 0; j < nom; ++j) { kmat(i, j) = w(i) * k_it(t(i), om(j)); }
    }

    return kmat;
  }

  nda::vector<double> build_k_it(double t, nda::vector_const_view<double> om) {

    int nom = om.size();

    auto kvec = nda::vector<double>(nom);
    for (int j = 0; j < nom; ++j) { kvec(j) = k_it(t, om(j)); }

    return kvec;
  }

  nda::vector<double> build_k_it(nda::vector_const_view<double> t, double om) {

    int nt = t.size();

    auto kvec = nda::vector<double>(nt);
    for (int i = 0; i < nt; ++i) { kvec(i) = k_it(t(i), om); }

    return kvec;
  }

  nda::vector<double> build_k_it(nda::vector_const_view<double> t, nda::vector_const_view<double> w, double om) {

    int nt = t.size();

    auto kvec = nda::vector<double>(nt);
    for (int i = 0; i < nt; ++i) { kvec(i) = w(i) * k_it(t(i), om); }

    return kvec;
  }

  std::tuple<double, double> geterr_k_it(fineparams &fine, nda::vector_const_view<double> t, nda::vector_const_view<double> om,
                                         nda::matrix_const_view<double> kmat) {

    auto _ = range::all;

    int nt   = fine.nt;
    int nom  = fine.nom;
    int p    = fine.p;
    int npt  = fine.npt;
    int npom = fine.npom;

    // Get fine composite Chebyshev time and frequency grids with double the
    // number of points per panel as the given fine grid

    fineparams fine2(fine.lambda, 2 * fine.p);
    auto ttst  = build_it_fine(fine2);
    auto omtst = build_rf_fine(fine2);
    int p2     = fine2.p;

    // Interpolate values in K matrix to finer grids using barycentral Chebyshev
    // interpolation, and test against exact expression for kernel.

    barycheb bc(p);
    barycheb bc2(p2);
    auto xc = bc2.getnodes(); // Cheb nodes on [-1,1]

    double ktru = 0, ktst = 0, errtmp = 0;

    // First test time discretization for each fixed frequency.

    double errt = 0;
    for (int j = 0; j < nom; ++j) {
      errtmp = 0;
      for (int i = 0; i < npt; ++i) { // Only need to test first half of matrix
        for (int k = 0; k < p2; ++k) {

          ktru = k_it(ttst(i * p2 + k), om(j));
          ktst = bc.interp(xc(k), kmat(range(i * p, (i + 1) * p), j));

          errtmp = max(errtmp, abs(ktru - ktst));
        }
      }
      errt = max(errt, errtmp / max_element(kmat(_, j)));
    }

    // Next test frequency discretization for each fixed time.

    double errom = 0;
    for (int i = 0; i < nt / 2; ++i) {
      errtmp = 0;
      for (int j = 0; j < 2 * npom; ++j) {
        for (int k = 0; k < p2; ++k) {

          ktru = k_it(t(i), omtst(j * p2 + k));
          ktst = bc.interp(xc(k), kmat(i, range(j * p, (j + 1) * p)));

          errtmp = max(errtmp, abs(ktru - ktst));
        }
      }
      errom = max(errom, errtmp / max_element(kmat(i, _)));
    }

    return {errt, errom};
  }

  nda::matrix<dcomplex> build_k_if(int nmax, nda::vector_const_view<double> om, statistic_t statistic) {

    int nom = om.size();

    if (statistic == Fermion) {

      // 2*n+1 should go from -2*nmax+1 to 2*nmax-1, so n goes from -nmax to nmax-1
      auto kmat = nda::matrix<dcomplex>(2 * nmax, nom);

      for (int n = -nmax; n < nmax; ++n) {
        for (int j = 0; j < nom; ++j) { kmat(nmax + n, j) = k_if(n, om(j), statistic); }
      }
      return kmat;
    } else {

      // 2*n+1 should go from -2*nmax to 2*nmax, so n goes from -nmax to nmax
      auto kmat = nda::matrix<dcomplex>(2 * nmax + 1, nom);

      for (int n = -nmax; n <= nmax; ++n) {
        for (int j = 0; j < nom; ++j) { kmat(nmax + n, j) = k_if(n, om(j), statistic); }
      }
      return kmat;
    }
  }

  nda::vector<double> build_dlr_rf(double lambda, double eps, statistic_t statistic, bool symmetrize) {

    // Get fine grid parameters

    auto fine = fineparams(lambda);

    // Get fine grids in frequency and imaginary time

    //auto t  = build_it_fine(fine);
    auto [t, w] = build_it_fine_wgt(fine);

    auto om = build_rf_fine(fine);

    // Get discretization of analytic continuation kernel on fine grid (the K matrix)

    //auto kmat = build_k_it(t, om);
    auto kmat = build_k_it(t, w, om);

    if (!(symmetrize && statistic == Boson)) { // Treat symmetrized bosonic case separately

      // Pivoted Gram-Schmidt on columns of K matrix to obtain DLR frequencies

      auto [q, norms, piv] = (symmetrize ? pivrgs_sym(transpose(kmat), eps) : pivrgs(transpose(kmat), eps));
      long r               = norms.size();
      std::sort(piv.begin(), piv.end()); // Sort pivots in ascending order

      auto omega = nda::vector<double>(r);
      for (int i = 0; i < r; ++i) { omega(i) = om(piv(i)); }

      return omega;

    } else { // Symmetrized bosonic case: enforce omega = 0 chosen as DLR frequency

      //auto kvec0 = build_k_it(t, 0.0); // K at zero frequency: K(tau, 0)
      auto kvec0 = build_k_it(t, w, 0.0); // K at zero frequency: K(tau, 0)

      // Pivoted Gram-Schmidt on columns of K matrix, augmented by vector K(tau, 0), to obtain DLR frequencies
      auto [q, norms, piv] = pivrgs_sym(transpose(kmat), kvec0, eps);
      long r               = norms.size();
      std::sort(piv.begin(), piv.end()); // Sort pivots in ascending order

      auto omega = nda::vector<double>(r);
      for (int i = 0; i < (r - 1) / 2; ++i) { // Fill in negative frequencies
        omega(i) = om(piv(i + 1) - 1);        // Shift by 1 to obtain pivots of original matrix, not augmented matrix
      }
      omega((r - 1) / 2) = 0.0;                   // Zero frequency is always chosen
      for (int i = (r - 1) / 2 + 1; i < r; ++i) { // Fill in positive frequencies
        omega(i) = om(piv(i) - 1);
      }

      return omega;
    }
  }

  nda::vector<double> build_dlr_rf(double lambda, double eps) { return build_dlr_rf(lambda, eps, statistic_t::Fermion, NONSYM); }

} // namespace cppdlr
