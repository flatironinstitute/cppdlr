#include "dlr_build.hpp"
#include "utils.hpp"
#include "dlr_kernels.hpp"

using namespace std;
using namespace nda;

// auto _ = range::all;

namespace cppdlr {

  fineparams::fineparams(double lambda, int p)
     :

       // All values have been chosen empirically to give fine discretization of
       // Lehmann kernel accurate to double machine precision

       lambda(lambda),
       p(p),
       npt(max((int)ceil(log(lambda) / log(2.0)) - 2, 1)),
       npo(max((int)ceil(log(lambda) / log(2.0)), 1)),
       nt(2 * p * npt),
       no(2 * p * npo) {}

  nda::vector<double> get_omfine(fineparams &fine) {

    int p   = fine.p;
    int npo = fine.npo;

    auto bc = barycheb(p);             // Get barycheb object for Chebyshev nodes
    auto xc = (bc.getnodes() + 1) / 2; // Cheb nodes on [0,1]

    // Real frequency grid points

    auto om = nda::vector<double>(fine.no);

    double a = 0, b = 0;

    // Points on (0,lambda)

    a = 0.0;
    for (int i = 0; i < npo; ++i) {
      b                                           = fine.lambda / pow(2.0, npo - i - 1);
      om(range((npo + i) * p, (npo + i + 1) * p)) = a + (b - a) * xc;
      a                                           = b;
    }

    // Points on (-lambda,0)

    om(range(0, npo * p)) = -om(range(2 * npo * p - 1, npo * p - 1, -1));

    return om;
  }

  nda::vector<double> get_tfine(fineparams &fine) {

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

  nda::matrix<double> get_kfine(fineparams &fine, nda::vector<double> &t, nda::vector<double> &om) {

    auto _ = range::all;

    int nt = fine.nt;
    int no = fine.no;

    nda::matrix<double> kmat(nt, no);

    for (int i = 0; i < nt / 2; ++i) {
      for (int j = 0; j < no; ++j) { kmat(i, j) = kfun(t(i), om(j)); }
    }

    kmat(range(nt / 2, nt), _) = kmat(range(nt / 2 - 1, -1, -1), range(no - 1, -1, -1));

    return kmat;
  }

  std::tuple<double, double> get_kfineerr(fineparams &fine, nda::vector<double> &t, nda::vector<double> &om, nda::matrix<double> kmat) {

    auto _ = range::all;

    int nt  = fine.nt;
    int no  = fine.no;
    int p   = fine.p;
    int npt = fine.npt;
    int npo = fine.npo;

    // Get fine composite Chebyshev time and frequency grids with double the
    // number of points per panel as the given fine grid

    fineparams fine2(fine.lambda, 2 * fine.p);
    auto ttst  = get_tfine(fine2);
    auto omtst = get_omfine(fine2);
    int p2     = fine2.p;

    // Interpolate values in K matrix to finer grids using barycentral Chebyshev
    // interpolation, and test against exact expression for kernel.

    barycheb bc(p);
    barycheb bc2(p2);
    auto xc = bc2.getnodes(); // Cheb nodes on [-1,1]

    double ktru, ktst, errtmp, errt, errom;

    // First test time discretization for each fixed frequency.

    errt = 0;
    for (int j = 0; j < no; ++j) {
      errtmp = 0;
      for (int i = 0; i < npt; ++i) { // Only need to test first half of matrix
        for (int k = 0; k < p2; ++k) {

          ktru = kfun(ttst(i * p2 + k), om(j));
          ktst = bc.interp(xc(k), kmat(range(i * p, (i + 1) * p), j));

          errtmp = max(errtmp, abs(ktru - ktst));
        }
      }
      errt = max(errt, errtmp / max_element(kmat(_, j)));
    }

    // Next test frequency discretization for each fixed time.

    errom = 0;
    for (int i = 0; i < nt / 2; ++i) {
      errtmp = 0;
      for (int j = 0; j < 2 * npo; ++j) {
        for (int k = 0; k < p2; ++k) {

          ktru = kfun(t(i), omtst(j * p2 + k));
          ktst = bc.interp(xc(k), kmat(i, range(j * p, (j + 1) * p)));

          errtmp = max(errtmp, abs(ktru - ktst));
        }
      }
      errom = max(errom, errtmp / max_element(kmat(i, _)));
    }

    return {errt, errom};
  }

} // namespace cppdlr
