#include "dlr_build.hpp"
#include "utils.hpp"
#include "dlr_kernels.hpp"

using namespace std;
using namespace nda;

auto _ = range::all;

namespace cppdlr {

  fineparams::fineparams(double lambda) : 

    // All values have been chosen empirically to give fine discretization of
    // Lehmann kernel accurate to double machine precision

    lambda(lambda),
    p(24), 
    npt(max((int) ceil(log(lambda)/log(2.0))-2,1)),
    npo(max((int) ceil(log(lambda)/log(2.0)),1)),
    nt(2*p*npt),
    no(2*p*npo) { }
   

  tuple<nda::vector<double>, nda::vector<double>>
    get_finegrids(fineparams &fine) {

      int p = fine.p;
      int npt = fine.npt;
      int npo = fine.npo;

      barycheb bc(p); // Get barycheb object for Chebyshev nodes

      auto xc = (bc.getnodes()+1)/2; // Cheb nodes on [0,1]


      // Imaginary time grid points

      nda::vector<double> t(fine.nt);

      double a,b;

      // Points on (0,1/2)

      a = 0.0;
      for (int i = 0; i < npt; ++i) {
        b = 1.0/pow(2.0,npt-i);
        t(range(i*p,(i+1)*p)) = a + (b-a)*xc;
        a = b;
      }

      // Points on (1/2,1) in relative format

      t(range(npt*p,2*npt*p)) = -t(range(npt*p-1,-1,-1));


      // Real frequency grid points
      
      nda::vector<double> om(fine.no);

      // Points on (0,lambda)
      
      a = 0.0;
      for (int i = 0; i < npo; ++i) {
        b = fine.lambda/pow(2.0,npo-i-1);
        om(range((npo+i)*p,(npo+i+1)*p)) = a + (b-a)*xc;
        a = b;
      }

      // Points on (-lambda,0)

      om(range(0,npo*p)) = -om(range(2*npo*p-1,npo*p-1,-1));

      return {t,om};

    }

  nda::matrix<double> get_kfine(fineparams &fine, nda::vector<double> &t,
      nda::vector<double> &om) {

    int nt = fine.nt;
    int no = fine.no;
    
    nda::matrix<double> kmat(nt,no);

    for (int i=0; i<nt/2; ++i) {
      for (int j=0; j<no; ++j) {
        kmat(i,j) = kfun(t(i),om(j));
      }
    }

    kmat(range(nt/2,nt),_) = kmat(range(nt/2-1,-1,-1),range(no-1,-1,-1));

    return kmat;
    
    // TODO: Add error checking
    
  }

} // namespace cppdlr
