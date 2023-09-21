#include <nda/nda.hpp>
#include <cppdlr/cppdlr.hpp>

using namespace cppdlr;
using namespace nda;

/**
* @brief Green's function which is a random sum of poles
*
* G_ij(iom_n) = sum_l c_ijl K(n,om_ijl) with random c_ijl, om_ijl
*
* @param[in] norb      Number of orbital indices
* @param[in] beta      Inverse temperature
* @param[in] n         Imaginary frequency evaluation point index
* @param[in] statistic Particle Statistic: Fermion or Boson
*
* @return Green's function evaluated at iom_n
*/
nda::matrix<dcomplex> gfun(int norb, double beta, int n, statistic_t statistic) {

  int npeak = 5;

  auto g    = nda::matrix<dcomplex>(norb, norb);
  g         = 0;
  auto c    = nda::vector<double>(npeak);
  double om = 0;
  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {

      // Get random weights that sum to 1
      for (int l = 0; l < npeak; ++l) {
        c(l) = (sin(1000.0 * (i + 2 * j + 3 * l + 7)) + 1) / 2; // Quick and dirty rand # gen on [0,1]
      }
      c = c / sum(c);
      c = beta * c;

      // Evaluate Green's function
      for (int l = 0; l < npeak; ++l) {
        om = sin(2000.0 * (3 * i + 2 * j + l + 6)); // Rand # on [-1,1]
        g(i, j) += c(l) * k_if(n, beta * om, statistic);
      }
    }
  }

  return g;
}


/**
* @brief Compare symmetrized and unsymmetrized DLR for interpolation and
* evaluation of matrix-valued Green's function in imaginary time
*/
int main() {

  double lambda = 1000;  // DLR cutoff
  double eps    = 1e-10; // DLR tolerance
  auto statistic = Fermion; // Fermionic Green's function

  double beta = 1000;  // Inverse temperature
  int nmaxtst = 10000; // # imag time test points

  int norb = 2; // Orbital dimensions

  // Get DLR frequencies
  auto dlr_rf     = build_dlr_rf(lambda, eps);
  auto dlr_rf_sym = build_dlr_rf(lambda, eps, Fermion, SYM);

  int r    = dlr_rf.size();
  int rsym = dlr_rf_sym.size();

  // Get DLR imaginary frequency object
  auto ifops = imfreq_ops(lambda, dlr_rf, statistic);
  auto ifops_sym = imfreq_ops(lambda, dlr_rf_sym, statistic, SYM);

  // Sample Green's function G at DLR imaginary frequency nodes
  auto const &dlr_if = ifops.get_ifnodes();
  auto const &dlr_if_sym = ifops_sym.get_ifnodes();


  auto g = nda::array<dcomplex, 3>(r, norb, norb);
  auto g_sym = nda::array<dcomplex, 3>(rsym, norb, norb);
  for (int i = 0; i < r; ++i) { g(i, _, _) = gfun(norb, beta, dlr_if(i), statistic); }
  for (int i = 0; i < rsym; ++i) { g_sym(i, _, _) = gfun(norb, beta, dlr_if_sym(i), statistic); }

  // DLR coefficients of G
  auto gc     = ifops.vals2coefs(beta, g);
  auto gc_sym = ifops_sym.vals2coefs(beta, g_sym);

  // Compute L infinity error
  auto gtru  = nda::matrix<dcomplex>(norb, norb);
  auto gtst  = nda::matrix<dcomplex>(norb, norb);
  auto gtst_sym  = nda::matrix<dcomplex>(norb, norb);
  double err = 0, err_sym = 0;
  for (int n = -nmaxtst; n < nmaxtst; ++n) {
    gtru = gfun(norb, beta, n, statistic);
    
    gtst = ifops.coefs2eval(beta, gc, n);
    gtst_sym = ifops_sym.coefs2eval(beta, gc_sym, n);
    
    err  = std::max(err, max_element(abs(gtru - gtst)));
    err_sym  = std::max(err_sym, max_element(abs(gtru - gtst_sym)));
  }

  // Print results
  std::cout << "Unsymmetrized DLR rank = " << r << std::endl;
  std::cout << "Symmetrized DLR rank = " << rsym << std::endl;
  std::cout << "L infinity error for unsymmetrized DLR = " << err << std::endl;
  std::cout << "L infinity error for symmetrized DLR = " << err_sym << std::endl;
}