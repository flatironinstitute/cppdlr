/** 
* @file dyson_ed_it.cpp
*
* @brief Tests for dyson_it class
*/

#include <gtest/gtest.h>
#include <nda/nda.hpp>
#include <nda/gtest_tools.hpp>
#include <cppdlr/cppdlr.hpp>

#include <nda/linalg/eigenelements.hpp>

using namespace cppdlr;
using namespace nda;

/**
* @brief Compute free Green's function for random real 3x3 Hamiltonian by exact
* diagonalization, and compare its upper-left 2x2 block with solution by
* Dyson equation in imaginary time
*
* Given a 3x3 Hamiltonian h, the free Green's function can be obtained by exact
* diagonalization. Then, by eliminating the final row, it can be shown that the
* upper-left 2x2 block of the Green's function can be obtained by solving the
* Dyson equation with Hamiltonian given by h(1:2,1:2), and self-energy given by
* h(1:2,3) * G0 * h(3,1:2), where G0 is the free Green's function with energy
* h(3,3). We use this to validate and demonstrate a 2x2 Dyson solver in
* imaginary time.
*/

TEST(dyson_it, dyson_vs_ed_real) {

  // --- Problem setup --- //

  // Set problem parameters
  double beta = 100;  // Inverse temperature
  double mu   = 0;    // Chemical potential
  int n       = 3;    // Number of sites for original Hamiltonian
  int norb    = 2;    // Number of sites for reduced Hamiltonian
  int ntst    = 1000; // Number of imaginary time test points

  // Set DLR parameters
  double lambda = 100;
  double eps    = 1.0e-14;

  // Get random nxn Hamiltonian w/ eigenvalues in [-1 1]
  auto a         = nda::matrix<double>(nda::rand<double>(std::array<int, 2>({n, n}))); // Random matrix
  a              = (a + transpose(a)) / 2;                                             // Make symmetric
  auto [eval, u] = nda::linalg::eigenelements(a);                                      // Random orthogonal matrix
  eval           = -1 + 2 * nda::rand<double>(std::array<int, 1>({n}));                // Random eigenvalues in [-1,1]
  auto h         = matmul(matmul(u, nda::diag(eval)), transpose(conj(u)));             // Random symmetric matrix

  // --- Build DLR --- //

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);
  int r      = itops.rank();

  // --- Set up and solve Dyson equation --- //

  // Get self-energy
  auto sig = nda::array<dcomplex, 3>(r, norb, norb);
  auto g33 = free_gf(beta, itops, mu, real(h(n - 1, n - 1)));
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < norb; j++) {
      for (int k = 0; k < norb; k++) { sig(i, j, k) = h(j, n - 1) * g33(i) * conj(h(k, n - 1)); }
    }
  }
  auto sigc = itops.vals2coefs(sig);

  // Solve Dyson equation
  auto dys = dyson_it(beta, itops, mu, h(range(norb), range(norb))); // Dyson imaginary time solver object
  auto g   = dys.solve(sig);                                         // Solve Dyson equation w/ given self-energy
  auto gc  = itops.vals2coefs(g);                                    // Get DLR coefficients

  // --- Compare with exact solution --- //

  // Get upper-left 2x2 block of exact solution
  auto gex  = make_regular(free_gf(beta, itops, mu, h)(_, range(norb), range(norb)));
  auto gexc = itops.vals2coefs(gex); // Get DLR coefficients

  // Get equispaced test grid in relative format
  auto ittst = eqptsrel(ntst);

  // Compare on test grid
  double t  = 0;
  auto gtst = nda::array<dcomplex, 3>(ntst, norb, norb);
  auto gtru = nda::array<dcomplex, 3>(ntst, norb, norb);

  for (int i = 0; i < ntst; ++i) {
    t             = ittst(i);
    gtst(i, _, _) = itops.coefs2eval(gc, t);
    gtru(i, _, _) = itops.coefs2eval(gexc, t);
  }

  // Check error
  std::cout << "Max error: " << max_element(abs((gtst - gtru))) << std::endl;
  EXPECT_LT(max_element(abs((gtst - gtru))), 1.0e-13);
}

/**
* @brief Compute free Green's function for random complex 3x3 Hamiltonian by exact
* diagonalization, and compare its upper-left 2x2 block with solution by
* Dyson equation in imaginary time
*
* Given a 3x3 Hamiltonian h, the free Green's function can be obtained by exact
* diagonalization. Then, by eliminating the final row, it can be shown that the
* upper-left 2x2 block of the Green's function can be obtained by solving the
* Dyson equation with Hamiltonian given by h(1:2,1:2), and self-energy given by
* h(1:2,3) * G0 * h(3,1:2), where G0 is the free Green's function with energy
* h(3,3). We use this to validate and demonstrate a 2x2 Dyson solver in
* imaginary time.
*/

TEST(dyson_it, dyson_vs_ed_cmplx) {

  // --- Problem setup --- //

  // Set problem parameters
  double beta = 100; // Inverse temperature
  double mu   = 0;   // Chemical potential
  int n       = 3;   // Number of sites for original Hamiltonian
  int norb    = 2;   // Number of sites for reduced Hamiltonian
  int ntst    = 100; // Number of imaginary time test points

  // Set DLR parameters
  double lambda = 100;
  double eps    = 1.0e-14;

  // Get random nxn Hamiltonian w/ eigenvalues in [-1 1]
  auto a = nda::matrix<dcomplex>(nda::rand<double>(std::array<int, 2>({n, n})) + 1i * nda::rand<double>(std::array<int, 2>({n, n}))); // Random matrix
  a      = (a + transpose(conj(a))) / 2;                                   // Make Hermitian
  auto [eval, u] = nda::linalg::eigenelements(a);                          // Random unitary matrix
  eval           = -1 + 2 * nda::rand<double>(std::array<int, 1>({n}));    // Random eigenvalues in [-1,1]
  auto h         = matmul(matmul(u, nda::diag(eval)), transpose(conj(u))); // Random Hermitian matrix

  // --- Build DLR --- //

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);
  int r      = itops.rank();

  // --- Set up and solve Dyson equation --- //

  // Get self-energy
  auto sig = nda::array<dcomplex, 3>(r, norb, norb);
  auto g33 = free_gf(beta, itops, mu, real(h(n - 1, n - 1)));
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < norb; j++) {
      for (int k = 0; k < norb; k++) { sig(i, j, k) = h(j, n - 1) * g33(i) * conj(h(k, n - 1)); }
    }
  }
  auto sigc = itops.vals2coefs(sig);

  // Solve Dyson equation
  auto dys = dyson_it(beta, itops, mu, h(range(norb), range(norb))); // Dyson imaginary time object
  auto g   = dys.solve(sig);                                         // Solve Dyson equation w/ given self-energy
  auto gc  = itops.vals2coefs(g);                                    // Get DLR coefficients

  // --- Compare with exact solution --- //

  // Get upper-left 2x2 block of exact solution
  auto gex  = make_regular(free_gf(beta, itops, mu, h)(_, range(norb), range(norb)));
  auto gexc = itops.vals2coefs(gex); // Get DLR coefficients

  // Get equispaced test grid in relative format
  auto ittst = eqptsrel(ntst);

  // Compare on test grid
  double t  = 0;
  auto gtst = nda::array<dcomplex, 3>(ntst, norb, norb);
  auto gtru = nda::array<dcomplex, 3>(ntst, norb, norb);

  for (int i = 0; i < ntst; ++i) {
    t             = ittst(i);
    gtst(i, _, _) = itops.coefs2eval(gc, t);
    gtru(i, _, _) = itops.coefs2eval(gexc, t);
  }

  // Check error
  std::cout << "Max error: " << max_element(abs((gtst - gtru))) << std::endl;
  EXPECT_LT(max_element(abs((gtst - gtru))), 1.0e-13);
}

/** 
* @brief Bethe lattice Green's function (semi-circular density of states) in
* imaginary time
* 
* Bethe lattice function is computed from the spectral representation by
* adaptive quadrature
*
* @param[in] c Quarter-bandwidth
* @param[in] beta Inverse temperature
* @param[in] t Imaginary time points in relative format
*
* @return Value of Green's function at t points
*/
nda::vector<double> g_bethe(double c, double beta, nda::vector_const_view<double> t) {

  int n      = 16;      // Quadrature order
  double tol = 1.0e-14; // Absolute tolerance

  // Define integrand (using omega symmetry, integral is int_0^2c sqrt(4*c^2 - om^2) *
  //(exp(-t*om) + exp(-(beta-t)*om))/(1+exp(-beta*om)) dom. Further, using
  //t symmetry, in relative format we can replace t by abs(t) and assume t in
  //[0, 1/2], so the integrand is numerically stable as written.

  // TODO: double check beta scaling
  double a = 0, b = 2 * c;
  auto f = [c, beta](nda::array<double, 1> om, double tt) -> nda::array<double, 1> {
    return -2 / M_PI * sqrt(4 * c * c - om * om) * (exp(-tt * beta * om) + exp(-(1 - tt) * beta * om)) / (1 + exp(-beta * om));
  };

  // Gauss-Legendre nodes and weights
  auto [xgl, wgl] = gaussquad(n);

  // Compute integral
  double tt = 0;
  auto g    = nda::vector<double>(t.size());
  for (int i = 0; i < t.size(); ++i) {
    tt   = abs(t(i));
    g(i) = adapgl<double>([tt, f](nda::array<double, 1> om) { return f(om, tt); }, a, b, tol, xgl, wgl);
  }

  return g;
}

/**
* @brief Compute Bethe lattice Green's function by solving Dyson equation with
* known self-energy
*/

TEST(dyson_it, dyson_bethe) {

  // --- Problem setup --- //

  // Set problem parameters
  double beta = 100;     // Inverse temperature
  double mu   = 0;       // Chemical potential
  double c    = 1.0 / 2; // Quarter-bandwidth
  int ntst    = 1000;    // Number of imaginary time test points

  // Set DLR parameters
  double lambda = 100;
  double eps    = 1.0e-14;

  // --- Build DLR --- //

  // Get DLR frequencies
  auto dlr_rf = build_dlr_rf(lambda, eps);

  // Get DLR imaginary time object
  auto itops = imtime_ops(lambda, dlr_rf);

  // --- Get true Green's function --- //

  auto gbethe = g_bethe(c, beta, itops.get_itnodes());

  // --- Set up and solve Dyson equation --- //

  // Get self-energy
  // TODO: this "make_regular" isn't ideal; any way to fix?
  auto sig = make_regular(c * c * gbethe);

  // Solve Dyson equation
  double h = 0;                            // Zero Hamiltonian
  auto dys = dyson_it(beta, itops, mu, h); // Dyson imaginary time object
  auto g   = dys.solve(sig);               // Solve Dyson equation w/ given self-energy

  // --- Compare with exact solution --- //

  auto gbethec = itops.vals2coefs(gbethe); // DLR coefficients
  auto gc     = itops.vals2coefs(g);

  auto ittst = eqptsrel(ntst); // Equispaced test grid in relative format

  // Compare on test grid
  double t  = 0;
  auto gtst = nda::vector<double>(ntst);
  auto gtru = nda::vector<double>(ntst);

  for (int i = 0; i < ntst; ++i) {
    t       = ittst(i);
    gtst(i) = itops.coefs2eval(gc, t);
    gtru(i) = itops.coefs2eval(gbethec, t);
  }

  // Check error
  std::cout << "Max error: " << max_element(abs((gtst - gtru))) << std::endl;
  EXPECT_LT(max_element(abs((gtst - gtru))), 1.0e-13);
}