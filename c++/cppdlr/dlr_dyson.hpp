#pragma once
#include <nda/nda.hpp>
#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_kernels.hpp>

#include <nda/linalg/eigenelements.hpp>

namespace cppdlr {

  /**
  * @class dyson_it
  * @brief Class for solving Dyson equation in imaginary time
  */

  template <typename T>
    requires(std::floating_point<T> || nda::MemoryMatrix<T>)
  class dyson_it {

    public:
    /**
    * @brief Constructor for dyson_it
    * @tparam T Scalar type of Hamiltonian
    * @param[in] beta Inverse temperature
    * @param[in] itops DLR imaginary time object
    * @param[in] mu Chemical potential
    * @param[in] h Hamiltonian
    *
    * \note Hamiltonian must either be a symmetric matrix, a Hermitian matrix,
    * or a real scalar.
    */
    dyson_it(double beta, imtime_ops &itops, double mu, T const &h) : beta(beta), itops(itops), mu(mu), h(h) {

      int r = itops.rank(); // DLR rank

      // Get free Green's function
      auto g0  = free_gf(beta, itops, mu, h);
      auto g0c = itops.vals2coefs(g0); // DLR coefficients

      // Get matrix of convolution by free Green's function
      g0mat = itops.convmat(beta, Fermion, g0c);

      // Get right hand side of Dyson equation
      if constexpr (std::floating_point<T>) { // If h is real scalar, rhs is a vector
        norb = 1;
        rhs  = g0;
      } else { // Otherwise, rhs is given by g0 w/ some indices transposed (for compatibility w/ LAPACK)
        norb = h.shape(0);
        rhs  = nda::array<get_value_t<T>, 3>(norb, r, norb);
        for (int i = 0; i < r; i++) {
          for (int j = 0; j < norb; j++) {
            for (int k = 0; k < norb; k++) { rhs(k, i, j) = g0(i, j, k); }
          }
        }
      }
    }

    /**
    * @brief Solve Dyson equation for given self-energy
    *
    * @tparam Tsig Scalar type of self-energy
    * @param[in] sig Self-energy at DLR imaginary time nodes
    *
    * @return Green's function at DLR imaginary time nodes
    *
    * \note Free Green's function (right hand side of Dyson equation) specified
    * at construction of dyson_it object
    */
    template <nda::MemoryArray Tsig> typename Tsig::regular_type solve(Tsig const &sig) const {

      int r     = itops.rank();          // DLR rank
      auto sigc = itops.vals2coefs(sig); // DLR coefficients of self-energy

      // Obtain Dyson equation system matrix I - G0 * Sig, where G0 and Sig are the
      // matrices of convolution by the free Green's function and self-energy,
      // respectively.
      auto sysmat = make_regular(nda::eye<double>(r * norb) - g0mat * itops.convmat(beta, Fermion, sigc));

      // Factorize system matrix
      auto ipiv = nda::vector<int>(r * norb);
      nda::lapack::getrf(sysmat, ipiv);

      // Solve Dyson equation
      auto g = Tsig(sig.shape());             // Could have complex self-energy w/ real G0; make g have type of sig
      g      = rhs;                           // Get right hand side of Dyson equation
      if constexpr (std::floating_point<T>) { // If h is scalar, g is scalar-valued
        nda::lapack::getrs(sysmat, g, ipiv);  // Back solve
        return g;
      } else { // Otherwise, g is matrix-valued, and need to transpose some indices after solve to undo LAPACK formatting
        auto g_rs = nda::matrix_view<get_value_t<Tsig>>(nda::reshaped_view(g, std::array<int, 2>({norb, r * norb})));
        nda::lapack::getrs(sysmat, g_rs, ipiv);
        auto g_return = nda::array<get_value_t<Tsig>, 3>(r, norb, norb);
        for (int i = 0; i < r; i++) {
          for (int j = 0; j < norb; j++) {
            for (int k = 0; k < norb; k++) { g_return(i, j, k) = g(k, i, j); }
          }
        }
        return g_return;
      }
    }

    private:
    double beta;       ///< Inverse temperature
    imtime_ops &itops; ///< imtime_ops object
    double mu;         ///< Chemical potential
    T const &h;        ///< Hamiltonian (must be symmetric or Hermitian matrix, or real scalar)

    int norb; ///< Number of orbital indices
    // TODO: can these conditionals be cleaned up? Why doesn't auto work
    // (shouldn't it be possible to deduce return type at compile time?)?
    typename std::conditional<std::floating_point<T>, nda::array<T, 1>, nda::array<get_value_t<T>, 3>>::type
       g0;                                  ///< Free Green's function at DLR nodes
    nda::matrix<nda::get_value_t<T>> g0mat; ///< Matrix of convolution by free Green's function
    typename std::conditional<std::floating_point<T>, nda::array<T, 1>, nda::array<get_value_t<T>, 3>>::type
       rhs; ///< Right hand side of Dyson equation
  };

  /**
  * @brief Compute free-particle imaginary time Green's function for a given
  * Hamiltonian
  *
  * The Green's function is computed by diagonalizing the Hamiltonian, and is
  * returned by its values at the DLR imaginary time nodes.
  *
  * @param[in] beta Inverse temperature
  * @param[in] it imtime_ops object
  * @param[in] mu Chemical potential
  * @param[in] h Hamiltonian
  *
  * @return Green's function at DLR imaginary time nodes
  *
  * \note Hamiltonian must either be a symmetric matrix, a Hermitian matrix,
  * or a real scalar.
  */
  template <typename T>
    requires(std::floating_point<T> || nda::MemoryMatrix<T>)
  // If h is scalar, return scalar-valued Green's function; if h is matrix,
  // return matrix-valued Green's function
  auto free_gf(double beta, imtime_ops const &itops, double mu, T const &h) {

    int r = itops.rank();

    if constexpr (std::floating_point<T>) { // If h is scalar, return scalar-valued Green's function
      auto g = nda::array<T, 1>(r);
      for (int i = 0; i < r; i++) { g(i) = -k_it(itops.get_itnodes(i), beta * (h - mu)); }
      return g;
    } else { // Otherwise, return matrix-valued Green's function

      int norb = h.shape(0);

      // Diagonalize Hamiltonian
      auto [eval, evec] = nda::linalg::eigenelements(h);

      // Get free Green's function
      auto g = nda::array<nda::get_value_t<T>, 3>(r, norb, norb);
      g      = 0;
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < norb; j++) { g(i, j, j) = -k_it(itops.get_itnodes(i), beta * (eval(j) - mu)); }
        g(i, _, _) = evec * matmul(g(i, _, _), transpose(conj(evec)));
      }

      return g;
    }
  }
} // namespace cppdlr