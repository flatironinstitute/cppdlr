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
// Authors: Jason Kaye

#pragma once
#include "nda/layout_transforms.hpp"
#include <nda/nda.hpp>
#include <cppdlr/dlr_imtime.hpp>
#include <cppdlr/dlr_kernels.hpp>

#include <nda/linalg/eigenelements.hpp>
#include <type_traits>

namespace cppdlr {

  /**
  * @class dyson_it
  * @brief Class for solving Dyson equation in imaginary time
  * @tparam Ht Type of Hamiltonian
  */

  // Type of Hamiltonian, and scalar type of Hamiltonian
  template <typename Ht, nda::Scalar Sh = std::conditional_t<std::floating_point<Ht>, Ht, get_value_t<Ht>>>
    requires(std::floating_point<Ht> || nda::MemoryMatrix<Ht>)
  class dyson_it {

    public:
    /**
    * @brief Constructor for dyson_it
    * @param[in] beta Inverse temperature
    * @param[in] itops DLR imaginary time object
    * @param[in] h Hamiltonian
    * @param[in] mu Chemical potential
    *
    * \note Hamiltonian must either be a symmetric matrix, a Hermitian matrix,
    * or a real scalar.
    */
    dyson_it(double beta, imtime_ops itops, Ht const &h, double mu = 0) : beta(beta), itops_ptr(std::make_shared<imtime_ops>(itops)) {
      // dyson_it object contains a shared pointer to the imtime_ops object
      // itops. This is done to avoid making a copy of itops, which is meant to
      // handle all imaginary time operations on the given DLR imaginary time
      // grid.

      int r    = itops_ptr->rank();           // DLR rank
      auto g0  = free_gf(beta, itops, h, mu); // Free Green's function (right hand side of Dyson equation
      auto g0c = itops_ptr->vals2coefs(g0);   // DLR coefficients of free Green's function

      // Get matrix of convolution by free Green's function
      g0mat = itops_ptr->convmat(beta, Fermion, g0c);

      // Get right hand side of Dyson equation
      if constexpr (std::floating_point<Ht>) { // If h is real scalar, rhs is a vector
        norb = 1;
        rhs  = g0;
      } else { // Otherwise, rhs is given by g0 w/ some indices transposed (for compatibility w/ LAPACK)
        norb = h.shape(0);
        rhs.resize(norb, r, norb);
        rhs = permuted_indices_view<nda::encode<3>({1, 2, 0})>(g0);
      }
    }

    /**
    * @brief Solve Dyson equation for given self-energy
    *
    * @tparam Tsig Type of self-energy
    * @param[in] sig Self-energy at DLR imaginary time nodes
    *
    * @return Green's function at DLR imaginary time nodes
    *
    * \note Free Green's function (right hand side of Dyson equation) specified
    * at construction of dyson_it object
    */
    // Tsig is type of sig. Tg is type of Green's
    // function, which is of type Tsig, with scalar type replaced by the common
    // type of the Hamiltonian's scalar type, and the self-energy's scalar type
    // (real if both are real, complex otherwise).
    template <nda::MemoryArray Tsig, nda::MemoryArray Tg = make_common_t<Tsig, Sh, nda::get_value_t<Tsig>>> Tg solve(Tsig const &sig) const {

      int r     = itops_ptr->rank();          // DLR rank
      auto sigc = itops_ptr->vals2coefs(sig); // DLR coefficients of self-energy

      // Obtain Dyson equation system matrix I - G0 * Sig, where G0 and Sig are the
      // matrices of convolution by the free Green's function and self-energy,
      // respectively.
      auto sysmat = make_regular(nda::eye<double>(r * norb) - g0mat * itops_ptr->convmat(beta, Fermion, sigc));

      // Factorize system matrix
      auto ipiv = nda::vector<int>(r * norb);
      nda::lapack::getrf(sysmat, ipiv);

      // Solve Dyson equation
      auto g    = Tg(sig.shape());                                                    // Declare Green's function
      g         = rhs;                                                                // Get right hand side of Dyson equation
      auto g_rs = nda::matrix_view<get_value_t<Tg>>(nda::reshape(g, norb, r * norb)); // Reshape g to be compatible w/ LAPACK
      nda::lapack::getrs(sysmat, g_rs, ipiv);                                         // Back solve
      if constexpr (std::floating_point<Ht>) {                                        // If h is scalar, g is scalar-valued
        return g;
      } else { // Otherwise, g is matrix-valued, and need to transpose some indices after solve to undo LAPACK formatting
        return permuted_indices_view<nda::encode<3>({2, 0, 1})>(g);
      }
    }

    private:
    double beta;                           ///< Inverse temperature
    std::shared_ptr<imtime_ops> itops_ptr; ///< shared pointer to imtime_ops object
    int norb;                              ///< Number of orbital indices

    typename std::conditional_t<std::floating_point<Ht>, nda::array<Sh, 1>, nda::array<Sh, 3>>
       rhs; ///< Right hand side of Dyson equation (in format compatible w/ LAPACK); vector if Hamiltonian is scalar, rank-3 array otherwise
    nda::matrix<Sh> g0mat; ///< Matrix of convolution by free Green's function
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
  * @param[in] h Hamiltonian
  * @param[in] mu Chemical potential
  *
  * @return Green's function at DLR imaginary time nodes
  *
  * \note Hamiltonian must either be a symmetric matrix, a Hermitian matrix,
  * or a real scalar.
  */
  template <typename Ht>
    requires(std::floating_point<Ht> || nda::MemoryMatrix<Ht>)
  // If h is scalar, return scalar-valued Green's function; if h is matrix,
  // return matrix-valued Green's function
  auto free_gf(double beta, imtime_ops const &itops, Ht const &h, double mu = 0) {

    int r = itops.rank();

    if constexpr (std::floating_point<Ht>) { // If h is scalar, return scalar-valued Green's function
      auto g = nda::array<Ht, 1>(r);
      for (int i = 0; i < r; i++) { g(i) = k_it(itops.get_itnodes(i), beta * (h - mu)); }
      return g;
    } else { // Otherwise, return matrix-valued Green's function

      int norb = h.shape(0);

      // Diagonalize Hamiltonian
      auto [eval, evec] = nda::linalg::eigenelements(h);

      // Get free Green's function
      auto g = nda::array<nda::get_value_t<Ht>, 3>(r, norb, norb);
      g      = 0;
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < norb; j++) { g(i, j, j) = k_it(itops.get_itnodes(i), beta * (eval(j) - mu)); }
        g(i, _, _) = matmul(evec, matmul(g(i, _, _), transpose(conj(evec))));
      }

      return g;
    }
  }
} // namespace cppdlr