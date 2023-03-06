#pragma once
#include <nda/nda.hpp>

namespace cppdlr {

  /**
   * Class responsible for all DLR imaginary time operations, including building
   * imaginary time grid and transformations
   */

  class imtime_ops {

    public:
    imtime_ops(double lambda, nda::vector_const_view<double> dlr_rf);

    /** Transform values of Green's function on DLR imaginary time grid to its DLR coefficients */

    // This method takes g of any type fulfilling the constraints of the NDA
    // MemoryArray concept. In particular, it can be an NDA vector, matrix,
    // array, or view of any of these. The notation T::regular_type indicates
    // that if a view is passed in, a vector/matrix/array (and not a view of
    // these) will still be returned. The only constraint is that the final
    // index of g be of dimension r. So, for a matrix-valued Green's function
    // G_ij(tau) with 1<=i,j<=n, we pass in an nxnxr array.

    template <nda::MemoryArray T> 
    typename T::regular_type vals2coefs(T const &g) {

      if (r != g.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // First, we reshape g to a matrix with second dimension r. Then we call a
      // vals->coefs function which operates on matrices. Then we return the
      // result to its original shape. 

      auto g_reshaped = nda::reshaped_view(g,std::array<long,2> {r,g.size()/r});
      auto gc = vals2coefs_mat(g_reshaped);
      // [Q] Why reshaped_view instead of reshape? But this works, and reshape doesn't compile.
      return nda::reshaped_view(gc,g.shape());

    }

    private:
    nda::matrix<double> vals2coefs_mat(nda::matrix_const_view<double> g);

    public:

    /** Transform DLR coefficients of Green's function to its values on DLR imaginary time grid */

    template <nda::MemoryArray T> 
    typename T::regular_type coefs2vals(T const &gc) {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // First, we reshape gc to a matrix with second dimension r. Then we call a
      // coefs->vals function which operates on matrices. Then we return the
      // result to its original shape. 

      auto gc_reshaped = nda::reshaped_view(gc,std::array<long,2> {r, gc.size()/r});
      auto g = coefs2vals_mat(gc_reshaped);
      // [Q] Why reshaped_view instead of reshape? But this works, and reshape doesn't compile.
      return nda::reshaped_view(g,gc.shape());

    }

    private:
    nda::matrix<double> coefs2vals_mat(nda::matrix_const_view<double> gc);

    public:

    /** Least squares fit of Green's function given by imaginary time samples, returning its DLR coefficients */
    nda::array<double, 3> data2coefs(nda::vector_const_view<double> t, nda::array_const_view<double, 3> g);

    /** Evaluate DLR expansion given by its DLR coefficients at an imaginary time point */
    
    template <nda::MemoryArray T>
    auto coefs2eval(T const &gc, double t) {

      if (r != gc.shape(0)) throw std::runtime_error("First dim of g != DLR rank r.");

      // Scalar-valued Green's functions are handled differently than matrix-valued Green's functions 
      if constexpr(T::rank==1) {
       return coefs2eval_vec(gc,t);
      } else {

      // First, we reshape gc to a matrix with second dimension r. Then we call a
      // coefs->value function which operates on matrices and returns a vector. Then we return the
      // result to its original shape (leaving out the final dimension, which is summed out). 

      auto gc_reshaped = nda::reshaped_view(gc,std::array<long,2> {r,gc.size()/r});

      std::array<long,T::rank - 1> shape_out;
      for (int i = 0; i < T::rank - 1; ++i) { shape_out[i] = gc.shape(i+1); }

      return reshape(coefs2eval_mat(gc_reshaped,t),shape_out);

      }

    }

    private:
    double coefs2eval_vec(nda::vector_const_view<double> gc, double t);
    nda::vector<double> coefs2eval_mat(nda::matrix_const_view<double> gc, double t);

    public:

    /** Get vector of evaluation of DLR expansion at a point **/
    nda::vector<double> get_kevalvec(double t);

    /** Evaluate DLR expansion given by its values on DLR imaginary time grid at an imaginary time point */
    nda::matrix<double> vals2eval(nda::array_const_view<double, 3> g);

    /** Access DLR imaginary time nodes*/
    nda::vector_const_view<double> get_itnodes() const {return dlr_it; };

    /** DLR rank **/
    int rank() const { return r; }

    private:
    int r;                      /// DLR rank
    nda::vector<double> dlr_rf; /// DLR frequencies
    nda::vector<double> dlr_it; /// DLR imaginary time nodes
    nda::matrix<double> cf2it;  /// Transformation matrix from DLR coefficients to values at DLR imaginary time nodes

    struct {
      nda::matrix<double> lu;
      nda::vector<int> piv;
    } it2cf; /// Transformation (represented by LU factors and pivots) from values at DLR imaginary time nodes to DLR coefficients
  };

} // namespace cppdlr
