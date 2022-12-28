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
    nda::array<double, 3> vals2coefs(nda::array_const_view<double, 3> g);

    /** Transform DLR coefficients of Green's function to its values on DLR imaginary time grid */
    nda::array<double, 3> coefs2vals(nda::array_const_view<double, 3> gc);

    /** Least squares fit of Green's function given by imaginary time samples, returning its DLR coefficients */
    nda::array<double, 3> data2coefs(nda::vector_const_view<double> t, nda::array_const_view<double, 3> g);

    /** Evaluate DLR expansion given by its DLR coefficients at an imaginary time point */
    nda::matrix<double> coefs2eval(nda::array_const_view<double, 3> gc, double t);

    /** Evaluate DLR expansion given by its values on DLR imaginary time grid at an imaginary time points */
    nda::matrix<double> vals2eval(nda::array_const_view<double, 3> gc);

    /** Access DLR imaginary time nodes*/
    nda::vector_const_view<double> get_itnodes() const;

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