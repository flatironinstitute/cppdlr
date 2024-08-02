(changelog)=

# Changelog

## Version 1.2.0

This update contains minor additional functionality, bug fixes, and optimizations.

We thank all contributors: Jason Kaye, Hugo U. R. Strand, Nils Wentzell

### New features
* Added imtime_ops::get_ipmat function to expose matrix of imaginary time inner product
* Implemented transpose of values -> coefficients operation in imtime_ops::vals2coefs 

### Optimizations
* imtime_ops::inner_prod function can return double, rather than always complex double

### Other changes
* Test tolerance adjustments
* Minor documentation updates


## Version 1.1.1

cppdlr version 1.1.1 is a patch-release that includes
cmake build system improvements, continuous integration improvements,
and a backward compatibilty fix.

We thank all contributors: Hugo U.R. Strand, Nils Wentzell

Find below an itemized list of changes in this release.

* Synchronize ghactions workflow with main
* Use nda 1.3.x branch by default
* Fix issue with Doxygen output directory
* simplified it2cf_zlu treatment in h5_read for bwd compat


## Version 1.1.0

This update to cppdlr adds significant new functionality, including symmetrized DLR grids, and improvements to ensure accuracy of DLR expansions in L^2(tau) norm (both described in the updated documentation).

We thank all contributors: Thomas Hahn, Alexander Hampel, Jason Kaye, Henri Menke, Hugo U. R. Strand, Nils Wentzell

### New features
* Symmetrized DLR grids, w/ tests to compare symmetrized/non-symmetrized grids
* Documentation updates: symmetrized grids, discussion of error
* Imaginary time quadrature weighting for DLR frequency selection to ensure accuracy in L^2(tau) norm, updated tests accordingly
* Function to compute imaginary time inner products of DLR expansions
* Program to print DLR ranks for various Lambda
* Tests to interpolate in imag time and measure error in imag freq, and vice versa
* Unit tests print errors via fmt
* Implementation of symmetrized pivoted Gram-Schmidt with tolerance or rank specified
* Print warning if user chooses epsilon dangerously small
* Expose fermionic and bosonic Matsubara frequency kernels to user directly

### Bug fixes
* Fixed bosonic kernel at i*om_n=0
* Fixed bug involving non-contiguous views in vals2coefs implementations
* Fixed bug in error checking of fine discretization of imag time kernel in geterr_k_it
* Fixed uninitialized niom value in dlr_imfreq constructor
* Fixed range bug in imtime_ops.interp_matrix_sym_bos

### Optimizations
* Make complex copy of it2cf.lu in imtime_ops constructor to avoid on-the-fly copies
* imtime_ops::convolve optimizations; use nda::matmul when possible

### Other changes
* Various code simplifications
* Test tolerance adjustments
* References to cppdlr preprint in documentation and readme
* Merge latest changes of app4triqs skeleton
* Add ubuntu-intel build to jenkins


## Version 1.0.0

This is the initial release of cppdlr, a C++ library implementing the discrete Lehmann representation of imaginary time Green's functions.

We thank all contributors: Thomas Hahn, Jason Kaye, Hugo U. R. Strand, Nils Wentzell
