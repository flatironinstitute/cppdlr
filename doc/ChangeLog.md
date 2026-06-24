(changelog)=

# Changelog

## Version 1.3.0

This update contains additional functionality, bug fixes, and build system improvements.

We thank all contributors: Thomas Hahn, Alexander Hampel, Jason Kaye, Henri Menke, Dylan Simon, Hugo U. R. Strand, Nils Wentzell, Yang Yu

### New features
* Add free functions `build_it2if` and `build_if2it` to construct the matrices mapping DLR imaginary time grid values to imaginary frequency grid values and vice versa
* Add `imtime_ops::convmat_inplace` builder, splitting `convmat` into separate allocation and in-place construction steps (#13)
* Support non-square matrix-valued Green's functions in `imtime_ops::convolve`
* Add serialization and deserialization support to `imfreq_ops` and `imtime_ops`

### Bug fixes
* Fix bug in the imaginary frequency to imaginary time matrix in the bosonic, symmetrized case

### Other changes
* Sign-canonicalize the DLR real-frequency and imaginary-frequency grids for deterministic grid selection
* Deprecate the redundant `statistic` parameter in `convolve`, `convmat`, and `convmat_inplace`; the new signatures omit it (#17)
* Update to account for interface changes in nda and remove `using namespace nda`
* Bump bundled `fmt` to version 12.0.0
* Various CMake fixes (`link_directories` for the `cppdlr_c` target, doc-build custom commands)
* Modernize and synchronize GitHub Actions and Jenkins CI configuration
* Add JOSS citation information and banner to the README
* Documentation and test tolerance updates


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
