(changelog)=

# Changelog

## Version 1.0.0

This is the initial release of cppdlr, a C++ library implementing the discrete Lehmann representation of imaginary time Green's functions.

We thank all contributors: Thomas Hahn, Jason Kaye, Hugo U. R. Strand, Nils Wentzell

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

### General
* Remove redundant sys.path insertion in conf.py.in
* Remove pip3 install command from Dockerfile
* Update requirements.txt with dependencies for documentation builds
* Disable notes about C++ ABI changes when using gcc
* Fix compiler warning, use abs from std namespace
* Update macos build instructions in Jenkinsfile and build.yml
* Set proper GNU install dirs also in env vars file
* Fix typo in Jenkinsfile
* Bump actions/checkout from 2 to 4
* Fix settings environment variables
* Allow manual dispatch and triggering action from other workflow
* Prevent unintentional parallelization in OpenBLAS
* Use ccache to speed up compilation
* Remove numpydoc sources from doc/sphinxext
* Raise shm size for docker run commands to comply with mkl requirements
* Skip image and binary files in replace_and_rename
* Use python3 instead of python2 in replace_and_rename.py script


### cmake
* Bump version number to 1.1.0
* Bump Version number to 3.3
* Remove C as a project language
* Update top-level CMakeLists.txt file with latest app4triqs skeleton
* Fix inclusion of targets file in cppdlr-config.cmake.in
* Do not define DEBUG macros for RelWithDebInfo builds
* Correct target file inclusion PATH to be absolute
* Only use GNUInstallDirs for LIBDIR
* Fix target inclusion directory to use GNUInstallDirs
* Use unstable branch of cpp2py
* Minor improvements in top-level CMakeLists.txt
* Consistently use GNUInstallDirs for install commands
* Use GNUInstallDirs to obtain installation directories
* Set policy 114 to new
* Set policy CMP0144 to new
* Run Debug checks also in RelWithDebInfo build mode

### actions
* Bump actions/cache restore/save to version 4

### ghactions
* Synchronize build file with app4triqs, add caching using ccache
* Use libc++ for clang builds

### jenkins
* Add ubuntu-intel build
* For sanitized build use RelWithDebInfo build mode

### build
* add packaging directory to cmake
* automatically set version in packaging