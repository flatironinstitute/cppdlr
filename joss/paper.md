---
title: 'cppdlr: Imaginary time calculations using the discrete Lehmann representation'
tags:
  - C++
  - quantum many-body systems
  - imaginary time Green's function
  - Matsubara Green's function
  - many-body Green's function methods
  - low-rank compression
authors:
  - name: Jason Kaye
    orcid: 0000-0001-8045-6179
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
    corresponding: true # (This is how to denote the corresponding author)
  - name: Nils Wentzell
    orcid: 0000-0003-3613-007X
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Hugo U. R. Strand
    orcid: 0000-0002-7263-4403
    affiliation: "3, 4"
affiliations:
 - name: Center for Computational Quantum Physics, Flatiron Institute, New York, NY, USA
   index: 1
 - name: Center for Computational Mathematics, Flatiron Institute, New York, NY, USA
   index: 2
 - name: School of Science and Technology, Örebro University, Örebro, Sweden
   index: 3
 - name: Institute for Molecules and Materials, Radboud University, 6525 AJ Nijmegen, the Netherlands
   index: 4
date: 22 September 2023 
bibliography: paper.bib

---

# Summary

Imaginary time Green's functions encode the static and dynamical response of quantum systems at thermal equilibrium to external perturbations, such as applied electromagnetic fields. They therefore represent a direct point of connection between theoretical calculations and experimental measurements. As a consequence, they appear routinely in quantum many-body calculations at finite temperature, both for model systems like the Hubbard model `[@hubbard63]`, and in ab-initio electronic structure calculations beyond density functional theory, e.g., using Hedin's GW method `[@hedin65,@golze19]`.
Highly compact and accurate representations of imaginary time Green's functions and related imaginary time-dependent response functions are therefore an important ingredient in the development of robust and efficient codes for quantum many-body calculations.
 However, obtaining such representations has traditionally been challenging, particularly for low temperature calculations, in which the functions develop steep gradients. 

In the past several years, significant progress has been achieved using low-rank approximations of the
spectral Lehmann representation, which is given by
$$G(\tau) = - \int_{-\infty}^\infty d\omega \,
\frac{e^{-\tau \omega}}{1 + e^{-\beta \omega}} \, \rho(\omega).$$ 
Here, $G(\tau)$ is a fermionic single-particle imaginary time Green's function, and
$\rho(\omega)$ is its corresponding spectral function, which encodes information
about the single-particle excitations of the underlying quantum many-body system.
The spectral function always exists, but
is typically not known. However, the existence of this integral representation
constrains the space of possible imaginary time Green's functions
to lie within the image of the integral operator, which is numerically low-rank,
enabling the construction of highly compact basis representations. The
intermediate representation (IR) was introduced first, and used the singular value
decomposition to obtain an orthogonal but non-explicit basis of imaginary time
Green's functions `[@shinaoka17,@chikano18]`. The recently-introduced discrete Lehmann representation (DLR) uses the
interpolative decomposition to obtain a non-orthogonal basis consisting of known
exponential functions `[@kaye22_dlr]`. The number of basis functions required in both representations is similar, and typically significantly less than the previous state-of-the art methods based on orthogonal polynomials `[@bohenke11,@gull18,@dong20]`.

The DLR's use of an explicit basis of simple functions makes many common operations,
including interpolation, integration, Fourier transform, and convolution, simple and
highly efficient. This has led to a variety of recent algorithmic advances, including in
reducing the size of the Matsubara frequency mesh in dynamical mean-field theory
calculations `[@sheng23]`, stabilizing the calculation of the single-particle self-energy via the Dyson
equation `[@labollita23]`, improving the efficiency of the imaginary time discretization in the
mixing Green's function of the Keldysh formalism `[@kaye23_eqdyson]`, and accelerating the
evaluation of imaginary time Feynman diagrams `[@kaye23_diagrams]`. It has also
yielded immediate applications in computational physics, for example in low-temperature
studies of superconductivity `[@cai22,@hou23,@johnston23]`. The DLR can be straightforwardly
integrated into existing algorithms and codes, often yielding significant
improvements in efficiency, accuracy, and algorithmic simplicity.

# Statement of need

`cppdlr` is a C++ library which constructs the DLR and implements its standard operations. The
flexible yet high-level interface of `cppdlr` makes it appealing for use both 
in small-scale applications and in
existing large-scale software projects.
The DLR has previously been implemented in other programming languages,
specifically in Python via `pydlr`, in Fortran via `libdlr`, and in Julia via
`Lehmann.jl` `[@kaye22_libdlr,@pydlr,@libdlr,@Lehmann.jl]`, as well as in the `sparse-ir` library implementing the IR `[@wallerberger23]`. `cppdlr` nevertheless provides a needed platform for
future developments. First, `cppdlr` is written in C++, a common language used
by many large projects in the quantum many-body physics community. Second, it
offers a high-level user interface simpler than that of `libdlr`, enabled by
the use of C++ templating and the `nda` library `[@nda]` for array types and BLAS/LAPACK
compatibility. These features have, for example, enabled the
implementation of the DLR in the TRIQS library `[@parcollet15]`
for quantum many-body calculations.

`cppdlr` is distributed under the Apache License Version 2.0 through a public Git repository `[@cppdlr_git]`. The project documentation `[@cppdlr_doc]` is extensive, containing background on the DLR, a user guide describing example programs packaged with the library, and application interface (API) reference documentation for all classes and functions. We envision `cppdlr` as a platform for future algorithmic developments involving the DLR, and as a go-to tool for applications employing the DLR.

# Acknowledgements

We are thankful for helpful discussions with Kun Chen, Olivier Parcollet, Malte Rösner, and Yann in 't Veld. H.U.R.S. acknowledges funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (Grant agreement No. 854843-FASTCORR). The Flatiron Institute is a division of the Simons Foundation. 
