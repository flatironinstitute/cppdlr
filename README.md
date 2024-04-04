# `cppdlr`: Imaginary time calculations using the discrete Lehmann representation

[![build](https://github.com/flatironinstitute/cppdlr/workflows/build/badge.svg?branch=main)](https://github.com/flatironinstitute/cppdlr/actions?query=workflow%3Abuild)

Authors: [Jason Kaye](https://users.flatironinstitute.org/~jkaye/), [Nils
Wentzell](https://github.com/Wentzell), and [Hugo U. R.
Strand](https://github.com/HugoStrand) (2023)

`cppdlr` is a C++ library implementing the discrete Lehmann representation (DLR) of
imaginary time single-particle Green's functions and other imaginary time
quantities. For more information on the DLR, see the references below.

Please see [the documentation](https://flatironinstitute.github.io/cppdlr/) for more information about `cppdlr`.

While we recommend using the latest stable release of this library, which is the default branch of this repository, the `main` branch is treated as a development branch and might have new (unstabilized) features. The most up-to-date documentation corresponding to this development branch can be found [here](https://flatironinstitute.github.io/cppdlr/main/).

Libraries implementing the DLR are available in other languages:

- Python, via [pydlr](https://github.com/jasonkaye/libdlr)
- Fortran, via [libdlr](https://github.com/jasonkaye/libdlr)
- Julia, via [Lehmann.jl](https://github.com/numericaleft/Lehmann.jl)

## References and citation

If you use `cppdlr` in your software or published research works, please mention
this, and cite the following references. Citations help to encourage the
development and maintainence of open-source scientific software.

- The original reference on the DLR: [J. Kaye, K. Chen, O. Parcollet, "Discrete Lehmann representation of imaginary time Green's functions", Phys. Rev. B, 105, 235115 (2022).](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.235115) \[[arXiv:2107.13094](https://arxiv.org/abs/2107.13094)\]
- The companion paper to `cppdlr`: [J. Kaye, H. U. R. Strand, N. Wentzell, "cppdlr: Imaginary time calculations using the discrete Lehmann representation", arXiv:2404.02334 (2024).](https://arxiv.org/abs/2404.02334)

## Contact

Please email jkaye@flatironinstitute.org with questions, or post an [issue](https://github.com/flatironinstitute/cppdlr/issues).

## License

`cppdlr` is licensed under the Apache License, Version 2.0, for more information see the `LICENSE` file.
