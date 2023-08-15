.. _main:

``cppdlr``: Imaginary time calculations using the discrete Lehmann representation
=================================================================================

``cppdlr`` is a C++ library implementing the discrete Lehmann representation (DLR) of
imaginary time single-particle Green's functions and other imaginary time
quantities. It can be downloaded `from GitHub <https://github.com/flatironinstitute/cppdlr>`_.

Getting started
---------------

  - Install the library by following the instructions on the :ref:`installation page<install>`.
  - If you are unfamiliar with the DLR, you can read the :ref:`background page<background>` for a
  brief introduction, or the references below for more detail. If you are
  familiar with the DLR, you should still take a brief look at the
  :ref:`background page<background>` for important information about the library.
  - To begin using `cppdlr` in your application, the :ref:`examples page<examples>` provides a
  good starting point.
  - If you need more specific information about `cppdlr` classes and functions,
  take a look at the :ref:`reference documentation page<documentation>`.

References
----------

If you use ``cppdlr`` in your software or published research works, please cite
our repository and one, or
both, of these references. Citations help to encourage the development and
maintainence of open-source scientific software.

- The original reference on the DLR: `Discrete Lehmann representation of imaginary time Green's functions, Jason Kaye, Kun Chen, and Olivier Parcollet, Phys. Rev. B 105, 235115, 2022.
  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.235115>`_
  [`arXiv:2107.13094 <https://arxiv.org/abs/2107.13094>`_]
- The companion paper to `libdlr <https://github.com/jasonkaye/libdlr>`_, which
  contains a briefer overview of the DLR: `libdlr: Efficient imaginary time calculations using the discrete
  Lehmann representation, Jason Kaye, Kun Chen, and Hugo U.R. Strand, Comput.
  Phys. Commun. 280, 108458, 2022.
  <https://www.sciencedirect.com/science/article/pii/S0010465522001771>`_
  [`arXiv:2110.06765 <https://arxiv.org/abs/2110.06765>`_]

Related libraries
-----------------

Libraries implementing the DLR are available in other languages:

- Python, via `pydlr <https://github.com/jasonkaye/libdlr>`_ 
- Fortran, via `libdlr <https://github.com/jasonkaye/libdlr>`_
- Julia, via `Lehmann.jl <https://github.com/numericaleft/Lehmann.jl>`_
 
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install.rst
   background.rst
   examples.rst
   documentation.rst
   issues.rst
