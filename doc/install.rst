.. highlight:: bash

.. _install:

Installation
============

Compiling from source
---------------------

.. note:: To guarantee reproducibility in scientific calculations we strongly recommend the use of a stable `release <https://github.com/TRIQS/triqs/releases>`_ of both TRIQS and its applications.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/cppdlr`` repository from GitHub::

     $ git clone https://github.com/flatironinstitute/cppdlr cppdlr.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir cppdlr.build && cd cppdlr.build

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake -DCMAKE_INSTALL_PREFIX=path_to_install_dir ../cppdlr.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Versions
--------

To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd cppdlr.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``cppdlr`` can be configured using CMake-options::

    cmake ../cppdlr.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_cppdlr         |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+

Dependencies
--------------------------

The dependencies of ``cppdlr`` are as follows:

* gcc version 12 or later OR clang version 15 or later
* BLAS/LAPACK
* hdf5
* mpi

``hdf5`` and ``mpi`` are required by the ``nda`` library, which is itself built
automatically with ``cppdlr``.

If you wish to build the documentation, the dependencies also include:

* libclang
* python

and the Python packages

* sphinx
* nbsphinx
* myst_parser
* sphinx_rtd_theme
* linkify-it-py
* clang
