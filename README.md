# NeoPZ

NeoPZ is a open-source, cross-platform, finite element library
written in C++ with support for many advanced finite element techniques.

## Features
- Discontinuous, *H1, H(div)* and *H(curl)* conforming approximation spaces
- Multiphysics support
- *hp*-adaptivity
- Non-uniform meshes with support for hanging nodes
- Non-linear geometrical mappings (curved elements with exact representation)
- Run-time defined geometrical refinement patterns
- Forward automatic differentiation

And much more!

## Requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.13.0+

### Optional external libraries
- [Paraview](https://www.paraview.org/download/) is recommended for display of VTK output files.

The usage of NeoPZ can be improved by linking against the following libraries:

- [log4cxx](https://logging.apache.org/log4cxx/latest_stable/), for logging information for debugging purposes.
- [Intel MKL](https://software.intel.com/en-us/mkl), for enabling sparse matrices solvers (in-house algorithms for skyline matrices are available, among other matrix storage formats).
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) for (experimental) support of the [BDDC  technique](https://epubs.siam.org/doi/abs/10.1137/S1064827502412887?journalCode=sjoce3).
- [TBB](https://github.com/oneapi-src/oneTBB), also used for experimental techniques such as BDDC
- [blaze](https://bitbucket.org/blaze-lib/blaze), needed for projects using [SBFEM](https://www.cies.unsw.edu.au/scaled-boundary-finite-element-method-2a)
- [LAPACK](http://www.netlib.org/lapack/), for eigenvalues computations in full or banded matrices. If enabled, it is also internally used occasionaly instead of in-house linear algebra algorithms.
- [BLAS](http://www.netlib.org/blas/), replaces some in-house linear algebra algorithms.
- [Boost](https://www.boost.org/), for building unit tests and experimental techniques.

## Configuration and Install
The NeoPZ library uses CMake for configuring and installing the library. As a simple example, on UNIX systems, this could be done as:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DUSING_LOG4CXX=ON -DUSING_MKL=ON ..
make -j && sudo make install
```
where some CMake options were given as an illustration. The following CMake options are available:

- `REAL_TYPE`: which floating point type to use for real data types (geometry)
- `STATE_TYPE`: which floating point type to use for the Finite Element matrices. It can be real or complex.
- `CMAKE_INSTALL_PREFIX`: where to install the NeoPZ library. Defaulted to `/opt/neopz` on UNIX systems.
- `USING_BLAZE`: Enable blaze library support
- `USING_BOOST`: Enable Boost libraries (`Boost::date_time`, `Boost::graph` and `Boost::unit_test_framework`) support
- `USING_BLAS` : Enable BLAS library support
- `USING_LAPACK`: Enable LAPACK library support
- `USING_MKL`: Enable MKL library support
- `USING_LOG4CXX`: Enable logging information through log4cxx support.
- `USING_METIS`: Enable Metis library support.
- `CMAKE_BUILD_TYPE`: `Release`, `Debug` or `RelWithDebugInfo`.
- `BUILD_PLASTICITY_MATERIALS`: to include support for classes modelling plastic materials.


*Internal options/deprecated options*
The remaining `USING_XXX` CMake options are for internal usage only.
The following options are listed for completeness:

- `BUILD_BROKEN`: if you are looking for adventures.
- `BUILD_PERF_TESTS`: for building a few performance tests. The tests are in need of a revision.
- `BUILD_TUTORIAL`: soon to be deprecated.Legacy code.
- `BUILD_PROJECTS`: soon to be deprecated. Legacy code.
- `BUILD_PUBLICATIONS`: build legacy projects with code used to generate results that were published.


For more info on how to install NeoPZ, please refer to the tutorial for [Linux](http://www.labmec.org.br/wiki/howto/roteiro_pzlinux_eng), [Windows](http://www.labmec.org.br/wiki/howto/pz_windows) or [MacOS](http://www.labmec.org.br/wiki/howto/pz_no_ubuntu_macosx_mavericks).


## NeoPZ documentation

A Doxygen documentation can be found 
[here](http://www.labmec.org.br/pz/arquivos-html/html/index.html).

## How to cite NeoPZ

Devloo, P.R. B., 1997. PZ: An object oriented environment
for scientific programming. Computer Methods in Applied Mechanics
and Engineering 150, 133–153.
https://doi.org/10.1016/s0045-7825(97)00097-2

```
@article{neopz97,
author = {Devloo, Philippe Remy Bernard},
year = {1997},
month = {12},
pages = {},
title = {PZ: An object oriented environment for scientific programming},
volume = {150},
booktitle = {Computer Methods in Applied Mechanics and Engineering},
doi = {https://doi.org/10.1016/s0045-7825(97)00097-2}
}
```
