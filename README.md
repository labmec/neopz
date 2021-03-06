# NeoPZ

NeoPZ is a open-source, cross-platform, finite element library
written in C++.

## Features
- Discontinuous, *H1, H(div)* and *H(curl)* conforming approximation spaces
- *hp*-adaptivity
- Non-linear geometric maps
- Refinement patterns

And much more!

## Requirements
- A C++ 17 compiler
- [CMake](https://cmake.org/download/) 3.11.0+
- A visualization software to display the VTK files outputted by NeoPZ.
 [Paraview](https://www.paraview.org/download/) is recommended.

### Optional external libraries
The usage of NeoPZ can be improved by linking against the following libraries:
- [log4cxx](https://logging.apache.org/log4cxx/latest_stable/), for logging information
- [Intel MKL](https://software.intel.com/en-us/mkl), for optimized linear system solvers
- etc.

## Installation
To install NeoPZ please refer to the tutorial for [Linux](http://www.labmec.org.br/wiki/howto/roteiro_pzlinux_eng), [Windows](http://www.labmec.org.br/wiki/howto/pz_windows) or [MacOS](http://www.labmec.org.br/wiki/howto/pz_no_ubuntu_macosx_mavericks).

## NeoPZ documentation

A Doxygen documentation can be found 
[here](http://www.labmec.org.br/pz/arquivos-html/html/index.html).

## Publication

Bernard Devloo, P.R., 1997. PZ: An object oriented environment
for scientific programming. Computer Methods in Applied Mechanics
and Engineering 150, 133–153.
https://doi.org/10.1016/s0045-7825(97)00097-2
