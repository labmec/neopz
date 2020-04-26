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
- A C++ 11 compiler
- [CMake](https://cmake.org/download/) 3.11.0+
- A visualization software to display the VTK files outputted by NeoPZ.
 [Paraview](https://www.paraview.org/download/) is recommended.
- [pthreads for Win32](http://sourceware.org/pthreads-win32/) 
(if running on Windows)

### Optional external libraries
The usage of NeoPZ can be improved by linking against the following libraries:
- [log4cxx](https://logging.apache.org/log4cxx/latest_stable/), for logging information
- [Boost](https://www.boost.org/) TODO
- [Intel MKL](https://software.intel.com/en-us/mkl), for optimized linear system solvers
- etc.

## Installation
To install NeoPZ run the following commands:
```shell
git clone https://github.com/labmec/neopz.git
mkdir neopz-build
cd neopz-build
ccmake ../neopz
```
This will open CMake configuration window, in which you can:
- Set the installation prefix path,
- Set the build type (Debug, Release, RelWithDebInfo),
- Toggle the external libraries to be linked against NeoPZ, 
- Choose to install executable targets under the ``BUILD_PROJECTS`` option.
- etc.

Then press ``c`` to 'configure', ``g`` to 'generate' the Makefile and install 
NeoPZ with:
```shell
make install
```
\* The last command may require ``sudo`` according to the chosen installation path.

## Using NeoPZ in your project

TODO

## NeoPZ documentation

A Doxygen documentation can be found 
[here](http://www.labmec.org.br/pz/arquivos-html/html/index.html).

## Publication

Bernard Devloo, P.R., 1997. PZ: An object oriented environment
for scientific programming. Computer Methods in Applied Mechanics
and Engineering 150, 133–153.
https://doi.org/10.1016/s0045-7825(97)00097-2
