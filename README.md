# NeoPZ — SBFEM crack-tip enrichments

This repository, on the `paper/gfem-sbfem-enrichments` branch, accompanies the paper:

N. Shauer, C. Song, C. A. Duarte, and P. R. B. Devloo,
*A Well-Conditioned and Optimally Convergent SBFEM-Enriched GFEM for
Two-Dimensional Fracture in Orthotropic Media and Interfaces*.
Manuscript submitted for publication.

**Generate material-dependent crack-tip functions and evaluate them in your own code.**

This research branch, `paper/gfem-sbfem-enrichments`, provides
[`TPZGFemEnrichment`](Pre/TPZGFemEnrichment.h), a local SBFEM generator for
two-dimensional elastic crack-tip fields. It supports homogeneous isotropic and
orthotropic materials, as well as material sectors for interface-crack problems.
**No separate GFEM repository is required.**

The generator builds a circular domain centered at the crack tip, computes its
SBFEM modes, and evaluates complete crack-tip functions and their Cartesian
derivatives. The functions include the radial factor, not just the boundary
eigenvector. Your application supplies the global mesh, partition of unity,
enrichment coefficients, integration, and any shifting or stabilization.

General library information is retained in the
[general NeoPZ README](docs/README-NeoPZ.md) and
[NeoPZ documentation](https://labmec.github.io/neopz).

## Contents

- [Installation](#installation)
- [Run the examples](#run-the-examples)
- [Run unit tests](#run-unit-tests)
- [Use in your own code](#use-in-your-own-code)
- [JSON configuration](#json-configuration)
- [Conventions and limitations](#conventions-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Docker image](#docker-image)

## Installation

### Dependencies

- A C++17 compiler, CMake >= 3.14, and Ninja or Make.
- BLAS/LAPACK for the SBFEM eigenproblem; macOS uses Apple's Accelerate framework.
- [nlohmann JSON](https://json.nlohmann.me/) >= 3.2.0 for configuration and the public API.
- Catch2 v3 for unit tests only. NeoPZ fetches it if not found, requiring network access during configuration.

Eigen, MKL, autodiff, and a separate GFEM library are **not required** for this
generator. The documented configuration uses `double` for both REAL and STATE.

**Ubuntu / Debian:**

```sh
sudo apt update
sudo apt install build-essential git cmake ninja-build libblas-dev liblapack-dev nlohmann-json3-dev
```

Package references: [JSON](https://packages.ubuntu.com/noble/nlohmann-json3-dev)
and [LAPACK](https://packages.ubuntu.com/noble/liblapack-dev).

**macOS**, using Apple's Command Line Tools and Homebrew:

```sh
xcode-select --install
brew install cmake ninja nlohmann-json
```

Skip the first command if the tools are already installed. No separate LAPACK
installation is needed for this macOS configuration. Other package managers
are described in the [JSON installation guide](https://json.nlohmann.me/integration/package_managers/).

### Obtain this branch

Replace `REPOSITORY_URL` with the clone URL of the repository or fork hosting
this research branch; it may not yet be available in upstream NeoPZ.

```sh
git clone --branch paper/gfem-sbfem-enrichments REPOSITORY_URL neopz
cd neopz
```

### Configure, build, and install

Run from the repository root:

```sh
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/../neopz_install" \
  -DREAL_TYPE=double -DSTATE_TYPE=double \
  -DUSING_LAPACK=ON \
  -DUSING_SBFEM_ENRICHMENT=ON \
  -DBUILD_SBFEM_ENRICHMENT_EXAMPLES=ON \
  -DBUILD_UNITTESTING=ON \
  -DBUILD_PLASTICITY_MATERIALS=OFF

cmake --build build --target pz tinyfad_generate sbfem_eigenvalues sbfem_evaluate --parallel 4
cmake --install build
```

`tinyfad_generate` is an existing helper needed by NeoPZ's installation rules,
not an enrichment dependency. The chosen prefix is writable without `sudo`.

`USING_SBFEM_ENRICHMENT` is OFF by default so ordinary NeoPZ builds do not acquire
a JSON requirement. Here it must be ON. The tutorials are independently enabled
with `BUILD_SBFEM_ENRICHMENT_EXAMPLES`. Set `BUILD_UNITTESTING=OFF` to omit the
unit tests and Catch2; the example checks remain available.

If Homebrew packages are not found, add `-DCMAKE_PREFIX_PATH="$(brew --prefix)"`
to configuration. Alternatively set `nlohmann_json_DIR` to the directory
containing `nlohmann_jsonConfig.cmake`. Use a fresh build directory when changing
compilers or dependency installations.

## Run the examples

Two executables are provided in [Examples/SBFemEnrichment](Examples/SBFemEnrichment).
The commands below assume the Ninja build above and execution from the
repository root. NeoPZ may print progress before the example's results.

### Selected eigenvalues for an isotropic crack

```sh
./build/Examples/SBFemEnrichment/sbfem_eigenvalues \
  Examples/SBFemEnrichment/isotropic.json --verify-isotropic
```

The program prints selected SBFEM eigenvalues `s` and displacement exponents
`lambda = -s`. Expect two leading homogeneous isotropic modes with `s≈-0.5`,
giving displacement behavior `xi^0.5`, where `xi=r/R`. The verification flag
checks their count and exponents; omit it for other materials or exponent ranges.
This reports the **selected crack-tip modes**, not the entire SBFEM spectrum.

### Values and Cartesian derivatives at a point

```sh
./build/Examples/SBFemEnrichment/sbfem_evaluate \
  Examples/SBFemEnrichment/isotropic.json 0.35 0.20 --verify-gradients
```

For each zero-based mode index `j`, the output is:

```text
j  phi_x  phi_y  dphi_x/dx  dphi_x/dy  dphi_y/dx  dphi_y/dy
```

The point `(0.35,0.20)` is in the generator's tip-centered coordinates.
`--verify-gradients` compares derivatives with centered finite differences.
Use this check away from the crack, tip, sector boundaries, and material
interfaces; never differentiate across a jump. Omit it for ordinary sampling.

Both programs accept a JSON file. With no arguments, they use `isotropic.json`
in the current working directory; evaluation defaults to `(0.35,0.20)`.
Use `--help` for syntax.

### Other materials

```sh
./build/Examples/SBFemEnrichment/sbfem_evaluate \
  Examples/SBFemEnrichment/orthotropic.json 0.35 0.20

./build/Examples/SBFemEnrichment/sbfem_eigenvalues \
  Examples/SBFemEnrichment/bimaterial.json
```

These are small demonstrations, not reproductions of every paper configuration.

## Run unit tests

```sh
cmake --build build --target \
  Test_TPZGFEMEnrichment Test_TPZGFEMMaterials Test_TPZGFEMRegression \
  TestSBFem TestSBFemHdiv --parallel 4

ctest --test-dir build -R 'Test_TPZGFEM|^TestSBFem|^SBFemExample' --output-on-failure
```

Coverage includes material setup, complex reconstruction, optional modal
transformations, Cartesian gradients, crack-face selection, angular cuts,
Williams exponent ranges, and defaults. The example checks verify leading
isotropic exponents, finite-difference gradients for the three supplied
materials, and rejection of evaluation at the tip. These are local tests, not
global GFEM convergence/conditioning studies; see [test details](UnitTest_PZ/TestGFemEnrichment/README.md).

For the complete NeoPZ suite, build all test targets before using unfiltered
`ctest`; additional tests may require additional features. Offline test builds
can use `-DFETCHCONTENT_SOURCE_DIR_CATCH2=/path/to/Catch2` with a v3 checkout.

## Use in your own code

An installed-package consumer links only `NeoPZ::pz`; include directories and
JSON are supplied transitively:

```cmake
cmake_minimum_required(VERSION 3.14)
project(MyEnrichmentCode LANGUAGES CXX)
find_package(NeoPZ CONFIG REQUIRED)
if(NOT NeoPZ_SBFEM_ENRICHMENT)
    message(FATAL_ERROR "Install NeoPZ with USING_SBFEM_ENRICHMENT=ON")
endif()
add_executable(my_enrichment main.cpp)
target_compile_features(my_enrichment PRIVATE cxx_std_17)
target_link_libraries(my_enrichment PRIVATE NeoPZ::pz)
```

A minimal `main.cpp`:

```cpp
#include "TPZGFemEnrichment.h"
#include <fstream>
#include <iostream>

int main()
{
    std::ifstream input("isotropic.json");
    if(!input) return 1;
    nlohmann::json config;
    input >> config;

    TPZGFemEnrichment enrichment;
    enrichment.LoadFromJSON(config); // Solve the local eigenproblem once.
    const auto eigenvalues = enrichment.SelectedEigenvalues();

    std::vector<REAL> point = {0.35, 0.20};
    std::vector<std::array<REAL,2>> phi;
    std::vector<std::array<REAL,4>> dphi;
    enrichment.ComputeEnrichmentWithDerivatives(
        TPZGFemEnrichment::EAny, point, phi, dphi);

    for(size_t j=0; j<phi.size(); ++j)
        std::cout << j << " s=" << eigenvalues[j]
                  << " phi_x=" << phi[j][0] << " phi_y=" << phi[j][1]
                  << " dphi_x/dx=" << dphi[j][0] << '\n';
}
```

Configure with `-DNeoPZ_DIR=/path/to/neopz_install/lib/cmake/neopz` and copy
`isotropic.json` into the program's working directory. See the complete
[`evaluate.cpp`](Examples/SBFemEnrichment/evaluate.cpp) for argument handling,
finite-value checks, and derivative verification.

The example folder is also an independent CMake project. After installation,
run from the repository root:

```sh
cmake -S Examples/SBFemEnrichment -B build-consumer -G Ninja \
  -DNeoPZ_DIR="$PWD/../neopz_install/lib/cmake/neopz"
cmake --build build-consumer --parallel 4
ctest --test-dir build-consumer --output-on-failure
./build-consumer/sbfem_evaluate build-consumer/isotropic.json 0.35 0.20
```

This consumer needs no source-tree headers or extra GFEM sources. You can copy
the folder into your own project.

### Array layout and coordinates

| Expression | Meaning |
|---|---|
| `phi[j][0]`, `phi[j][1]` | x and y displacement components of real function `j` |
| `dphi[j][0]`, `dphi[j][1]` | x and y derivatives of its x component |
| `dphi[j][2]`, `dphi[j][3]` | x and y derivatives of its y component |
| `SelectedEigenvalues()[j]` | Original eigenvalue associated with function `j` |
| `getNEigenvalues()` | Number of real output functions, not distinct eigenvalues |

Derivatives are with respect to physical Cartesian coordinates, not `xi` or
`eta`; they are not strains. Tensor shear strain for a field is
`0.5*(dphi[j][1]+dphi[j][2])`.

The tip is at `(0,0)`. For a physical tip `(x_tip,y_tip)` with the same coordinate
axes, pass `(x-x_tip,y-y_tip)`. Orient the crack/material axes consistently in
JSON. If rotating the coordinate frame, transform values and gradients
consistently. The evaluator does not assign GFEM coefficients to components.

Generate once, then reuse the object at many points: evaluation does not solve
the eigenproblem again. `GetEnrichmentFunction()` returns a callback; the owner
must outlive all callback calls.

## JSON configuration

JSON describes the **local generator**, not the global GFEM mesh or loading.
Keys with spaces, spelling, and case matter. The examples use strict JSON:
numeric radians, no `pi` expressions, comments, or trailing commas.

The complete [isotropic configuration](Examples/SBFemEnrichment/isotropic.json):

```json
{
  "problem_type": "Elastic",
  "radius": 1.0,
  "angles": [-3.141592653589793, -1.570796326794897, 0.0, 1.570796326794897, 3.141592653589793],
  "sector material ids": [1, 1, 1, 1],
  "skeleton porder": 6,
  "nref_skeleton": 1,
  "materials": [
    {
      "matid": 1,
      "tension state": "plane strain",
      "isotropic": {"young": 1.0, "poisson": 0.3}
    }
  ],
  "boundary_conditions": [],
  "eigenvalue exponent ranges": [[0.01, 0.9]],
  "complex mode reconstruction": "representative_parts"
}
```

### Geometry and discretization

| Key | Meaning |
|---|---|
| `problem_type` | Use `"Elastic"` for the documented two-component evaluator |
| `radius` | Positive local-circle radius `R`; **not** the global geometrical enrichment radius |
| `angles` | Increasing sector endpoints spanning `2*pi`; first/last rays are the crack faces |
| `sector material ids` | One material ID per consecutive angular interval |
| `skeleton porder` | Boundary interpolation order `q`; default 6 |
| `nref_skeleton` | Binary skeleton refinement levels; default 1 |

Four initial sectors with one refinement give eight boundary elements of order
6. The boundary is circular, not polygonal. Coincident endpoint nodes at `-pi`
and `pi` remain separate to represent the crack faces. The example crack lies
along negative x; shift all angular endpoints consistently to incline it.

### Materials

Each sector ID needs a corresponding material entry. Use consistent units; the
supplied properties are dimensionless demonstrations.

- **Isotropic:** `isotropic.young` and `isotropic.poisson` are `E` and `nu`.
- **Plane state:** `tension state` is exactly `"plane stress"` or `"plane strain"`.
- **Orthotropic:** specify `Ex`, `Ey`, `nuxy`, and `Gxy` inside `orthotropic`.
  Its optional nonzero two-component `direction` defines the material x-axis
  via `atan2(direction[1],direction[0])` in the generator's frame. See the
  [45-degree example](Examples/SBFemEnrichment/orthotropic.json).
- **Bimaterial:** use different IDs for different sectors. In the
  [supplied example](Examples/SBFemEnrichment/bimaterial.json), `[1,1,2,2]`
  assigns material 1 below `y=0` and material 2 above it: `E=1` and `10`,
  respectively, with `nu=0.3` and plane stress in each.

Use IDs 1 and 2 as in the examples. Internal defaults `-1`, `4`, `5`, and
auxiliary IDs generated with offset 1000 must not collide with physical IDs.
Advanced keys `skeletonMatid` and `point_singularity` change internal skeleton
and tip IDs; leave them alone unless you understand the mesh construction.

`boundary_conditions: []` is sufficient for these local mode-generation
examples. It does **not** mean the eventual global GFEM problem has no BCs;
these tutorials do not solve that loaded boundary-value problem.

### Mode selection and complex reconstruction

`eigenvalue exponent ranges` selects **open** intervals of `abs(real(s))`.
The stored SBFEM eigenvalue is `s`; evaluation uses `lambda=-s` in `xi^lambda`.
The default `[[0.01,0.9]]` includes the leading `lambda≈0.5` modes. For isotropic
Williams fields, `[[1.4,1.6]]` selects the next pair and
`[[0.4,0.6],[1.4,1.6]]` selects both pairs, subject to discretization accuracy.

`complex mode reconstruction` accepts:

- `"representative_parts"` (default): take real/imaginary parts of one complex representative to obtain two real functions.
- `"conjugate_sum"`: use explicit conjugate sums/differences. The scale or coordinates can differ; equality of spaces does not imply identical entries or conditioning.

For complex modes, `SelectedEigenvalues()` repeats the chosen representative
for its two real output functions, rather than necessarily returning both
conjugate signs. For repeated real eigenvalues, signs, scaling, and modal mixing
are also arbitrary. Output fields are **not automatically separated Mode I/II
or normalized to unit stress intensity**.

Optional modal changes are OFF by default. Existing JSON keys
`symmetry adapted modes` and `blockwise symmetry adapted modes` form modal
combinations; the latter requires explicit exponent ranges. Physical decoupling
requires the relevant reflection symmetry, not present for arbitrary rotated
orthotropic or interface cracks. `energy orthonormal modes` and
`energy normalize adapted modes` are optional experimental settings, not
requirements. The three primary transformation options are mutually exclusive;
consult the [implementation](Pre/TPZGFemEnrichment.cpp) and tests before using them.

## Conventions and limitations

- Do not evaluate derivatives at the singular tip `(0,0)`. The example rejects it and restricts points to the local disk.
- Near crack faces, approach from the desired side with a small signed offset and `EAny`. The internal `EUpper`/`ELower` families select angular-cut endpoints, not material IDs or general y half-planes. For `[-pi,pi]`, forced values use the first (`-pi`) and last (`+pi`) endpoint, respectively.
- Use a fresh object per configuration: repeated loading is not reliable yet. Avoid casual copying; mesh ownership/lifecycle hardening remains future work.
- Core input validation is incomplete; invalid cards may trigger NeoPZ diagnostics rather than recoverable exceptions. The examples add basic guards, not a complete schema validator.
- The documented evaluator is for 2D elasticity. An existing scalar setup path does not make scalar enrichment evaluation supported by this API.
- These tutorials do not implement a global GFEM solve, DSGFEM, condition-number estimation, or crack propagation.

## Troubleshooting

| Symptom | Check |
|---|---|
| Cannot find JSON | Install its development package; set `nlohmann_json_DIR` or `CMAKE_PREFIX_PATH` |
| Enrichment feature is OFF | Reconfigure, rebuild, and install NeoPZ with `USING_SBFEM_ENRICHMENT=ON` |
| Consumer cannot find NeoPZ | Set `NeoPZ_DIR` to installed `lib/cmake/neopz`; capitalization matters |
| JSON version mismatch | Match the version used to compile NeoPZ or rebuild it; JSON types occur in compiled public signatures |
| Undefined symbols/header mismatch | Do not mix installations; rebuild consumers after updating NeoPZ |
| Example cannot open JSON | Paths are relative to the current working directory, not the executable |
| No selected modes | Check exponent ranges, geometry, and material setup |
| Different signs/mixed fields across platforms | Compare modal spaces and derivatives, not individual eigenvector signs/order |

Keep the installation at its configured prefix, or reinstall at a new prefix
rather than copying it. For reproducibility record the commit, compiler, JSON
and BLAS/LAPACK versions, build options, and input JSON files.

## Docker image

We plan to provide a Docker image in the near future with NeoPZ, dependencies,
and the enrichment examples preconfigured. **An image is not published yet.**
Its registry address and instructions will be added here when available;
source installation is the current route.

## Citation

For the SBFEM-generated GFEM enrichments developed in this branch, please cite:

N. Shauer, C. Song, C. A. Duarte, and P. R. B. Devloo,
*A Well-Conditioned and Optimally Convergent SBFEM-Enriched GFEM for
Two-Dimensional Fracture in Orthotropic Media and Interfaces*.
Manuscript submitted for publication.

For NeoPZ itself, see P. R. B. Devloo, *PZ: An object oriented environment for
scientific programming*, Computer Methods in Applied Mechanics and Engineering
150 (1997), 133–153.
[DOI: 10.1016/S0045-7825(97)00097-2](https://doi.org/10.1016/S0045-7825(97)00097-2).
