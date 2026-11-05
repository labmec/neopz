# SBFEM crack-tip enrichment generator tests

`TPZGFemEnrichment` is an optional part of `NeoPZ::pz`. It builds a local circular
SBFEM discretization and evaluates selected crack-tip functions and Cartesian
derivatives. No separate GFEM repository, Eigen, or autodiff is required.

Configure, build, test, and install with:

```sh
cmake -S . -B build -DUSING_LAPACK=ON -DUSING_SBFEM_ENRICHMENT=ON -DBUILD_UNITTESTING=ON -DCMAKE_INSTALL_PREFIX=/path/to/neopz_install
cmake --build build --parallel
ctest --test-dir build -R 'Test_TPZGFEM|TestSBFem' --output-on-failure
cmake --install build
```

The feature is OFF by default. When enabled, it requires nlohmann JSON >= 3.2.0
(provide `nlohmann_json_DIR` or `CMAKE_PREFIX_PATH` if needed). The installed
NeoPZ package finds and exports this dependency transitively, requiring the same
JSON version used to build NeoPZ because JSON types occur in the compiled API.
Consumers include `TPZGFemEnrichment.h` and link only `NeoPZ::pz`:

```cmake
find_package(NeoPZ CONFIG REQUIRED)
if(NOT NeoPZ_SBFEM_ENRICHMENT)
    message(FATAL_ERROR "Install NeoPZ with USING_SBFEM_ENRICHMENT=ON")
endif()
target_link_libraries(my_program PRIVATE NeoPZ::pz)
```

## Tests transferred from GFEM

- `Test_TPZGFEMEnrichment`: construction and agreement between evaluation APIs.
- `Test_TPZGFEMMaterials`: isotropic limit, orthotropic and bimaterial materials.
- `Test_TPZGFEMRegression`: complex reconstruction, optional modal transforms,
  Cartesian finite-difference gradients, crack-face evaluation and rotated
  angular cuts, Williams exponent ranges, defaults, and option validation.

Space-equivalence tests fit a constant invertible modal transformation and check
it at independent points and on gradients; they do not require identical
eigensolver signs/phases. Gradient checks use a centered step of 1e-6 and relative
tolerance 2e-6; held-out modal-space checks use 1e-8. These are local tests, not
global GFEM convergence/conditioning tests. Existing NeoPZ Catch2 infrastructure
is reused; no separate runner or duplicate `test_library` is introduced.

## Compatibility scope of the initial transfer

The generator implementation and defaults are preserved from GFEM commit
`6aacc36`; this transfer does not change existing NeoPZ SBFEM algorithms. The
validated application is 2D elasticity with double REAL/STATE and LAPACK.
Use a fresh generator per configuration and keep it alive while using its
callback. Existing limitations concerning repeated loading, raw mesh ownership,
malformed geometry input, and scalar evaluation were not redesigned in this
port. Such robustness changes should have dedicated tests and a separate review.
Always rebuild downstream applications against matching installed headers and
library after updating NeoPZ.
