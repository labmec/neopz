# NeoPZ Codebase Atlas

**Analyzed commit:** `develop` @ `6ffd38b12` (2026-06-12) in `labmec/neopz` (local clone `NeoPZ_divfree`).
**Working-tree caveat:** checkout is `SemiHybridElasticity` @ `4de234fae` = develop + 3 commits touching exactly 5 files
(`Material/Elasticity/TPZHybridElasticity2D.cpp`, `Matrix/TPZSYSMPMatrix.h`, `Matrix/pzmatrix.h`, `Pre/TPZH1ApproxCreator.{h,cpp}`);
claims about those files are cross-checked against `git show develop:<path>`.
**Repo default branch** is `main`; `develop` is the integration branch and the review canon (user instruction).
**Evidence markers:** `[repo]` verified first-hand at the pin; `[agent]` located by a read-only explorer and cited but not yet independently re-verified; `[run]` runtime observation of installed builds. Details and page-level citations live in `ai-analysis/wiki/`.

---

## 1. What this repository is

NeoPZ is a single-library C++17 finite element environment (`project(PZ)`, one `add_library(pz)` target — CMakeLists.txt:8,52 [repo])
developed by LabMEC/Unicamp since the 1990s (Devloo, CMAME 1997). It provides, per README.md:14-21 [repo]: discontinuous, H1-, H(div)- and
H(curl)-conforming approximation spaces, multiphysics, hp-adaptivity, hanging nodes, runtime-defined refinement patterns, exact curved
geometry, and forward automatic differentiation. It is a library only — applications live in downstream repos
(here: `../divfreebubbles`, see §7).

## 2. Module map

Single library; top-level directories are source groups appended to `pz` in dependency-ish order (CMakeLists.txt:320-343 [repo]).
File counts measured on the working tree [agent, spot-checked].

| Module | .h/.cpp | Responsibility (one line) | Key types |
|---|---|---|---|
| `Util/` | 36/24 | Foundation containers & helpers: vectors, stacks, chunk vectors, ref-counted pointer, logging | `TPZVec`, `TPZManVector`, `TPZStack`, `TPZChunkVector`, [[TPZAutoPointer]], `pzlog` |
| `Common/` | 53/7 | Numeric type config (REAL/STATE), error macros, element-type enum, thread pool | `pzreal.h`, `DebugStop`, `MElementType`, `TPZThreadPool` |
| `Save/` | 13/11 | Object persistence with ClassId registry + versioned chunk translators | `TPZSavable`, `TPZStream`, `TPZPersistenceManager` → [[persistence]] |
| `Integral/` | 11/9 | Quadrature rules per topology (long double precision internally) | `TPZIntPoints`, `tpzgaussrule` → [[quadrature]] |
| `Solvers/` | 27/16 | Linear & eigen solver wrappers incl. Pardiso/MUMPS bindings | `TPZMatrixSolver`, `TPZStepSolver`, `TPZPardisoSolver` → [[matrix-and-solvers]] |
| `Matrix/` | 32/29 | Matrix storage zoo: dense, banded, skyline, Yale sparse (sym/nonsym) | `TPZFMatrix`, `TPZSkylMatrix`, `TPZFYsmpMatrix` → [[matrix-and-solvers]] |
| `Topology/` | 11/10 | Master-element combinatorics: sides, permutations, transforms | `pztopology::TPZTriangle` … → [[topology-module]] |
| `Geom/` | 15/11 | Reference→physical geometric maps per topology + blend maps | `pzgeom::TPZGeoQuad`, `pznoderep`, `tpzgeoblend` → [[geometric-mappings]] |
| `SpecialMaps/` | 20/19 | Exact curved maps (arcs, spheres, tori, cylinder, NACA, quadratic els) | `TPZArc3D`, `TPZCylinderMap`, `TPZQuadSphere` |
| `Shape/` | 26/22 | Shape functions as **static, non-virtual** per-topology classes; H1/HDiv/HCurl families | `pzshape::TPZShapeQuad`, `TPZShapeHDiv*`, `TPZShapeHCurl*`, `TPZShapeData` → [[shape-functions]] |
| `Refine/` | 12/11 (+71 `.rpt`) | h-refinement via runtime refinement patterns (data-driven) | `TPZRefPattern`, `TPZRefPatternDataBase` → [[refinement-hanging-nodes]] |
| `External/` | 10/17 | Vendored renumbering/partitioning: Sloan, (R)CM, METIS wrapper, Boost graph | `TPZRenumbering`, `TPZCutHillMcKee`, `TPZSloanRenumbering` |
| `Material/` | 199/178 | Weak forms & constitutive models; variadic mixin base + per-physics dirs + big legacy layer | [[material-system]]: `TPZMaterial`, `TPZMatBase<TVar,Ifaces...>`, `TPZBndCond`; `needrefactor/` (19+108 files) |
| `Mesh/` | 65/61 | Geometric & computational meshes; element hierarchies for all space families | [[TPZGeoMesh]], [[TPZCompMesh]], `TPZCompEl`, `TPZInterpolatedElement`, [[TPZCompElHDiv]], `TPZCompElHCurl`, `TPZMultiphysicsCompMesh`, `TPZSubCompMesh`, `TPZCondensedCompEl` |
| `Analysis/` | 13/13 | Drives the solve sequence: renumber → assemble → solve → post-process | [[TPZAnalysis]], `TPZLinearAnalysis`, `TPZEigenAnalysis` |
| `Post/` | 29/26 | Post-processing/output: legacy graph meshes (DX/MV/V3D/VTK) + modern VTK generator | [[post-processing-vtk]]: `TPZVTKGenerator`, `TPZGraphMesh`, `TPZVTKGeoMesh` |
| `Frontal/` | 8/9 | Frontal (out-of-core style) solver machinery | `TPZFront`, `TPZFrontMatrix` |
| `StrMatrix/` | 24/24 | "Structural matrices": assembly strategies binding mesh↔storage↔parallel scheme | [[structural-matrices]]: `TPZStructMatrix(T)`, `TPZSkylineStructMatrix`, `TPZSSpStructMatrix`, `TPZEquationFilter` |
| `Pre/` | 34/31 | Pre-processing: mesh readers/generators, analytic solutions, **approximation-space creators**, hybridization, MHM controllers | [[mesh-io-generators]], [[approx-space-creators]]: `TPZGmshReader`, `TPZGenGrid2D/3D`, `TPZAnalyticSolution`, `TPZCreateApproximationSpace`, `TPZApproxCreator` → `TPZHDivApproxCreator`/`TPZH1ApproxCreator`, `TPZHybridizeHDiv`, `TPZMHM*MeshControl` |
| `SubStruct/` | 13/12 | Substructuring / BDDC (Dohrmann) domain decomposition | `TPZDohrStructMatrix`, `tpzdohrsubstruct` |
| `Random/` | 5/1 | Random number generators (uniform/normal/constrained) | `TPZRandom` |
| `Optimization/` | 1/0 | Stub: stochastic search | `TPZStochasticSearch` |
| `Exception/` | 3/3 | Typed exceptions | `TPZConvergenceException` |
| `PerfUtil/`, `PerfTests/` | — | Perf measurement toolkit + standalone benchmarks (stale; "in need of a revision" per own README [agent]) | `run_stats_table`, `SubStruct/substruct.cpp` |
| `Publications/` | 3 pairs | Paper-companion code: `hdiv2dpaper201504`, `hdiv3dpaper201504`, `hdivCurvedJCompAppMath` [repo] | — |
| `UnitTest_PZ/` | 33 suites | Catch2 v3.3.2 test suites (§6) | — |

Quirks worth knowing [repo]: 4 double-extension `.h.h` files carry template bodies (`Geom/pznoderep.h.h`, `Mesh/TPZGeoElement.h.h`,
`Mesh/pzgeoelrefless.h.h`, `Mesh/tpzgeoelrefpattern.h.h`); two naming eras coexist (`pz*.h` older, `TPZ*.h` newer); `Refine/RefPatterns/`
holds 71 `.rpt` refinement-pattern *data* files loaded at runtime.

## 3. Dependency layering (first pass)

CMake `add_subdirectory` order (CMakeLists.txt:320-343 [repo]) approximates the layering; all code links into one `pz` target, so
"dependencies" are include-level, not link-level:

```
foundation:      Util → Common → Save
numerics:        Integral, Solvers, Matrix
reference layer: Topology → Geom → SpecialMaps → Shape → Refine
discretization:  External, Material → Mesh
orchestration:   Analysis, Post, Frontal, StrMatrix, Pre, SubStruct
extras:          Random, Optimization, Exception
```

Observed include-level couplings [repo]: `Mesh/pzcmesh.h` includes `pzcreateapproxspace.h` (Pre) — i.e. **Mesh ↔ Pre are mutually
entangled** (the computational mesh owns a `TPZCreateApproximationSpace` member); `Analysis/TPZAnalysis.h` includes StrMatrix + Solvers +
External (renumbering); `TPZGeoMesh` and `TPZCompMesh` hold mutual `fReference` pointers (pzgmesh.h:55, pzcmesh.h:49 [repo]). A true
include-graph pass is deferred to Phase 5.

## 4. Build system

[repo, verified] CMake ≥ 3.14 (`CMakeLists.txt:3`; README says 3.13+ — minor doc mismatch), C++17 enforced (`:13-14`), single shared lib
(static on Windows). Type system via `cmake/StandardPZSettings.cmake`: `REAL_TYPE` (geometry scalar) and `STATE_TYPE` (FE scalar) each
float/double/long double, plus `REAL_TYPE=pzfpcounter` (op-counting instrumented scalar); complex state = `CSTATE`. Generated
`Common/pz_config.h` embeds source dir, refpattern dir, git branch/revision/date.
24 `option()` flags [repo count]: only Threads is required; optional `USING_MKL` (forces LAPACK; Pardiso), `USING_MUMPS`, `USING_LAPACK`
(Accelerate on Apple), `USING_METIS`, `USING_TBB`, `USING_OMP`, `USING_LOG4CXX` (→`PZ_LOG`), `USING_BOOST`, `USING_EIGEN`,
`USING_UMFPACK`, `USING_BLAZE`, perf instrumentation (`LIKWID/PAPI/LIBNUMA`), and build toggles (`BUILD_UNITTESTING`, `BUILD_PERF_TESTS`,
`BUILD_PUBLICATIONS`, `BUILD_PROJECTS`, `BUILD_PLASTICITY_MATERIALS`, `BUILD_DOCS`/`BUILD_SPHINX_DOCS`) [agent, options spot-verified].
Install exports a proper CMake package: `find_package(NeoPZ)` → `NeoPZ::pz` (install layout `<prefix>/pz/{include,lib}` +
`<prefix>/lib/cmake/neopz`) [agent; confirmed in practice by divfreebubbles' cache `NeoPZ_DIR=.../NeoPZ_install/lib/cmake/neopz` [repo]].

## 5. Entry points a user actually touches

Typical downstream flow (confirmed in divfreebubbles targets [repo]):
1. **Geometry**: `TPZGeoMesh` from `TPZGmshReader` (in-tree .msh parser — gmsh is *not* linked), `TPZGenGrid2D/3D`, or
   `TPZGeoMeshTools::CreateGeoMeshOnGrid`.
2. **Spaces**: either manual per-space `TPZCompMesh` construction + `TPZMultiphysicsCompMesh`, or the modern
   `TPZHDivApproxCreator`/`TPZH1ApproxCreator` layer (Pre/) with `ProblemType{EDarcy,EElastic,EStokes}`,
   `HybridizationType{ENone,EStandard,EStandardSquared,ESemi}`, optional rigid-body spaces + condensation
   (Pre/TPZApproxCreator.h:15-16,42-68 [repo]).
3. **Physics**: a `TPZMatBase<STATE, Interfaces...>` material per material-id + `CreateBC` boundary conditions ([[material-system]]).
4. **Solve**: `TPZLinearAnalysis` + a `TPZStructMatrix` flavor + `TPZStepSolver` (direct ELDLt/ELU/Cholesky or CG/GMRES) — or external
   Pardiso/MUMPS via sparse struct-matrices.
5. **Output**: `TPZVTKGenerator` (modern, NGSolve-derived — Post/TPZVTKGenerator.h:1-6 [repo]) or legacy `TPZGraphMesh` family;
   errors via `TPZAnalysis::PostProcessError` against `SetExact`.

## 6. Tests, CI, docs (condensed; full review in Phase 7)

[agent, key items spot-verified] Catch2 v3.3.2 (FetchContent), 33 suites under `UnitTest_PZ/`, ~180+ TEST_CASEs. Notably *mathematical*
suites: `TestDeRham` (De Rham complex exactness via rank/kernel with SVD), `TestMesh/TestHDiv.cpp` (De Rham checks incl. under face/node
permutations, side-shape continuity, order checks), `TestTopology` (constant div/curl reproduction per topology), `TestHCurl`,
`TestHDivApproxSpaceCreator` (parametrized: HDiv family × Darcy/Elastic × mesh type × hybridization × condensation),
`TestSBFem` (convergence), `TestSolverComparison` (MUMPS vs Pardiso), `TestMultithreading` (parallel==serial assembly).
Known-failing tests marked `[!shouldfail]` (SVD, some skyline ops). CTest granularity = one test per suite executable
(`catch_discover_tests` disabled). CI = 5 GitHub Actions workflows; macOS job is the always-on test gate; Linux + MKL jobs run only when a
prebuilt `neopz-deps` image exists; no Windows CI, no coverage tooling. Docs = Doxygen + Sphinx/Breathe → labmec.github.io/neopz,
published from `main` only. No LICENSE/CONTRIBUTING files in tree [agent].

## 7. Downstream application landscape (Session 2, 2026-07-06)

The library's user base is visible in `~/GitHub`: ~20 research repos embed their own NeoPZ copy (all near-develop, few local changes). The five most recently active were surveyed read-only (`wiki/apps/`, spot-verified):

| App (active) | Domain | NeoPZ surface it stresses |
|---|---|---|
| Iterative-Saddle_Point (2026-03) | Uzawa/augmented-Lagrangian mixed Darcy & Stokes | Matrix/Pardiso layer as user API; manual Stokes hybridization + condensation |
| GFEM (2025-12) | fracture GFEM enrichment | `TPZCompElH1` subclassing (`ComputeShape` override); SBFem eigenmodes; custom StrMatrix/MatRed |
| ErrorEstimation (2025-11) | a-posteriori estimators + hp-adaptivity | mesh cloning/reconstruction, patch solves, closed adaptive loops, quarter-point/SBFem/NACA geometry |
| wann (2025-09) | 3D/2D/1D reservoir–wellbore Darcy → ANN data | connect surgery (`AddDependency`), cylinder+blend maps, directional refinement |
| MixedElasticity (2025-05) | Hellinger–Reissner tensor elasticity | 3–7-field multiphysics, weak/strong symmetry, MHM controllers + submeshes, frontal solver |

Cross-cutting: manual space construction is first-class API everywhere; the connect layer is a public research surface; MatRed-style reduction reimplemented independently 3×; refinement patterns used for adaptivity, grading, *and* space construction; ~20 downstream material subclasses vs 1 element subclass; app→lib same-name-class migration is the growth mechanism (details: [[apps-overview]]).

## 7b. Companion application repo: `../divfreebubbles` (branch `3DKernelHdiv`; the Session-1 vehicle)

[repo/agent] Research app repo ("div-free bubbles" Laplace per README — stale; actual scope now spans mixed/hybrid Darcy and elasticity).
Own support lib `divfree/`: materials (`TPZMatDivFreeBubbles`, `TPZMixedDarcyH1`, …), creators (`TPZH1HybridApproxCreator` — derives from
NeoPZ's `TPZH1ApproxCreator`, one of the 5 delta files; `TPZMHMGeoMeshCreator`, `TPZMHMHDivApproxCreator`), custom H(div) elements with
duplicated connects (excluded from build), Schur-complement solvers (`TPZMatRedSolver` modes `EDefault/EDarcyHDiv/EDarcyH1Hybrid/EMHMSparse`,
`TPZSparseMatRed`). 8 built targets (`iter_elast`, `dupl_connects2`, `MHM_HDivConstant`, `dFreeBubbles1el`, `2frac`, `hpc4`,
`voronoi_mixed_elas`, `semiHybrid_elas`). `divfree_config.h` bakes absolute `MESHDIR`. Catch2 tests exist behind `BUILD_TESTS=OFF`.
Vertical slices selected for this assessment: see [[flow-iter-elast]] (mandated), [[flow-dupl-connects]], [[flow-mhm-hdivconstant]],
[[flow-dfreebubbles-1el]], [[flow-unit-test-hdiv-creator]] (pages created in Phase 2).

## 8. Concepts queued for reference research (Phase 3)

[[h1-space]], [[hdiv-space]], [[hcurl-space]], [[de-rham-complex]], [[mixed-methods]], [[hybridization]] (standard / squared / semi),
[[static-condensation]], [[mhm]], [[sbfem]], [[refinement-hanging-nodes]], [[geometric-mappings]], [[piola-transformations]],
[[quadrature]], [[assembly]], [[error-estimation-convergence]], [[vtk-output]].

## 9. Corrections & uncertainty ledger

- C1 [fixed]: explorer report placed `pzcreateapproxspace.h` in `Mesh/`; actually `Pre/pzcreateapproxspace.h` [repo]. (`Mesh/pzcmesh.h:18`
  includes it, which likely caused the confusion — and reveals the Mesh↔Pre coupling noted in §3.)
- C3 [fixed, Session 2]: Session-1 claims that `Material/needrefactor/` "still compiles into `pz`" were wrong — no `add_subdirectory(needrefactor)` exists and `libpz.dylib` contains none of its symbols [repo, verified by nm]. The Material row's 199h/178cpp counts include needrefactor *files on disk*, not compiled code. Residual risk = header/name shadowing + a few out-of-library targets including its headers.
- Open questions tracked in `wiki/log.md` (OQ1–OQ5): approx-space copy semantics, installed-build revision stamp, `[!shouldfail]` scope,
  `voronoi_mixed_elas` null-`gAnalytic` suspicion, missing LICENSE.
- All `[agent]` claims above are treated as *located, pending re-verification*; they get promoted to `[repo]` (with line citations) as the
  relevant phases touch them.
