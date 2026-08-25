---
type: flow
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - elasticity
  - hybridization
---

# Flow: iter_elast — hybrid H1 elasticity benchmark (mandated slice)

Driver: `divfreebubbles/targets/iter_elast.cpp` (read in full [repo]). Binary: `build/targets/iter_elast` (built Jun 30; runtime = installed NeoPZ @ `852a5116c` per install `pz_config.h` [run]).

## Shallow trace (Phase 2 answers)
1. **Problem**: 2D linear elasticity on unit square, homogeneous analytic case (`TElasticity2DAnalytic`, E=1, ν=0), discretized with *hybridized-squared H1* spaces; purpose = benchmark iterative Schur-complement solving vs sparse direct (time/memory sweep over 50²…400² quad meshes at fixed p) (:166-337).
2. **Start**: `main` loops `idivs={50,…,400}` (:189-192).
3. **Mesh**: `CreateGeoMesh<pzshape::TPZShapeQuad>` → `TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, {0,0,0},{1,1,1}, matIds, nDivs, EQuadrilateral, createBoundEls)` (:204-212, 417-466) → [[TPZGeoMesh]], [[mesh-io-generators]].
4. **Geometric elements**: structured quads + boundary 1D elements, all `matIds` boundary→`EBoundary`, volume→`EDomain`.
5. **Computational mesh**: `TPZH1HybridApproxCreator hdivCreator(gmesh)` (app-side, [[divfree-support-lib]]) → `CreateApproximationSpace()` returns `TPZMultiphysicsCompMesh` (:218-229) → [[TPZCompMesh]].
6. **Space selection**: `SetProbType(EElastic)`, `SetDefaultOrder(iorder)`, `SetExtraInternalOrder(2)`, `SetShouldCondense(false)`, `SetHybridType(EStandardSquared)` (:220-226) → [[approx-space-creators]], [[hybridization]]. Post-creation: `ComputeOrthogonalizingRestraints(*cmesh, geltogel, HybridData())`, `HybridizeLowOrderFluxes`, `GroupAndCondenseElements` (:231-233) → [[static-condensation]]. *(These three are app-side extensions over the develop-delta base `TPZH1ApproxCreator` — deep semantics = Phase 4 target.)*
7. **Materials**: `TPZHybridElasticity2D` (develop-delta file) on `EDomain` with forcing+exact from `TElasticity2DAnalytic`; single Dirichlet BC (type 0) on all boundaries via `CreateBC` + `SetForcingFunctionBC` (:470-494) → [[material-system]].
8. **Assembly**: iterative path delegates entirely to `TPZMatRedSolver<STATE>(an, EDarcyH1Hybrid).Solve(...)` (:287-289); direct path `TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>` (MKL) or `TPZSSpStructMatrixMumps` (32 threads) + `an.Assemble()` (:293-321) → [[structural-matrices]], [[assembly]].
9. **Solve**: iterative = Schur/matrix-reduction inside `TPZMatRedSolver` ([[divfree-support-lib]]); direct = `TPZStepSolver::SetDirect(ELDLt)` + `an.Solve()` (:315-332) → [[matrix-and-solvers]], [[TPZAnalysis]].
10. **Results/errors**: only equation counts + wall-times + peak RSS, appended to `results_Elastic2D_memory_time.txt`; `PostProcessError` block **commented out** (:391-410).
11. **VTK**: `TPZVTKGenerator` block **commented out** (:361-383) → no visualization leg in this slice.
12. **Central classes**: `TPZH1HybridApproxCreator` → `TPZH1ApproxCreator` (delta), `TPZHybridElasticity2D` (delta), `TPZMultiphysicsCompMesh`, `TPZLinearAnalysis` (`RenumType::ENone`!), `TPZMatRedSolver`, `TPZSSpStructMatrix(Mumps)`, `TPZStructMatrixOR`.
13. **Partially RESOLVED (Phase 4)**: `EStandardSquared` = double hybridization, fully traced — see [[hybridization]] (second interface/Lagrange layer; only 2nd-level skeleton global; atomic meshes = skeleton-flux HDivStandard + broken-H1 volume(+wraps+2nd-level primal skeleton); elastic multipliers {−1,−1,−1,1} with right interface reset to +1). Space construction matches [[avancini-2025-double-hybrid-elasticity]] structurally (H1-primal variant). Solver-mode question RESOLVED: `EDarcyH1Hybrid` vs `EElasticityH1Hybrid` share the sign-flip branch but select different preconditioner block sizes (`ord` vs `2(ord+1)−3`) — iter_elast's choice is a consequential mislabel for the benchmark ([[finding-matred-solver-mode-mislabel]]). Reduction anatomy fully traced (K00 = Lagrange-level-1 connects, Cholesky via Pardiso/MUMPS; matrix-free Schur CG with block-diagonal ELU preconditioner, 500 iter cap / 1e-10) — see [[matrix-and-solvers]] and `ALGORITHM_NOTES.md` §5. Still open: "orthogonalizing restraints" semantics (app-side `ComputeOrthogonalizingRestraints` — expert/maintainer question).

## Runtime trace [run @ 852a5116c(+), 2026-07-02, pOrder=1 iterative sweep from scratch cwd]
- Sweep completed exit 0; per-idiv console shows the reduction anatomy: e.g. idiv=50: `Number of equations = 168200` (full) → `condensed = 19600`; `NUMBER OF EQUATIONS: Full problem = 19600, High Order Flux = 4900, Linear Flux = 14700`; `Time Assembling SparseMatRed 78 ms`; `Assembling block diagonal 1 ms`; `Decomposing K00... 39 ms`; CG iteration log with residuals ~0.27→9e-11.
- **CG iterations = 19 at idiv=50 AND idiv=400** (residual ~9e-11 both) — mesh-size-independent iteration count; contraction ≈ 0.3/iter. Strong evidence the K00-block reduction acts as an (apparently) spectrally robust preconditioner for this problem family. (Deep mechanics: pending `TPZMatRedSolver` trace.)
- Output-table semantics (matches `results_Elastic2D_memory_time.txt` row `iterative 1 50 78 75 296304`): t1 = SparseMatRed assembly ms; t2 ≈ K00 decomposition + CG ms; third value = `getPeakMemoryMB()`.
- **Memory-unit bug (app-side)**: `getPeakMemoryMB` divides `ru_maxrss` by 1024 with a comment "ru_maxrss is in KB on Linux" (iter_elast.cpp:44-50) — on macOS `ru_maxrss` is in **bytes**, so values printed as "MB" are actually KiB on Mac (observed "Memory usage: 296304 MB" ≈ 289 true MB at 50²; 1.01205e7 "MB" ≈ 9.65 GB at 400²). Cross-platform 1024× distortion in published benchmark tables if mixed. → [[finding-rusage-memory-units]].
- PZ_LOG active in installed build: run creates a `LOG/` dir (46 log files) in cwd; config read from `neopz_install/pz/include/Util/log4cxx.cfg`.
- `results_Elastic2D.txt` rows are written without newline separators (cosmetic app bug).

## Quirks noted (app-side, shallow)
- File-scope `std::ofstream` globals open `results_*.txt` in CWD at static-init (two truncate; :53-55) — **never run from `build/targets/`**.
- The `exactSol` lambda fed to `an.SetExact` is a 3D scalar sinh-harmonic (zero at z=0), unrelated to the elasticity exact solution; harmless while error block is commented, but a trap if re-enabled (:88-154 vs :275) — [hypothesis: leftover from a Darcy variant of the benchmark].
- `hdivfamily` variable declared but unused by the H1 creator (:175).

Related: [[hybridization]] · [[static-condensation]] · [[approx-space-creators]] · [[divfree-support-lib]] · [[flow-dupl-connects]]
