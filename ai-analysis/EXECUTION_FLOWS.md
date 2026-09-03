# NeoPZ Execution Flows

**Version 1 (Phase 2 — shallow discovery).** Deepened in Phase 4. Analyzed commit `develop @ 6ffd38b12`; app repo `../divfreebubbles @ 3DKernelHdiv`; runtime binaries link installed NeoPZ stamped `852a5116c` (`neopz_install/pz/include/Common/pz_config.h` [repo]).
Full per-slice pages: `wiki/flows/`.
**Session-2 scope note:** the five slices below are all divfreebubbles/unit-test flows (HDiv/hybridization-weighted). Statically traced flow *sketches* from five further downstream apps — adaptive estimator loops, GFEM enrichment builds, Uzawa outer iterations, 3D/2D/1D coupled builds, MHM-controller elasticity — now live in `wiki/apps/` (one page per app + [[apps-overview]]); they were not executed, so they carry app-repo evidence class, not [run].

## The canonical NeoPZ pipeline (as observed across all five slices)

```
        geometry                 spaces                     physics              system                solve                output
TPZGeoMesh  ──────►  TPZCompMesh / TPZMultiphysicsCompMesh  ──────►  TPZMaterial+BCs  ──►  TPZLinearAnalysis  ──►  StepSolver/    ──►  VTK writers /
(GmshReader,          (ApproxSpace factory  or                        (per matid,          + TPZStructMatrix       Schur(MatRed)/      PostProcessError
 GenGrid, Tools)       TPZ*ApproxCreator [+hybridize,                  weak form)           (storage+threads)      Pardiso/MUMPS       (SetExact)
                       +condense, +substructure])
```

Two space-construction idioms coexist:
- **Manual** (older, `flow-dfreebubbles-1el`): per-field cmeshes with `TPZNullMaterial` placeholders + `ApproxSpace().SetAllCreateFunctions*` + `AutoBuild`, combined via `TPZMultiphysicsCompMesh`; even fully manual per-element `new TPZCompElKernelHDiv<...>` with neighbor-walking for wrap elements.
- **Creator-driven** (current, all other slices): `TPZHDivApproxCreator`/`TPZH1ApproxCreator` (+MHM/app-side derivatives) encapsulating multi-mesh construction, hybridization (`ENone/EStandard/EStandardSquared/ESemi`), rigid-body enrichment and condensation.

## Slice summaries

| Slice | Problem | Space idiom | Solve | Output legs | Page |
|---|---|---|---|---|---|
| iter_elast (mandated) | 2D hybrid-squared H1 elasticity, analytic homogeneous | app-side `TPZH1HybridApproxCreator` + orthogonalizing restraints + low-order-flux hybridization + group/condense | `TPZMatRedSolver(EDarcyH1Hybrid)` Schur vs direct MKL/MUMPS LDLt, 32 thr | timing+memory tables only (error/VTK commented) | [[flow-iter-elast]] |
| dupl_connects2 | 2D/3D mixed Darcy, `EHDivConstant` | lib `TPZHDivApproxCreator`, `ESemi`, condense=on | `TPZMatRedSolver(EDarcyHDiv)` vs direct | timing tables (error/VTK commented) | [[flow-dupl-connects]] |
| MHM_HDivConstant | Darcy on polygonal partition (quadtree import) | app-side MHM geo+approx creators, rigid-body on, substructures + condense | direct LDLt (ESemi/`EMHMSparse` branch = dead & **drifted**) | cmesh-as-VTK + txt dumps (always-on debug) | [[flow-mhm-hdivconstant]] |
| dFreeBubbles1el | Darcy on 1-element gmsh mesh | manual flux+pressure meshes, factory `SetAllCreateFunctionsHDiv` | direct LDLt | **error computation + multiphysics VTK active** | [[flow-dfreebubbles-1el]] |
| Unit tests | linear/constant representation grid; De Rham rank/kernel | lib creators incl. lib MHM creator; basis-level Gram matrices | direct / SVD (LAPACK) | Catch2 assertions (+VTK smoke) | [[flow-unit-test-hdiv-creator]] |

## Cross-slice observations (feeding Phases 4–7)

1. **The 5-file develop delta sits on the hot path of the mandated slice**: `TPZH1ApproxCreator` (base of the app creator) and `TPZHybridElasticity2D` (material). Runtime evidence must carry the `[run @ 852a5116c]` label (install stamp; possibly stale vs dylib rebuild — nuance logged).
2. **Solver-mode naming mismatch**: iter_elast passes `EDarcyH1Hybrid` for an elasticity problem though `EElasticityH1Hybrid` exists (`divfree/TPZMatRedSolver.h:15`) — either the reduction structure is problem-agnostic (naming debt) or a real mis-selection → Phase 4 must read `TPZMatRedSolver::Solve`.
3. **App-repo source drift**: `MHM_HDivConstant` references removed enumerator `EMHMSparse`; current tree likely doesn't rebuild that target (binary is Mar 27). OQ6.
4. **Benchmark-first hygiene**: error/VTK legs are commented out in the two benchmark drivers; analytic lambdas left inconsistent (3D solution in 2D runs; forcing not matching exact solution; dead assignments inside lambdas). These are app-side, not library, issues — but they mean **the benchmarks currently validate performance, not correctness** (correctness legs live in dFreeBubbles1el + unit tests).
5. **Renumbering disabled** (`RenumType::ENone`) in both benchmarks — deliberate for the Schur solver? Bandwidth effects on direct path? → Phase 6 question.
6. **Same-name classes in lib and app** (`TPZMHMHDivApproxCreator` in `Pre/` and in `divfree/`) — which one a target gets depends on include paths; migration-in-progress pattern → Phase 5 risk note.
7. **Outputs land in CWD** across drivers (results txt, mphysics2.txt, cmesh_multi.vtk, gmesh.vtk) — any run of these binaries must use a scratch CWD (engagement rule already in place).

## What each output means (v1 + Phase 4 run evidence)

- `results_*_memory_time.txt`: `<mode> <pOrder> <idiv> <t1> <t2> <mem>` — confirmed by a fresh run [run @ 852a5116c(+), 2026-07-02]: for `iterative` rows t1 = SparseMatRed assembly (stdout "Time Assembling SparseMatRed 78 ms" ↔ t1=78), t2 ≈ K00 decomposition + CG; `mem` = `ru_maxrss/1024` — **KiB on macOS mislabeled as MB** (1024× overstatement; correct on Linux) → [[finding-rusage-memory-units]].
- Run anatomy (iter_elast, p=1, 50²): 168,200 full eqs → 19,600 condensed → split "High Order Flux" 4,900 + "Linear Flux" 14,700; K00 factorized; **CG converges in 19 iterations at 50² and at 400²** (mesh-independent count, contraction ≈0.3/iter) — strong empirical robustness signal for the reduction; mechanism documented after the MatRedSolver trace.
- `postprocessHdiv.txt` (1el): `TPZAnalysis::PostProcessError`-style norms vs the (effectively linear) exact solution — indices/norms documented in Phase 4.
- `cmesh_multi.vtk` / `gmesh.vtk`: geometric visualization (mesh + matids/partition), not solution fields.
