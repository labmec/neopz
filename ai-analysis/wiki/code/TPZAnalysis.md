---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - analysis
  - solver
---

# TPZAnalysis — the solve orchestrator

## Responsibility
"Implements the sequence of actions to perform a finite element analysis" (header @brief [repo]): owns the [[TPZCompMesh]], equation renumbering, a [[structural-matrices|structural matrix]], a solver ([[matrix-and-solvers]]), the solution vector, exact-solution hooks, and post-processing ([[post-processing-vtk]], [[error-estimation-convergence]]).

## Key facts (verified [repo], Analysis/TPZAnalysis.h)
- Members: `TPZGeoMesh* fGeoMesh`, `TPZCompMesh* fCompMesh`, `TPZGraphMesh* fGraphMesh[3]` (one per dim), `TPZSolutionMatrix fSolution`, `TPZSolver* fSolver`, scalar/vector/tensor post-process name tables (lines 61-80).
- Renumbering options: `RenumType {ENone, EDefault(Metis-or-Sloan), ESloan, ECutHillMcKee, ECutHillMcKeeFast, EMetis}` (lines 47-54) — bandwidth/fill reduction before equation numbering.
- Built-in preconditioner factory: `Precond::{Jacobi, BlockJacobi, Element, NodeCentered}` (lines 30-45).
- Subclasses: `TPZLinearAnalysis` (linear static; the one used by divfreebubbles), `TPZEigenAnalysis` (+quadratic), `pznonlinanalysis`, `pztransientanalysis`, `pzmganalysis` (multigrid), substructure/frontal variants [agent paths].

## Assemble/Solve mechanics (verified [repo Analysis/TPZLinearAnalysis.cpp:35-180])
- `Assemble()` dispatches real/complex → `AssembleT<TVar>`: if no struct matrix set, defaults to `TPZSpStructMatrix` (MKL) else nonsym skyline, with console notice; if no solver, defaults to LU; **matrix reuse**: if solver already holds a right-sized matrix it's zeroed and re-assembled in place, else `fStructMatrix->CreateAssemble(fRhs)` builds it (:57-90). RHS sized by `ComputeNumberofLoadCases()`.
- `Solve()` → `SolveT<TVar>`: guards rhs size, respects `NReducedEquations()` (equation-filter path), computes residual norm, delegates to the `TPZMatrixSolver` (:128-180).

## Canonical use (observed in `divfreebubbles/targets/iter_elast.cpp:274-337` [repo])
`TPZLinearAnalysis an(cmesh, RenumType::ENone); an.SetExact(...); an.SetStructuralMatrix(matskl); an.SetSolver(step); an.Assemble(); an.Solve();` then error via `an.PostProcessError(...)` and VTK via `TPZVTKGenerator` (or legacy `DefineGraphMesh/PostProcess`).

## Related
[[TPZCompMesh]] · [[structural-matrices]] · [[matrix-and-solvers]] · [[assembly]] · [[post-processing-vtk]] · [[error-estimation-convergence]]

## Open questions
- `TPZAnalysis` is `TPZSavable` — is a full analysis actually serializable in practice? (persistence coverage is thin [agent]) → Phase 7.
- Ownership of `fSolver`/struct-matrix (raw pointers + clones?) → Phase 5.
