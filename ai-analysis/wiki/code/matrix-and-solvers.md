---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - linear-algebra
  - solvers
---

# Matrix/ + Solvers/ — storage & linear solvers

## Responsibility
`Matrix/` is the storage zoo (all deriving from `TPZMatrix<TVar>` over `TPZBaseMatrix`); `Solvers/` wraps direct and iterative solution strategies consumed by [[TPZAnalysis]].

## Matrix storage [repo paths; hierarchy details to verify Phase 5/6]
- `Matrix/pzmatrix.h` — abstract `TPZMatrix<TVar>` (**in the 5-file develop delta** — the delta makes `MultiplyByScalar` virtual per commit messages; cross-check `git show develop:Matrix/pzmatrix.h` before citing internals).
- `Matrix/pzfmatrix.h` — `TPZFMatrix<TVar>` dense column-major workhorse (also the RHS/solution container).
- `Matrix/pzskylmat.h` (`TPZSkylMatrix` symmetric skyline + in-house Cholesky/LDLt), `pzskylnsymmat.h` (nonsym skyline), `pzbndmat.h`/`pzsbndmat.h` (banded), `pzsfulmat.h` (sym full), `pzblock.h` (block indexing), `pzblockdiag.h` (block diagonal).
- Sparse: `TPZYSMPMatrix.h` (Yale/CSR nonsym), `TPZSYSMPMatrix.h` (sym CSR; **delta file**), Pardiso-backed variants (`TPZYSMPPardiso.h`, `TPZSYSMPPardiso.h`), MUMPS variants (`TPZSYSMPMumps.h` etc.), `TPZEigenSparseMatrix.h` (Eigen/Accelerate bridge) [agent].
- `Matrix/TPZMatrixWindow.h` (windowed views), `TPZTensor.h` (plasticity tensors) [agent].

## Solvers [repo paths]
- `Solvers/TPZSolver.h` → `TPZMatrixSolver<TVar>` (holds the matrix via [[TPZAutoPointer]]).
- `Solvers/pzstepsolver.h` — `TPZStepSolver`: SetDirect(ELU/ECholesky/ELDLt) or iterative CG/GMRES/Jacobi/SSOR with optional preconditioner (another `TPZMatrixSolver`); used everywhere downstream.
- `Solvers/TPZPardisoSolver.h` (MKL), `TPZMumpsSolver.h` (MUMPS) [agent].
- Eigen stack (Session 2, verified [✓ `Solvers/EigenSolvers/` listing]): `TPZEigenSolver` base (targets, npairs, generalised-vs-standard) → `TPZLinearEigenSolver` (Ax=λx / Ax=λBx), `TPZLapackEigenSolver` (dense/banded LAPACK), `TPZKrylovEigenSolver(+Base)` (Arnoldi projection), `TPZQuadEigenSolver` (quadratic EVP via shift-invert Krylov), `TPZSpectralTransform` (shift / shift-and-invert), `TPZEigenSort`. Analysis drivers: `TPZEigenAnalysis` (A/B matrices ↔ `TPZMatGeneralisedEigenVal`), `TPZQuadEigenAnalysis` (K/L/M ↔ `TPZMatQuadraticEigenVal`) — STATE and CSTATE instantiated; primary consumers are the complex electromagnetics materials ([[material-system]]). Note [[sbfem]] bypasses this stack (direct `dgeev_`/blaze).
- Renumbering lives in `External/` (Sloan, Cuthill-McKee, METIS, Boost) selected via `RenumType` in [[TPZAnalysis]] (TPZAnalysis.h:48-54 [repo]).

## TPZMatRed (verified [repo Matrix/pzmatred.h:23-79])
2×2 block substructuring container `[K00 K01; K10 K11]`, side-matrix storage templated (`TPZFMatrix` or `TPZVerySparseMatrix`), holds a `TPZMatrixSolver` for K00, tracks `fK01IsComputed/fIsReduced` state, and is **rigid-body-mode aware** (`fMaxRigidBodyModes`, `fNumberRigidBodyModes`) — floating-subdomain support built into the reduction core (feeds [[static-condensation]] and [[mhm]]). Note: `CopyFrom` lacks the self-assignment guard (`if (from)` only, :66-79) that was added to `TPZSYsmpMatrix` post-pin — same latent-bug family as [[finding-hybridelasticity2d-missing-rhs-at-pin]] notes.

## Decomposition state machine
`TPZMatrix` carries a decomposition flag (`ENoDecompose/ELU/ECholesky/ELDLt`) so repeated `Solve` reuses factors [pattern known from usage; verify]. In-house factorizations coexist with LAPACK/BLAS replacements when `USING_LAPACK` (README.md:38 [repo]).

## Related
[[structural-matrices]] · [[TPZAnalysis]] · [[matrix-and-solvers]]-consumers: [[flow-iter-elast]] (`TPZMatRedSolver` app-side Schur), [[static-condensation]] (`TPZMatRed`, in Matrix/ [agent: `pzmatred.h`])

## Open questions
- `TPZMatRed` (library) vs divfreebubbles `TPZMatRedSolver`/`TPZSparseMatRed`: which reduction machinery is lib vs app? → Phase 4 (iter_elast slice). Established so far [repo]: `divfree/TPZMatRedSolver.h:15` enum `ProblemOrigin {EDarcyHDiv, EElasticityHDiv, EDarcyH1Hybrid, EElasticityH1Hybrid}` — no `EDefault`/`EMHMSparse` (older drivers still reference them → app-side drift, OQ6).
- Thread-safety of shared matrix objects across solver threads — Phase 5/6.
