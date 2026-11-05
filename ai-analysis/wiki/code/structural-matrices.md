---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - assembly
  - parallel
---

# StrMatrix/ — structural matrices (assembly strategies)

## Responsibility
A `TPZStructMatrix` binds three choices: (1) which global matrix *storage* to create, (2) how to *assemble* it from element matrices (serial / threaded variants), (3) which equations enter (equation filter). [[TPZAnalysis]] delegates `Assemble()` to it → [[assembly]].

## Key files
- `StrMatrix/TPZStructMatrix.h` [repo]: type-agnostic base ("Describes the type-agnostic interface… `TPZStructMatrixT` is the one structural matrices should inherit from", lines 17-25); virtuals `Clone()` + `Create()`; ctors take `TPZCompMesh*` raw or `TPZAutoPointer` (lines 55-68); move ops deleted; destructor non-default due to incomplete-type deletion concerns (`@orlandini` comment, lines 30-46 — candid in-code engineering note).
- `StrMatrix/TPZStructMatrixT.h` — typed layer (`<STATE>/<CSTATE>`).
- `StrMatrix/TPZStrMatParInterface.h` — parallel-interface base (virtual base of `TPZStructMatrix` [repo:25]).
- Parallel assembly schemes [agent]: `pzstrmatrixor.h` (`TPZStructMatrixOR` — "owner rule"? classic thread-per-color?), `pzstrmatrixot.h` (`TPZStructMatrixOT`), `TPZStructMatrixOMPorTBB.h`, `pzstrmatrixflowtbb.h` (TBB flow-graph). Multiple coexisting strategies — inventory + benchmark relevance in Phases 5/6.
- Storage-specific concrete classes [agent]: `pzskylstrmatrix.h` (skyline), `TPZSpStructMatrix.h` / `TPZSSpStructMatrix.h` (sparse nonsym/sym; MKL Pardiso-backed variants; MUMPS variant `TPZSSpStructMatrixMumps` used by divfreebubbles `iter_elast.cpp:297`), `pzfstrmatrix.h` (full), `pzbstrmatrix.h` (band), `TPZFrontStructMatrix` (frontal), `pzbdstrmatrix.h` (block-diagonal, used for preconditioners).
- `StrMatrix/TPZEquationFilter.h` — restrict assembly/solution to an equation subset (used with iterative solvers and by `TPZMatRedSolver`-style reductions [agent]).

## Validation signal
`UnitTest_PZ/TestStruct` + `TestMultithreading` assert parallel == serial matrices and known-matrix assembly [agent] → Phase 7.

## Related
[[assembly]] · [[TPZAnalysis]] · [[matrix-and-solvers]] · [[TPZCompMesh]] · [[static-condensation]]

## OR vs OT — RESOLVED (Phase 5 sweep [agent, structure spot-verified])
- **OR** (`pzstrmatrixor.cpp`): producer/consumer. Workers pull elements one-by-one under a mutex (`NextElement`, :829-844), compute `CalcStiff` into fresh per-iteration `ek/ef`, hand results to a **single consumer thread** that alone writes the global matrix (:714) — no matrix locking needed, but the consumer serializes scatter.
- **OT** (`pzstrmatrixot.cpp`): **graph coloring**. Precomputed `fElSequenceColor`/`fElBlocked`; threads take work via `fCurrentIndex->fetch_add(1)` (:695), keep stack-local ek/ef, wait on a condition variable until their blocking predecessor completed (:821), then scatter **concurrently without atomics** — safe because coloring guarantees no DOF overlap (:848-872).
- Materials are shared mutable across threads in both → [[finding-thread-shared-materials]]. Parallel==serial equality unit-tested.

## Open questions
- How equation filters interact with condensed meshes and connect sequence numbers.
- OMPorTBB/flow-TBB variants: maintained or experimental? (PerfTests stale; no benchmark evidence.)
