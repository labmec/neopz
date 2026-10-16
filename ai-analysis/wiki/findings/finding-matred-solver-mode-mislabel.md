---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - divfreebubbles
  - solver
  - benchmarking
---

# App repo: iter_elast benchmarks elasticity with the Darcy-shaped preconditioner (mode mislabel with measurable effect)

## Repository evidence [agent trace, spot-checkable]
- `targets/iter_elast.cpp:287` passes `TPZMatRedSolver<STATE>::EDarcyH1Hybrid` for a 2D elasticity problem; `EElasticityH1Hybrid` exists (`divfree/TPZMatRedSolver.h:15`).
- `fProblemOrigin` is consulted in 4 places (`divfree/TPZMatRedSolver.cpp:102,125,130,140`): the H1-hybrid **sign-flip branch treats both modes identically** (:102-109 — pure label there), but the **preconditioner block size differs**: Darcy `bsize = ord`, elasticity `bsize = 2*(ord+1)-3` (:123-134) feeding the `TPZBlockDiagonal` CG preconditioner (:148-151). `nstate` is set but never used (dead).
- 3D path allows only `EDarcyHDiv`; anything else `DebugStop()`s (:140-146).

## Runtime relevance
[run @ 852a5116c(+)]: the observed mesh-independent 19-iteration CG counts were obtained with the Darcy blocking (bsize=1·ord at p=1). The benchmark's `t2` column *is* the quantity this mislabel perturbs; solution correctness is unaffected (any SPD preconditioner leaves the CG fixed point unchanged).

## Assessment
Classification: **confirmed bug — application repo**, consequential for *measurements* (preconditioner-shape mismatch in a solver benchmark), not for correctness. Severity minor. Caveat: switching to the elastic `bsize=3` requires `bsize | nEqHigh` in `TPZBlockDiagonal::Initialize` — worth checking divisibility before the one-line fix.

## What would resolve
Re-run the sweep with `EElasticityH1Hybrid` and compare iteration counts/`t2`. One-line change at `iter_elast.cpp:287`.

## Related
[[flow-iter-elast]] · [[matrix-and-solvers]] · [[divfree-support-lib]] · [[finding-rusage-memory-units]] (same benchmark's memory column)
