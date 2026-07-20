---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: major
evidence-commit: 6ffd38b12
tags:
  - neopz
  - material
  - elasticity
---

# At pinned develop, TPZHybridElasticity2D::Contribute omits the body-force RHS

## Repository evidence
`git diff develop HEAD -- Material/Elasticity/TPZHybridElasticity2D.cpp` [repo]: the working tree **adds** (i.e., develop @ 6ffd38b12 **lacks**) the entire forcing-function contribution to `ef` in `Contribute(...)` — evaluation of `fForcingFunction`/`fAnisotropicForcingFunction` at `datavec[0].x` and the `ef(2*in,col) += weight*(floc[0]*phi(in,0) - …fPreStress…)` loops. Commit history: fix landed as `2f6f7982d` "Added the RHS computation" (between develop tip and origin/develop `852a5116c`).

## Reference evidence
Weak form of elasticity requires ∫ f·v on the RHS; omitting it makes any problem with nonzero body force silently wrong (solution of the homogeneous equation with the given BCs). No reference dispute — this is arithmetic completeness, not a modeling choice.

## Assessment
- **Classification: confirmed implementation bug at the pinned commit — already fixed upstream** (origin/develop ≥ `2f6f7982d`). Not a live defect for users tracking origin/develop; a real defect for anyone pinned at/before `6ffd38b12`.
- Blast radius: only `TPZHybridElasticity2D` (hybrid 2D elasticity); problems with zero body force (e.g. iter_elast's homogeneous case if its ForceFunc ≡ 0) produce correct-looking results, which is exactly why it could slip in — no in-tree test exercises this material with nonzero f (→ Phase 7 gap: [[error-estimation-convergence]]).
- Also in the same delta: `TPZMatrix::MultiplyByScalar` made virtual (pzmatrix.h:190→193) and a **self-assignment guard** added to `TPZSYsmpMatrix::CopyFrom` (`from && from != this`, TPZSYSMPMatrix.h:42-45) — the guard implies a real self-copy path existed at the pin (candidate minor finding; verify caller in Phase 5).

## What would resolve/validate
A regression test: hybrid elasticity 2D with manufactured solution having nonzero body force, asserting convergence — none exists at the pin. Recommend in Phase 7.

## Related
[[material-system]] · [[flow-iter-elast]] · [[hybridization]] · [[error-estimation-convergence]]
