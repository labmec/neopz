---
type: source
status: reviewed
authority: strong
scope: paper
updated: 2026-07-02
confidence: medium
tags:
  - neopz
  - hybridization
  - stokes
  - hdiv
---

# Semi-hybrid-mixed method for Stokes–Brinkman–Darcy with H(div) velocities (Carvalho, Devloo et al. 2024)

## Bibliographic reference
P.G.S. Carvalho, P.R.B. Devloo, et al., "A semi-hybrid-mixed method for Stokes–Brinkman–Darcy flows with H(div)-velocity fields", Int. J. Numer. Methods Eng. (2024). DOI [10.1002/nme.7363](https://onlinelibrary.wiley.com/doi/10.1002/nme.7363). Related talk: "A Semi-Hybrid approximation of the Stokes equations using H(div) spaces" (Devloo, Oden Institute seminar).

## Why this source matters
Defines **semi-hybridization** — the published meaning behind `HybridizationType::ESemi` and the duplicated-connects machinery (`TPZCompElHDivDuplConnects*`, app + lib) exercised by [[flow-dupl-connects]] and the current SemiHybrid* research line.

## Claims extracted
- Velocity in H(div): normal-component continuity kept **strong** ("taken for granted"); **tangential continuity imposed weakly** by a Lagrange multiplier playing the role of tangential traction; pressure discontinuous, divergence-compatible with velocity [abstract].
- "Semi": only part of the interface continuity is moved to multipliers (vs full hybridization which breaks all continuity) → smaller multiplier space, keeps local conservation.
- Stokes–Darcy coupling is natural in this setting; in the Darcy region the weak tangential condition can be dropped.

## Applicability to NeoPZ
Semantics for `ESemi` in [[approx-space-creators]]/[[hybridization]]; rationale for duplicated connects (the duplicated face functions become the weakly-constrained part); context for `TPZMatRedSolver` reductions (skeleton multiplier block = Schur target) in [[flow-dupl-connects]] and [[flow-iter-elast]].

## Limits of applicability
Paper treats Stokes/Brinkman; dupl_connects applies ESemi to *Darcy* where the paper says tangential weak continuity is unnecessary — the code's ESemi-for-Darcy may hybridize the *normal-trace/flux* structure differently. Trace before judging (Phase 4).

## Related wiki pages
[[hybridization]] · [[hdiv-space]] · [[flow-dupl-connects]] · [[divfree-support-lib]]

## Open questions
- Exactly which connects are duplicated for ESemi in `TPZHDivApproxCreator` at the pin, and which material glues them (`TPZLagrangeMultiplier`?).
