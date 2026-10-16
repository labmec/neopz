---
type: source
status: reviewed
authority: strong
scope: paper
updated: 2026-07-02
confidence: medium
tags:
  - neopz
  - elasticity
  - hybridization
---

# Primal double-hybrid elasticity with H(div)–L² spaces (Avancini et al. 2025) + related analysis

## Bibliographic reference
G. Avancini, N. Shauer, H.L. Oliveira, P.R.B. Devloo, "A primal double-hybrid FEM for 3D compressible and incompressible elasticity using H(div)-L2 spaces", Comput. Methods Appl. Mech. Engrg. (2025), art. 118295 (ADS bibcode 2025CMAME.44618295A; [ResearchGate](https://www.researchgate.net/publication/391164187)).
Related analysis (same community): G. Taraschi, M.R. Correa, "Numerical analysis of a locking-free primal hybrid method for linear elasticity with H(div)-conforming stress recovery", [arXiv:2601.21635](https://arxiv.org/html/2601.21635) (2026) — (u, m, p) primal-hybrid + pressure; inf-sup and locking-free coercivity proofs (fetched abstract + structure).

## Why this source matters
Published counterpart of the elasticity-hybridization line the working tree is actively developing (branch SemiHybridElasticity; delta files `TPZH1ApproxCreator`, `TPZHybridElasticity2D`) and the mandated slice [[flow-iter-elast]].

## Claims extracted
- Displacements approximated in **H(div)** (normal component continuous by construction); **tangential continuity imposed weakly** via a Lagrange multiplier = shear stress on edges/faces; pressure in L², De Rham-compatible with the displacement space → stable for compressible through incompressible regimes [2025 abstract].
- "Double-hybrid": two multiplier levels in the formulation (verify exact second hybridization against the paper body in Phase 4).
- Taraschi–Correa: existence/uniqueness under kernel-coercivity + inf-sup (constants independent of λ → locking-free); element-wise recovery gives H(div)-conforming, locally equilibrated, weakly symmetric stress.

## Applicability to NeoPZ
Interprets `HybridizationType::EStandardSquared` and `HybridizationData::SetProblemHybridH1` semantics ([[approx-space-creators]], [[hybridization]]); explains why iter_elast/semiHybrid_elas hybridize *elasticity* and condense to skeleton systems; grounds Phase 4 invariants (multiplier order pairing, SPD-ness of condensed system, incompressible-limit behavior).

## Limits of applicability
**Mapping is hypothesis-level until Phase 4**: iter_elast drives a `TPZH1Hybrid…` (H1 displacement?) creator while the 2025 paper uses H(div) displacements — the code may implement the paper's method, a primal-H1 double hybrid variant, or both selected by ProblemType/family. Do not attribute the paper's stability claims to the code path before tracing.

## Related wiki pages
[[hybridization]] · [[flow-iter-elast]] · [[approx-space-creators]] · [[mixed-methods]] · [[hdiv-space]]

## Open questions
- Which exact spaces does `EStandardSquared` produce for ProblemType::EElastic at the pin? (Phase 4 trace of `TPZH1ApproxCreator` + `HybridizationData::SetProblemHybridH1`.)
- Does "orthogonalizing restraints" (app-side) correspond to a construction in the 2025 paper or to unpublished work?
