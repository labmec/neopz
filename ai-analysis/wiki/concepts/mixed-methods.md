---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
---

# Mixed methods (saddle-point formulations)

**Idea.** Approximate two fields simultaneously (e.g. Darcy: flux σ ∈ H(div) and pressure p ∈ L²) yielding a saddle-point system; stability requires inf-sup-compatible space pairs; payoff = locally conservative fluxes / direct stress approximation.

**In NeoPZ.** Multiphysics machinery: atomic cmeshes per field combined in `TPZMultiphysicsCompMesh` ([[TPZCompMesh]]); combined-space materials (`TPZMatCombinedSpacesT`): `DarcyFlow/TPZMixedDarcyFlow`, `Elasticity/TPZMixedElasticityND`, `TPZHybridMixedElasticityUP` ([[material-system]]); spaces built by `TPZHDivApproxCreator` with `ProblemType::{EDarcy,EElastic}` ([[approx-space-creators]]); Lagrange-multiplier levels on connects order the condensation. App-side slices: [[flow-dupl-connects]], [[flow-dfreebubbles-1el]], hpc4 (3D mixed elasticity, SPE10-like).

**Resolved (Sessions 1–2).** Mixed elasticity symmetry = **weak symmetry with a rotation/skew multiplier space** (scalar 2D / 3-vector 3D), confirmed both in the creator (`TPZHDivApproxCreator.cpp:603-606`, Session 1) and at scale downstream: [[app-mixed-elasticity]] runs 3-, 5- and 7-field Hellinger–Reissner couplings (stress rows as H(div) vectors with `NStateVariables=dim`, plus rigid-body multiplier spaces), and also builds *strongly*-symmetric Johnson–Mercier tensors at app level. The p>0 L² pressure realization is broken-H1, p=0 is `TPZCompElDisc` ([[discontinuous-l2-dg]]). Combined-space `Contribute` receives datavec in mesh-vector order ([[multiphysics-composition]]).
**Still open.** Flux×pressure order pairing per `HDivFamily` (RT-like vs BDM-like classification); local conservation not directly unit-tested.

**Reference anchors.** Boffi–Brezzi–Fortin (canonical); Devloo et al. mixed-elasticity papers (multiphysics + weak symmetry); Arnold's stress-element literature as contrast.

Related: [[hdiv-space]] · [[hybridization]] · [[static-condensation]] · [[de-rham-complex]] · [[error-estimation-convergence]]
