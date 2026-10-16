---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - hdiv
---

# H(div)-conforming spaces

**Idea.** Vector fields whose divergence is in L²; conformity = continuity of the *normal component* across faces. Used for fluxes in mixed formulations (Darcy, mixed elasticity) with local conservation properties. Classic families: Raviart–Thomas, BDM; NeoPZ builds its own hierarchical family.

**In NeoPZ.** Elements [[TPZCompElHDiv]] (+Bound/Collapsed/DuplConnects/Kernel variants); shapes `Shape/TPZShapeHDiv*`; flavors `HDivFamily {EHDivStandard, EHDivConstant, EHDivKernel, EHDivOptimized}` (Shape/TPZEnumApproxFamily.h:5 [repo], default EHDivStandard) selected via [[approx-space-creators]]; `fExtraInternalPOrder` gives hdiv+/hdiv++ enriched internal order (TPZApproxCreator.h:58-59 [repo]).
- *EHDivStandard*: full hierarchical H(div) family (Devloo-group construction).
- *EHDivConstant*: flavor with constant divergence per element (supports rigid-body-mode condensation; related to recent Devloo et al. papers) — hypothesis, verify Phase 3/4.
- *EHDivKernel*: divergence-free subspace (curl of potentials) — the divfreebubbles topic.
- *EHDivOptimized*: unknown semantics — appears in `TestHDivApproxSpaceCreator` GENERATE grid [repo:155]; research Phase 3/4.

**Invariants to check (Phase 4).** Normal-trace continuity incl. orientation sign consistency under face permutations (`drham_permute_check` tests exist [agent]); div maps onto the pressure space exactly ([[de-rham-complex]]); mapping to physical elements (contravariant [[piola-transformations]] or NeoPZ variant); inf-sup stability of chosen flux×pressure pairs ([[mixed-methods]]).

**Reference evidence (Phase 3).** Construction: [[devloo-group-shape-construction]] (JCAM 2013 — geometry-based vectors × hierarchical H1 scalars ⇒ conforming traces by construction). Conformity/inf-sup/Piola expectations: [[boffi-brezzi-fortin-2013]]. Flavors & divergence-order variants are published research ([[devloo-hdiv-variants-accuracy]]); semi-hybrid use: [[carvalho-2024-semi-hybrid-stokes]].

Related: [[TPZCompElHDiv]] · [[mixed-methods]] · [[de-rham-complex]] · [[piola-transformations]] · [[hybridization]] · [[flow-dupl-connects]]
