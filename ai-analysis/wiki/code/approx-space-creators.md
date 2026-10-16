---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - approximation
  - hybridization
---

# Approximation-space creators

## Responsibility
Two layers that decide *which computational element type* is instantiated on each geometric element, and orchestrate multi-mesh (multiphysics) space construction, hybridization and condensation.

## Layer 1 — element factory (verified [repo])
`Pre/pzcreateapproxspace.h` — `TPZCreateApproximationSpace` (Devloo, since 2009; header @brief "Administer the creation of approximation spaces"): a table of 8 function pointers `TCreateFunction fp[8]` (one per element topology: point…cube) plus flags (`fCreateHybridMesh`, `fCreateLagrangeMultiplier`, `fCreateWithMemory`) and space-family flavors `HDivFamily/H1Family/HCurlFamily` (pzcreateapproxspace.h:27-53). Styles: `EContinuous, EDiscontinuous, EHDiv, EHCurl, EMultiphysics, ESBFem, …` (line 50). Every [[TPZCompMesh]] embeds one (pzcmesh.h:18 include).
**Correction C1:** an explorer report placed this file in `Mesh/`; it is in `Pre/` (see log).
**OQ1:** copy ctor/`operator=` appear not to copy the family flags / style (lines 58-75, tail unread) — verify in Phase 5.

## Layer 2 — problem-level creators (verified [repo])
`Pre/TPZApproxCreator.h` — abstract base `TPZApproxCreator`: holds `HybridizationType {ENone, EStandard, EStandardSquared, ESemi}`, `ProblemType {ENone, EElastic, EDarcy, EStokes}`, material map, default p-order, `fExtraInternalPOrder` (hdiv+/hdiv++), `fShouldCondense`, `fIsRBSpaces` (rigid-body/constant enrichment enabling full internal condensation), and a nested `HybridizationData` struct managing wrap/interface/Lagrange material ids (TPZApproxCreator.h:15-16,38-100).
Concrete: `Pre/TPZHDivApproxCreator.{h,cpp}` (mixed H(div)×L² spaces), `Pre/TPZH1ApproxCreator.{h,cpp}` (**in the 5-file develop delta** — hybrid H1 spaces; cite only after `git show develop:` cross-check), MHM variants `Pre/TPZMHMHDivApproxCreator.h` / `TPZMHMH1ApproxCreator.h` [agent].
Older/parallel machinery: `Pre/TPZHybridizeHDiv.h` (procedural H(div) hybridization used by e.g. divfreebubbles `2frac`), MHM mesh controllers `Pre/TPZMHMeshControl.h`, `TPZMHMixedMeshControl.h`, `TPZMHMixedHybridMeshControl.h` [agent] → [[mhm]].

## Downstream extension
divfreebubbles derives `TPZH1HybridApproxCreator` (app-side) adding `ComputeOrthogonalizingRestraints`, `HybridizeLowOrderFluxes`, `GroupAndCondenseElements` (used in `targets/iter_elast.cpp:218-233` [repo]) → [[divfree-support-lib]], [[flow-iter-elast]].

## Related
[[TPZCompMesh]] · [[hybridization]] · [[static-condensation]] · [[mixed-methods]] · [[TPZCompElHDiv]] · [[material-system]] · [[mhm]]

## Open questions
- Division of labor between `TPZHybridizeHDiv` (older) and `TPZApproxCreator` hybridization (newer): duplication or complementary? → Phase 4/5.
- Meaning and math of `EStandardSquared` ("squared" hybridization — Lagrange multiplier hybridized twice?) and `ESemi` (semi-hybridization) → [[hybridization]] research, Phase 3, vs Devloo-group papers.
