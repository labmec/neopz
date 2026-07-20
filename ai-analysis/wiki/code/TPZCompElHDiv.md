---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - hdiv
  - hcurl
  - elements
---

# TPZCompElHDiv & vector-space element families

## Responsibility
Computational elements implementing H(div)- and H(curl)-conforming spaces on the mesh: face/edge connects shared between neighbors carry the normal/tangential trace continuity; internal connects carry bubbles.

## Key files [repo paths; internals to be verified in Phases 2/4]
- `Mesh/pzelchdiv.h` — `TPZCompElHDiv<TSHAPE>`: main H(div) element (template over topology).
- `Mesh/pzelchdivbound2.h` — boundary-side H(div) element (`TPZCompElHDivBound2`).
- `Mesh/TPZCompElHDivCollapsed.h` — collapsed-dimension H(div) elements (fracture flow).
- `Mesh/TPZCompElHDivDuplConnects*.h` — duplicated-connect variant (semi-hybridization support).
- `Mesh/TPZCompElKernelHDiv.h`, `TPZCompElKernelHDiv3D.h` — "kernel H(div)": divergence-free subspace elements built from curls/potentials → the divfreebubbles research topic.
- `Mesh/TPZCompElHCurl.h` (+ `TPZCompElHCurlFull` etc.), `Mesh/TPZHCurlEquationFilter.h` — H(curl) family.
- `Mesh/TPZCompElDisc.h` — discontinuous (L²) elements; `Mesh/pzelctemp.h` (`TPZIntelGen`) — generic H1 interpolated element.
- Shape-side counterparts in `Shape/TPZShapeHDiv*.h`, `Shape/TPZShapeHCurl*.h` → [[shape-functions]].

## Space "families" (flavors) [repo]
`Shape/TPZEnumApproxFamily.h` defines `HDivFamily {EHDivStandard, EHDivConstant, EHDivKernel}` (names per usage in divfreebubbles/UnitTests; exact enumerators to re-verify), `H1Family`, `HCurlFamily`. Family selection flows through [[approx-space-creators]] → element constructors → shape classes. `EHDivConstant` = constant-divergence flavor (RT0-like divergence structure with higher-order trace? → to pin down in Phase 3/4 vs Devloo-group papers).

## Continuity & mapping mechanics (verified, Phase 4)
Conformity = shared face connects + a three-part orientation protocol: (1) `fSideOrient[face] = Reference()->NormalOrientation(side)` fixed at construction (`pzelchdiv.cpp:49-53`); (2) signs folded into master directions in the Shape layer (`TPZShapeHDiv.cpp:104`); (3) facet-DOF permutation gather from corner-node ids (`HDivPermutation`, `TPZShapeHDiv.cpp:407-459`) so neighbor facet functions match. Master→physical map = **contravariant Piola with |detJ| convention** applied in `ComputeShape` (`pzelchdiv.cpp:1032-1033` [repo read]); FAD branch for curved-element derivative exactness (`:979-1031`). Full trace: [[piola-transformations]]. Unit tests `drham_check`/`drham_permute_check`/`sideshape_continuity` validate exactly this protocol.

## Related
[[shape-functions]] · [[topology-module]] · [[hdiv-space]] · [[hcurl-space]] · [[de-rham-complex]] · [[piola-transformations]] · [[approx-space-creators]] · [[mixed-methods]]

## Session 2 note
The sibling families are now documented in [[element-families]] (H1/HCurl/discontinuous/interfaces). Contrast worth remembering: HDiv orientation = explicit `fSideOrient` sign array; HCurl = implicit node-id transform ids; H1 = shared connects only. HCurl's covariant transform confirmed (`TransformShape`) — see [[piola-transformations]].

## Open questions
- Relationship between `EHDivKernel` elements and `TPZCompElKernelHDiv` (same thing? flavor vs class split?).
