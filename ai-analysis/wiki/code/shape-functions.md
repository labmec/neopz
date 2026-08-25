---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - shape-functions
---

# Shape/ — shape-function engine

## Responsibility
Computes hierarchical shape functions and their derivatives on master elements, for all supported topologies and space families (H1, H(div), H(curl), L²). Doxygen module note: implemented as **static classes with no virtual calls, for efficiency** [agent, consistent with headers].

## Key files [repo paths]
- H1 per topology: `Shape/pzshapelinear.h`, `pzshapequad.h`, `pzshapetriang.h`, `pzshapecube.h`, `pzshapetetra.h`, `pzshapeprism.h`, `pzshapepiram.h`, `pzshapepoint.h` (namespace `pzshape`).
- Newer unified drivers: `Shape/TPZShapeH1.h`, `TPZShapeHDiv.h`, `TPZShapeHDivConstant.h`, `TPZShapeHDivKernel2D.h`, `TPZShapeHCurl.h`, `TPZShapeHCurlNoGrads.h`, `TPZShapeDisc` (?), with `Shape/TPZShapeData.h` as the state carrier (orders, connect ids, precomputed side transforms).
- `Shape/pzgenericshape.h` — generic composition machinery.
- `Shape/TPZEnumApproxFamily.h` — family enums used by [[approx-space-creators]].

## Design notes (verified in Phase 4 for H(div))
- Hierarchical (not Lagrangian-nodal) bases: connect-ordered blocks (vertex/edge/face/internal functions), enabling variable p per connect → hp machinery ([[refinement-hanging-nodes]], [[hp-adaptivity]]).
- Side shape functions restricted to sides support conformity checks (`sideshape_continuity` test [agent]).
- **Verified**: H(div) vector shapes = scalar H1 shape × constant master direction (`TPZShapeHDiv.cpp:345-355`), directions built by Topology (`ComputeHDivDirections`, cached in `TPZShapeData.fHDiv.fMasterDirections`); Shape layer outputs master-element values only — the Piola map lives in the element layer ([[piola-transformations]]). Exactly the published construction ([[devloo-group-shape-construction]]).
- **Verified**: `TPZShapeHDivConstant` derives from `TPZShapeHCurlNoGrads` — per facet one RT0 divergence carrier + divergence-free curls/rotated gradients (`TPZShapeHDivConstant.cpp:129-215`); known FAD-branch inconsistency → [[finding-hdivconstant-fad-index]].
- Orientation: `fSideOrient` signs folded into master directions + facet permutation gather (`HDivPermutation`) — see [[TPZCompElHDiv]].

## Related
[[topology-module]] · [[TPZCompElHDiv]] · [[geometric-mappings]] · [[quadrature]] · [[h1-space]] · [[hdiv-space]] · [[hcurl-space]]

## Open questions
- Where derivatives are mapped to physical space (axes/jacobian application) — shape layer or element layer? (`TPZMaterialData.axes` suggests element layer.)
- Legacy per-topology `Chebyshev`-based orthogonal polynomials vs newer `TPZShapeH1` path: which is live for which element class?
