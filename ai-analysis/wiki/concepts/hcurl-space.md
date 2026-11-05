---
type: concept
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - hcurl
---

# H(curl)-conforming spaces

**Idea.** Vector fields with curl in L²; conformity = continuity of *tangential* components across faces/edges. Standard for electromagnetics (edge/Nédélec elements).

**In NeoPZ (Session-2 trace, [[element-families]]).** `TPZCompElHCurl<TSHAPE> : TPZIntelGen` with **no vertex connects** (one connect per non-vertex side); shapes `TPZShapeHCurl` / `TPZShapeHCurlNoGrads`. `HCurlFamily {EHCurlStandard, EHCurlNoGrads}` is a **live runtime switch** (unlike H1Family): `NConnectShapeF`/`InitMaterialData`/`ComputeShape` branch per family; NoGrads filters gradient fields from the standard basis (kernel-oriented; reused by `TPZShapeHDivConstant`, [[hdiv-space]]).
**Covariant Piola confirmed structurally**: `TransformShape` applies `phi = axesᵀ·J⁻ᵀ·phî` (comment "applies covariant piola transform", `TPZCompElHCurl.cpp:564-578`); curl mapped separately (3D `J·curl̂/detJ`, 2D `curl̂/detJ`) — the counterpart of the contravariant trace in [[piola-transformations]].
**Orientation protocol differs from HDiv**: implicit, via corner-node-id transform ids parameterizing `ComputeHCurlDirections` — no explicit sign array (HDiv uses `fSideOrient`). Hanging nodes: dedicated vector-valued `RestrainSideT` (L2 trace projection; DebugStops when small-side order < large-side order).
**Boundaries of support**: pyramid and point elements unavailable (`HCURL_EL_NOT_AVAILABLE` → DebugStop); effective HCurl order can exceed the nominal order on quad-type sides (`MaxOrder` override). Materials: `Material/Electromagnetics/` — all **CSTATE** (complex): `TPZWgma`/`TPZAnisoWgma`/`TPZPeriodicWgma` (waveguide modal analysis via the generalised/quadratic eigen mixins), `TPZScattering(+Src)`, PML decorators `TPZMatPML<TMAT>` ([[material-system]]).

**Usage note (Session 2).** None of the five surveyed downstream apps uses H(curl) ([[apps-overview]] §1) — in-tree users are the electromagnetics materials + `TestHCurl`/`TestDeRham`; downstream H(curl) work lives in older lines (WGMAResearch). Coverage claims about H(curl) should lean on the unit suites, not application evidence.

**Reference anchors.** De Siqueira–Devloo–Gomes (construction); Nédélec families via Boffi–Brezzi–Fortin / Monk.

Related: [[element-families]] · [[hdiv-space]] · [[de-rham-complex]] · [[shape-functions]] · [[topology-module]] · [[piola-transformations]]
