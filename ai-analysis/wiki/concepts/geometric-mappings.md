---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - geometry
---

# Geometric mappings (master → physical)

**Idea.** Each element is the image of a reference (master) element under a map X(ξ); FE integrals are pulled back via the Jacobian. Linear/multilinear maps for straight elements; higher-order or exact maps for curved geometry.

**In NeoPZ.** Per-topology map classes in `Geom/` (`pzgeom::TPZGeoQuad` etc.) plugged into element templates (`TPZGeoElRefLess<TGeo>`); *blend* maps `tpzgeoblend.h` (transfinite blending of curved boundary reps into element interiors); `SpecialMaps/` exact maps (arc, ellipse, sphere, torus, cylinder, NACA airfoil, quadratic elements). README headline: "non-linear geometrical mappings (curved elements with exact representation)" [repo:README.md:18]. `TPZGeoEl::Jacobian/GradX` deliver the metric; `TestGeometry` (`gradx_tests`) and `TestBlend` (semicircle comparisons) validate [agent].

**Traced (Session 2, [[geometry-refinement-maps]]).** `Jacobian` QR-factors GradX into square `jac` + orthonormal `axes` in the *base class* (non-virtual) — the 2D-in-3D mechanism. Blend maps discover curved neighbors at `Initialize` and add Gordon–Hall deviations weighted by per-topology blend factors; BC elements on curved sides automatically become blend elements. **Children inherit exact maps**: `TPZGeoElMapped` stores child corners in the eldest ancestor's parametric space and composes through it ("if the coarse grid map is consistent, then so will all refined meshes"). `TPZChangeEl` retrofits curvature onto imported meshes (`ChangeToCylinder/Arc3D/GeoBlend/QuarterPoint`) — the workflow wann uses on gmsh wells ([[app-wann]]) and ErrorEstimation on NACA profiles ([[app-error-estimation]], custom `TPZBlendNACA` subclass).
**Still open.** Integration-order adequacy for curved maps ([[quadrature]]); curved × vector-space validation ([[piola-transformations]], Phase-7 gap).

**Reference anchors.** Gordon–Hall blending; Devloo-group curved H(div) paper (hdivCurvedJCompAppMath); Ern–Guermond ch. on geometry.

Related: [[TPZGeoMesh]] · [[topology-module]] · [[piola-transformations]] · [[quadrature]]
