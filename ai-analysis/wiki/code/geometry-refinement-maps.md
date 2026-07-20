---
type: code
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - geometry
  - refinement
  - curved-maps
---

# Geometry layer: element hierarchy, refinement patterns, nonlinear maps

Session-2 deep dive (agent trace, load-bearing lines re-verified [✓]). Code-level companion to [[TPZGeoMesh]], [[geometric-mappings]], [[refinement-hanging-nodes]].

## 1. Element hierarchy (policy-based template stack)
`TPZGeoEl` (`Mesh/pzgeoel.h:41`, abstract: X/GradX in REAL and Fad flavors :554-563, Divide/genealogy seams, CreateBCGeoEl) → `TPZGeoElRefLess<TGeo>` (`pzgeoelrefless.h:31`: owns `TGeo fGeo` + neighbor ring; forwards X/GradX to the Geom policy; `Divide` DebugStops) → either `TPZGeoElRefPattern<TGeo>` (`tpzgeoelrefpattern.h:34`: `fSubEl` + `TPZAutoPointer<TPZRefPattern>`; general `Divide` :356-509) or `TPZGeoElement<TGeo,TRef>` (uniform; forwards to `TRef::` static tables). `TPZGeoMesh::CreateGeoElement(..., reftype)` picks uniform vs pattern (`pzgmesh.cpp:1333+`). The `.h.h` files are out-of-line template bodies. `Jacobian` is non-virtual in the base: QR-factors 3×dim `GradX` into square `jac` + orthonormal `axes` (`pzgeoel.h:542-551`) — the mechanism that lets 2D elements live in 3D.

## 2. Refinement: uniform tables vs runtime patterns
- Uniform `TPZRef*` (e.g. `Refine/pzreftriangle.cpp`): pure static data (son corners, mid-side coords, son→father transforms `buildt`, `fatherside`). `NewMidSideNode` reuses a neighbor's midnode if one exists (:153-179). **Pyramid caveat** [✓ `pzrefpyram.h:27`, `.cpp:336-347`]: `NSubEl=10` = 6 pyramids + 4 tets — uniform refinement of pyramids introduces tetrahedra.
- `TPZRefPattern` (`Refine/TPZRefPattern.h:77`) *is* a small `TPZGeoMesh` (element 0 = father, rest = partition) + precomputed side transforms/permutations. `.rpt` format documented at `.h:37-72` (nodes block, elements block, father first). `TPZGeoElRefPattern::Divide` checks neighbor side-pattern compatibility before lazily adopting the uniform pattern (:369-408), delegates node creation to the pattern, and errors loudly on incompatibility (`CreateMidSideNodes` DebugStops if an existing neighbor midnode is >1e-2 off, `TPZRefPattern.cpp:648-674`).
- Matching tools (`Refine/TPZRefPatternTools.cpp`): `GetCompatibleRefPatterns` (:28), `PerfectMatchRefPattern` (:193,437) driven by `SidesToRefine` (:951-999), and `RefineDirectional` (:1001) — the driver wann and GFEM use for well-heel/crack-tip grading ([[app-wann]], [[app-gfem]]); MixedElasticity uses hand-built patterns as macro-element space constructors ([[app-mixed-elasticity]]).
- **Global state**: `gRefDBase` (`TPZRefPatternDataBase.cpp:31`) maps type→patterns and id→pattern; `.rpt` library loads from `PZ_REFPATTERN_DIR` (configure-baked path — see [[finding-build-config-gaps]] relocatability note). Deserializing a refined mesh **re-resolves pattern ids against the live DB** [✓ `tpzgeoelrefpattern.h.h:20-35`] — a saved mesh is unreadable without the same DB populated (persistence coupling, [[persistence]]).

## 3. Genealogy → hanging nodes (the geometry/computation bridge)
`TPZInterpolatedElement::Divide` divides geometry first, then per new element `CreateMidSideConnect` (`pzintel.cpp:652`) asks `EqualLevelElementList` (share connect) or `LowerLevelElementList(1)` (:701,741) — the latter delegating to `TPZGeoElSide::LowerLevelCompElementList2` which **walks `Father2()/StrictFather()` ancestry** (`pzgeoelside.cpp:931`). When a coarser neighbor exists, `RestrainSide` builds the L2-projection dependency ([[TPZConnect]]); the only geometric input is `SideTransform3` (`pzgeoelside.cpp:682+`), which accumulates transforms up the refinement tree via `BuildTransform2`.

## 4. Nonlinear & special maps (two mechanisms + inheritance under refinement)
- **Analytic maps** (`SpecialMaps/`): `TPZArc3D` (circle fit through 3 points, closed-form X/GradX), `TPZCylinderMap<TGeo>` (cylindrical corner coords + rotation), `TPZEllipse3D`, `TPZWavyLine`, tori/spheres. `IsLinearMapping()==false` routes construction to mapped/blend paths.
- **Isoparametric quadratics** (`TPZQuadraticTrig/Quad/…`): quadratic Lagrange shapes over stored midside nodes.
- **`TPZGeoBlend<TGeo>`** (`Geom/tpzgeoblend.{h,cpp}`): Gordon–Hall transfinite blending — linear map + Σ blendFactor·(curved-side map − chord) (`tpzgeoblend.cpp:518+`); discovers curved neighbors in `Initialize` via `SetNeighbourInfo` (`.cpp:71`). Copy ctor DebugStops (`.h:65-72`) — beware mesh-clone paths.
- **Children inherit exact maps** via `TPZGeoElMapped<TBase>` (`Mesh/tpzgeoelmapped.h:29`): stores child corners **in the eldest ancestor's parametric space** [✓ intent comment `.h:24-27`] and composes `X = Xfather(KsiBar(ksi))` — "if the coarse grid map is consistent, then so will all refined meshes". Routing: `CreateGeoElement` dispatches nonlinear elements to `CreateGeoElementMapped` (`pzgeoelrefless.h.h:427-428`); `CreateBCGeoEl` falls back to blend BC elements on curved sides (:326-338).
- **`TPZChangeEl`** (`SpecialMaps/tpzchangeel.h`): in-place surgery — `ChangeToQuadratic/GeoBlend/Arc3D/Cylinder/QuarterPoint` (quarter-point = fracture singularity resolution; used downstream in [[app-error-estimation]], cylinder+blend in [[app-wann]]).

## 5. Sides & neighbors
`TPZGeoElSide` (`pzgeoelside.h:86`) = (element, side); neighbors form **singly-linked circular lists** spliced by `SetConnectivity` (`pzgeoelside.cpp:441-488`); `BuildConnectivities` does bulk discovery. Side-to-side transforms: element-local `SideToSideTransform`, cross-element `NeighbourSideTransform`, tree-walking `SideTransform3` — the conformity substrate for restraints and interfaces. Caveat: `TPZGeoElSideIndex::operator bool()` returns true for an *invalid* side (`pzgeoelside.h:55-58`) — inverted-looking semantics.

## 6. Mesh building (`Pre/`)
`TPZGenGrid2D/3D` (structured, `MMeshType` element choice), `TPZExtendGridDimension` (extrusion), `TPZGmshReader`: physical-name→matid maps per dimension (`TPZGmshReader.h:111-123`), optional tag remap, `InsertElement` constructs uniform or refpattern elements with the physical id as matid (`.cpp:580-649`); undefined-tag elements skipped unless opted in.

Related: [[TPZGeoMesh]] · [[geometric-mappings]] · [[refinement-hanging-nodes]] · [[TPZConnect]] · [[topology-module]] · [[mesh-io-generators]]
