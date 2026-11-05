---
type: code
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - elements
  - h1
  - hcurl
  - discontinuous
---

# Element families beyond H(div): H1, H(curl), discontinuous, interfaces

Session-2 deep dive (agent trace, load-bearing lines re-verified [✓]). Companion to [[TPZCompElHDiv]] — fills in the families Session 1 under-covered.

## H1 continuous — `TPZCompElH1<TSHAPE>`
`TPZCompElH1 : TPZIntelGen<TSHAPE> : TPZInterpolatedElement` (`Mesh/TPZCompElH1.h:11-12`, `pzelctemp.h:19-20`). **One connect per topological side including corners** (`fConnectIndexes[TSHAPE::NSides]`, `NConnects()==NSides`); side↔connect maps forward to the topology. `SetSideOrder` (`TPZCompElH1.cpp:164-213`) propagates order changes to block sizes and neighbor integration rules; `EffectiveSideOrder` = max over contained sub-sides (:216-243) — the face≥edge order rule. Shape work delegates to `TPZShapeH1<TSHAPE>` (corner × generating-function blend + internal). Hanging nodes: inherits the generic scalar `RestrainSide` — full support.
**H1Family subtlety** [✓ `pzcreateapproxspace.cpp:111-121`]: `fh1fam` is stored but never branched on at runtime — `EH1WidePrism` is resolved *at creation* by instantiating `TPZCompElH1<TPZShapeWidePrism>` instead of `<TPZShapePrism>`. Only prisms have a variant; for all other topologies the family is a no-op.

## H(curl) — `TPZCompElHCurl<TSHAPE>`
Same `TPZIntelGen` skeleton, different DOF model: **no vertex connects** (`fConnectIndexes[NSides−NCornerNodes]`, `TPZCompElHCurl.h:24-25`); connects self-built in `CreateHCurlConnects` (`.cpp:598-618`). `HCurlFamily` **is a live runtime switch** (`fhcurlfam`): `NConnectShapeF`/`InitMaterialData`/`ComputeShape` branch to `TPZShapeHCurl` vs `TPZShapeHCurlNoGrads` with deliberate no-default switches (`.cpp:215-228,309-318,374-382`).
- **Covariant Piola** applied in `TransformShape` [✓ `.cpp:564-578`, comment "applies covariant piola transform"]: `phi = axesᵀ·J⁻ᵀ·phî`; curl transformed separately (`TransformCurl`, 3D: `J·curl̂/detJ`; 1D/2D: `curl̂/detJ`) — closing the "H(curl) covariant trace not yet done" item in [[piola-transformations]] at the structural level.
- **Orientation: implicit via node-id transform ids** (`TPZShapeHCurl::Initialize` → `GetTransformId` → `ComputeHCurlDirections`), vs HDiv's explicit `fSideOrient` sign array — the two vector families solve the same problem by different protocols.
- Own vector-valued `RestrainSideT` (L2 trace projection, `.cpp:620+`; DebugStops if small-side order < large-side order :671-672) — hanging nodes supported via a dedicated path.
- Dead-ends: **pyramid and point unavailable** [✓ `HCURL_EL_NOT_AVAILABLE` → DebugStop, `pzcreateapproxspace.cpp:726-728`]; map-clone ctor DebugStops ("never tested, better safe than sorry", `TPZCompElHCurl.cpp:59`).

## Discontinuous — `TPZCompElDisc`
Derives **directly from `TPZInterpolationSpace`** (not `TPZInterpolatedElement`): **a single connect for the whole element** [✓ `TPZCompElDisc.cpp:377-383`]. Modal/orthogonal basis about the element center (`TPZShapeDisc`, types {ETensorial, EOrdemTotal, …Full}), normalizing constant `fConstC`, optional evaluation in global X coords, and **appendable external shape functions** (`fExternalShape`) — the enrichment hook. No side connects ⇒ no restraints; continuity is weak (interfaces). Created via `SetAllCreateFunctionsDiscontinuous()` (all topologies → `TPZCompElDisc::CreateDisc`). Downstream subclass: `TPZCompElDiscScaled` (element-size scaling, [[app-mixed-elasticity]]).

## How mixed meshes realize L2 pressure (nuance worth remembering)
[✓ `TPZHDivApproxCreator::CreateL2Space`, `TPZHDivApproxCreator.cpp:465-478`]: for p>0 the "L2" mesh is **broken-H1** — `SetAllCreateFunctionsContinuous()` + `ApproxSpace().CreateDisconnectedElements(true)`, where disconnection is achieved by `ResetReference()` right after each element's creation (`pzcreateapproxspace.cpp:232-234`, flag `fCreateHybridMesh`), so neighbors can't share connects; for p=0 (and always for `EHDivConstant`) it is a genuine order-0 `TPZCompElDisc`. There is **no** `EDisconnected` enum. `TPZL2Projection(/CS/HDiv/HCurl)` are *materials* (projection weak forms), orthogonal to the element choice.

## Interface elements & the DG path
- Single-space: `TPZInterfaceElement : TPZCompEl` (`TPZInterfaceEl.h:29`) stores left/right `TPZCompElSide`s + center normal; requires the material to implement `TPZMatInterfaceSingleSpace` (dynamic_cast, `.cpp:256-257`) and calls `ContributeInterface(data, dataleft, dataright, …)` — classic DG jump/flux terms. Auto-created by `TPZInterpolationSpace::CreateInterfaces` where a neighbor is discontinuous (`pzinterpolationspace.cpp:760-849`).
- Multiphysics: `TPZMultiphysicsInterfaceElement` + `TPZMatInterfaceCombinedSpaces` (`TPZMultiphysicsInterfaceEl.cpp:335-336,424`) — the glue of [[hybridization]] and of downstream Lagrange couplings ([[app-wann]], [[app-iterative-saddle-point]]).
- **The DG recipe is compositional**, not a dedicated creator: discontinuous (or broken-H1) space → `AutoBuild` → `TPZCreateApproximationSpace::CreateInterfaceElements` (`pzcreateapproxspace.cpp:1243-1268`) → interface-capable material. `AutoBuildContDisc` supports mixed continuous+discontinuous partitions. (No `SetAllCreateFunctionsDiscontinuousReferred` exists.)

## Dispatch summary (`Pre/pzcreateapproxspace.{h,cpp}`)
8-slot `std::function` table `fp[8]` per topology; `CreateCompEl` switches on `gel->Type()` (:1059-1091); style tracked in `fStyle {ENone, EContinuous, EDiscontinuous, EHDiv, EHCurl, EMultiphysics, EMultiphysicsSBFem, ESBFem, ECustom}`. Full `SetAllCreateFunctions*` inventory: Continuous(+WithMem), Discontinuous, HDiv(+DuplConnects, +Pressure — the latter `#ifndef STATE_COMPLEX`), HCurl(+WithMem), SBFem(+Multiphysics, LAPACK-gated), MultiphysicElem(+WithMem), custom table (`SetCreateFunctions`). Family flavors (`fh1fam/fhdivfam/fhcurlfam`) are captured into the creation lambdas. An abandoned `SetAllCreateFunctionsHDivFull` block sits commented at :927-969.

Related: [[TPZCompElHDiv]] · [[shape-functions]] · [[approx-space-creators]] · [[TPZConnect]] · [[discontinuous-l2-dg]] · [[h1-space]] · [[hcurl-space]]
