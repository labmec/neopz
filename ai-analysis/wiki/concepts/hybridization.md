---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - hybridization
---

# Hybridization (standard / "squared" / semi)

**Idea.** Break inter-element continuity of a space and reimpose it weakly via Lagrange multipliers living on the mesh skeleton. Benefits: block-diagonal element problems → [[static-condensation]] to a skeleton system; multipliers often have physical meaning (trace pressure / displacement). Classic theory: Cockburn–Gopalakrishnan–Lazarov unified hybridization; HDG as descendant.

**In NeoPZ (Phase 4, traced & verified).** First-class citizen: `HybridizationType {ENone, EStandard, EStandardSquared, ESemi}` (Pre/TPZApproxCreator.h:15 [repo]); `HybridizationData` manages wrap/interface/Lagrange matids (allocated in strides above max mesh matid, `TPZApproxCreator.cpp:38-58`) and multiplier signs `fMultipliers{left,right,2nd-left,2nd-right}`: H1 path Darcy {1,1,1,−1}, Elastic {−1,−1,−1,1} (`TPZApproxCreator.cpp:780-795` [repo, read]); HDiv path branches per type, e.g. Darcy+ESemi {1,−1,1,−1} (`:798-830` [agent]). Geometry: wrap geoels on every interior facet + interface geoels + Lagrange geoels registered in `fInterfaces` (`AddHybridizationGeoElements`, `TPZApproxCreator.cpp:60-263`); glue = `TPZMultiphysicsInterfaceElement` + `TPZLagrangeMultiplierCS` (`:337-368`).

Verified semantics per type [agent traces, key lines spot-verified]:
- **EStandard (H1)**: broken H1 volume + wrap comp-els *sharing volume connects* + skeleton flux space (HDivStandard on `fLagrangeMatId`, order = default) — one multiplier level; `EAvSol`-level connects explicitly kept out of condensation (`develop:TPZH1ApproxCreator.cpp:758-767`).
- **EStandardSquared (H1, iter_elast)**: literally hybridization², via `AddHybridSquareGeoElements` (`TPZApproxCreator.cpp:578-777`): second interface pair + second Lagrange layer; second-level *primal* multiplier space built inside the L2/H1 atomic mesh (`develop:TPZH1ApproxCreator.cpp:309-332`, Lagrange level `EHybFlux`); condensation groups then absorb the first-level flux+Lagrange DOFs into volume groups (`AssociateElements` numloops=2 + interface-connect propagation, `develop:…:816-881`) — **only the second-level skeleton stays global**. Matches the "double-hybrid" concept of [[avancini-2025-double-hybrid-elasticity]].
- **ESemi (HDiv, dupl_connects)**: requires `EHDivConstant/EHDivOptimized` (`TPZHDivApproxCreator.cpp:80-83`); flux mesh built with `TPZCompElHDivDuplConnects*` (each facet connect split into even=constant-flux + odd=higher-order); `SemiHybridizeDuplConnects` (`:1239-1299`) rebinds **only the even/constant connect, on the sideOrient==−1 side of each interior facet**, to the wrap element — higher-order facet functions remain strongly continuous; multiplier submesh order 0. Matches the semi-hybridization of [[carvalho-2024-semi-hybrid-stokes]] structurally (which trace continuity is weakened differs: here it's the *constant normal flux* that becomes multiplier-mediated — variant, not textbook copy).
- Elasticity-in-HDiv needs a rotation space for weak symmetry (scalar in 2D / 3-vector in 3D, `TPZHDivApproxCreator.cpp:603-606`); condensed HDivConstant elasticity is guarded (`:85-89` — "singular K00" comment) unless ESemi or rigid-body spaces supply the missing constant modes.

**Invariants to check.** Multiplier space order vs trace order (matching for stability); transmission conditions assemble with correct signs (left/right interface materials); condensed skeleton system SPD-ness (iter_elast solves it with LDLt/Schur-CG [repo]).

**Reference evidence (Phase 3).** Frame: [[cockburn-2009-unified-hybridization]] (multiplier = trace unknown; condensed SPD skeleton system; variants are legitimate design choices). `EStandardSquared` ↔ primal *double-hybrid* elasticity with H(div)–L² pair and weak tangential continuity via shear-traction multipliers ([[avancini-2025-double-hybrid-elasticity]], CMAME 2025 — mapping hypothesis-level until Phase 4). `ESemi` ↔ *semi-hybridization*: strong normal continuity kept, tangential/partial coupling weak via traction multiplier, realized with duplicated connects ([[carvalho-2024-semi-hybrid-stokes]], IJNME 2024).

Related: [[static-condensation]] · [[mixed-methods]] · [[approx-space-creators]] · [[flow-iter-elast]] · [[mhm]]
