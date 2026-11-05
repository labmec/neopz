---
type: concept
status: reviewed
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - hdiv
  - hcurl
---

# Piola transformations (vector-basis mapping) — RESOLVED for H(div)

**Idea (reference).** Mapping vector shape functions master→physical must preserve the conforming trace: H(div) uses the *contravariant* Piola map σ_phys = (1/detJ)·J·σ̂ (preserves normal traces & divergence pairing); H(curl) the *covariant* map ([[boffi-brezzi-fortin-2013]] Ch.2).

## What NeoPZ actually does (verified, Phase 4)
The H(div) pipeline **implements the contravariant Piola transform in a split factorization** [repo, verified first-hand]:
1. **Shape layer emits master-element quantities only**: each vector shape = scalar H1 shape × constant master direction vector (`TPZShapeHDiv.cpp:345-355` — `phi = φ_H1·v̂`, `div = ∇̂φ_H1·v̂`); master directions come from Topology with identity gradx (`TPZShapeHDiv.cpp:83-92`), where the topology routine is itself commented "contravariant piola mapping" (`tpztriangle.cpp:1064-1068`) [agent, spot-verified].
2. **Element layer applies the map pointwise**: `gradx.MultAdd(phiMaster, fDeformedDirections, …, 1./fabs(detjac))` and `divphi *= 1/fabs(detjac)` (`Mesh/pzelchdiv.cpp:1032-1033` [repo, read]) — i.e. σ_phys = (1/|detJ|)·J·σ̂, div_phys = div̂/|detJ|.
3. **Sign convention — NeoPZ variant**: uses **|detJ|**, delegating orientation signs to `fSideOrient` (from `TPZGeoEl::NormalOrientation`, `pzelchdiv.cpp:49-53`) folded into the master directions (`TPZShapeHDiv.cpp:104`); facet-DOF neighbor compatibility via topology permutation gather (`HDivPermutation`, `TPZShapeHDiv.cpp:407-459`) [agent, lines cited].
4. **Curved elements**: gradx/detjac evaluated pointwise (no affine shortcut); optional **FAD branch** (`fNeedsDeformedDirectionsFad`, `pzelchdiv.cpp:979-1031` [repo:1026-1030 read]) seeds ∂/∂x via jacinv and re-applies the same Piola map to get exact physical-space derivatives of the mapped basis. Algebraic div scaling is exact for general smooth maps (Piola identity), so no hidden affine assumption [agent derivation note].

`TPZShapeHDivConstant` (constant-divergence family): per facet one RT0 function carries the (constant) divergence; all other functions are divergence-free curls from `TPZShapeHCurlNoGrads` / rotated H1 gradients (`TPZShapeHDivConstant.cpp:129-215`) — matches the flavor semantics hypothesized in [[hdiv-space]].

## Residual expert-validation items (kept open deliberately)
- |detJ| ⊕ `fSideOrient` composition on *all* refinement/orientation configurations (the `NormalOrientation` father-walk was not exhaustively traced) — derivation or targeted test would close it.
- Intentionality of using algebraic divergence (not the computed-but-unused `divphiFad`) on curved elements — maintainer confirmation.
- H(curl) covariant map: **structurally confirmed in Session 2** — `TPZCompElHCurl::TransformShape` applies `phi = axesᵀ·J⁻ᵀ·phî` ("applies covariant piola transform", `TPZCompElHCurl.cpp:564-578`), curl mapped separately (3D `J·curl̂/detJ`); orientation is implicit via node-id transform ids (no `fSideOrient` analog) — see [[element-families]]. Sign-composition scrutiny at HDiv depth remains open.
- Related finding: [[finding-hdivconstant-fad-index]] (FAD branch facet-count inconsistency).

Related: [[hdiv-space]] · [[TPZCompElHDiv]] · [[shape-functions]] · [[geometric-mappings]] · [[devloo-hdiv-variants-accuracy]]
