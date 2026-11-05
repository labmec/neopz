---
type: finding
status: reviewed
updated: 2026-07-02
confidence: high
classification: confirmed-bug
severity: minor
evidence-commit: 6ffd38b12
tags:
  - neopz
  - hdiv
  - shape-functions
---

# TPZShapeHDivConstant: FAD branch indexes facet kernel-function counts inconsistently with the REAL branch

## Repository evidence (verified first-hand)
`Shape/TPZShapeHDivConstant.cpp` 3D `Shape(...)`: the REAL branch iterates per facet `i` with `data.fHCurl.fNumConnectShape[nedges + i]` kernel functions (:191); the FAD overload uses `data.fHCurl.fNumConnectShape[nedges]` for **every** facet (:304) — facet 0's count applied to all facets.

## Failure scenario
Trigger requires all of: 3D element, `HDivFamily::EHDivConstant`, the FAD path (`fNeedsDeformedDirectionsFad`, i.e. typically curved/nonlinear geometry), **and per-facet kernel counts that differ** (variable face orders — hp settings). Then the FAD shape tensor mis-partitions functions across facets (wrong count per facet ⇒ shifted `countKernel` association); totals may still hit `count == nshape` only if the sum coincidentally matches, otherwise `DebugStop()` (:211-212 analogue in FAD tail). With uniform face orders the two expressions coincide → benign in all uniform-p runs, which is why tests pass.

## Classification & severity
**Confirmed code inconsistency** (two branches computing the same partition differently — one must be wrong); practical severity *minor* because the triggering combination (curved × EHDivConstant × non-uniform face order × FAD) appears unexercised in-tree (no unit test combines them — Phase 7 gap). Mathematically it is a latent wrong-derivatives bug, not a stability subtlety.

## What would resolve
Unit test: 3D HDivConstant element with mixed face orders + `fNeedsDeformedDirectionsFad=true`, compare FAD shape values against REAL branch / finite differences. Fix is one token (`[nedges]` → `[nedges + i]`) pending maintainer confirmation of intent.

## Related
[[piola-transformations]] · [[hdiv-space]] · [[shape-functions]] · [[TPZCompElHDiv]]
