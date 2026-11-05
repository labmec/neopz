---
type: code
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - neopz
  - topology
---

# Topology/ — master-element combinatorics

## Responsibility
Defines the *combinatorial* structure of each reference element ("sides" = vertices, edges, faces, volume), side dimensions, side-to-side closure relations, parametric transforms between sides, and node/face permutation machinery. This layer underlies geometry ([[geometric-mappings]]), shape functions ([[shape-functions]]) and conformity of vector spaces ([[TPZCompElHDiv]]).

## Key files [repo paths]
`Topology/tpzpoint.h`, `tpzline.h`, `tpztriangle.h`, `tpzquadrilateral.h`, `tpztetrahedron.h`, `tpzcube.h`, `tpzprism.h`, `tpzpyramid.h` (namespace `pztopology`), plus `TPZTopologyUtils.h`.

## The "side" abstraction (NeoPZ-specific, load-bearing)
Every topological entity of an element is a numbered *side* (e.g. quadrilateral: 4 vertex sides + 4 edge sides + 1 face side = 9). Sides index: connects in [[TPZCompMesh]], neighbor lists in [[TPZGeoMesh]], side transforms (`TPZTransform`), side integration rules ([[quadrature]]), and restraint construction. Unit test `TestTopology` validates face-orientation data structures, transform projections, and constant div/curl reproduction per topology [agent].

## Notable
- Permutation tables for sides support orientation-independent conformity — validated by `drham_permute_check` style tests [agent].
- The pyramid is present as a topology; its H(div)/shape support has historically been special (mixed families sometimes exclude it) — check which families support pyramids (Phase 4).

## Related
[[shape-functions]] · [[geometric-mappings]] · [[TPZGeoMesh]] · [[quadrature]] · [[hdiv-space]]

## Open questions
- Exact encoding of face orientations (local-to-global side permutation id) and where vector-shape sign flips are applied.
- `TPZTransform` side-to-side composition rules — read `TPZTopologyUtils` in Phase 4.
