---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - mesh
  - geometry
---

# TPZGeoMesh — geometric mesh

## Responsibility
Container for the *geometry and topology* of the discretization: geometric elements, nodes, and elementwise/nodal BC identifiers — no degrees of freedom (those live in [[TPZCompMesh]]). Holds the refinement genealogy (father/son element trees) that underpins [[refinement-hanging-nodes]].

## Key files
- `Mesh/pzgmesh.h` / `.cpp` — the mesh: `TPZAdmChunkVector<TPZGeoEl*> fElementVec`, `TPZAdmChunkVector<TPZGeoNode> fNodeVec`, mesh dimension `fDim`, name, and a back-pointer `TPZCompMesh* fReference` (pzgmesh.h:48-70 [repo]).
- `Mesh/pzgeoel.h` — abstract `TPZGeoEl` (element = topology + material id + node indices + neighbor connectivity via sides).
- `Mesh/pzgeoelside.h` — `TPZGeoElSide`: the (element, side) pair used for all neighborhood traversal; neighbors form circular linked lists along shared sides. *(mechanism [agent]; verify in Phase 2 traces)*
- `Mesh/pzgeoelrefless.h(.h)` — `TPZGeoElRefLess<TGeo>`: concrete element without refinement capability; template body in the unusual `.h.h` companion file.
- `Mesh/tpzgeoelrefpattern.h(.h)` — `TPZGeoElRefPattern<TGeo>`: refinable element driven by [[refinement-hanging-nodes|refinement patterns]].
- `Geom/*` + `SpecialMaps/*` — the per-topology map classes `TGeo` plugged into the element templates → [[geometric-mappings]].

## Notable design facts [repo]
- `GMESHNOMATERIAL -9999` sentinel for "no material" (pzgmesh.h:24).
- Interface-material map keyed by material-id pairs (pzgmesh.h:72-75) — supports DG/interface constructions.
- Mutual raw-pointer reference with the computational mesh (`fReference` both ways; pzcmesh.h:49) — the "reference" mechanism used to associate geo↔comp elements. Lifetime/ownership implications → Phase 5.
- Storage is chunked (`TPZAdmChunkVector`) so element/node pointers survive vector growth; free slots are recycled.

## Related
[[TPZCompMesh]] · [[geometric-mappings]] · [[refinement-hanging-nodes]] · [[topology-module]] · [[mesh-io-generators]]

## Open questions
- Exact semantics of `ResetReference/LoadReferences` in multi-cmesh (multiphysics) workflows — needed for the Phase 2 flow traces.
- Ownership: who deletes `TPZGeoEl*` — mesh destructor? (Phase 5 lifetime review.)
