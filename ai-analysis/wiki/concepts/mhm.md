---
type: concept
status: draft
updated: 2026-07-02
confidence: low
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - multiscale
---

# MHM — Multiscale Hybrid-Mixed method

**Idea.** Multiscale method: hybridize on a coarse skeleton; solve local (fine-mesh) problems per coarse cell that upscale fine-scale behavior into the coarse system. Developed in the Brazilian FEM community (Harder–Paredes–Valentin line; Devloo-group implementations/extensions).

**In NeoPZ.** `Pre/TPZMHMeshControl.h`, `TPZMHMixedMeshControl.h`, `TPZMHMixedHybridMeshControl.h` (controller generation), `Pre/TPZMHMHDivApproxCreator.h`/`TPZMHMH1ApproxCreator.h` (creator generation) [agent]; substructuring via `TPZSubCompMesh` ([[static-condensation]]). App-side: `TPZMHMGeoMeshCreator` + `TPZMHMHDivApproxCreator` in [[divfree-support-lib]], exercised by [[flow-mhm-hdivconstant]] (polygonal coarse cells from quadtree import, `EHDivConstant` family, rigid-body spaces, `PutinSubstructures`/`CondenseElements`, `EMHMSparse` Schur solver) [agent/repo pending trace].

**Invariants to check (Phase 4).** Coarse-skeleton flux continuity; local-problem well-posedness (constant/rigid-body handling per subdomain); upscaled system SPD-ness; consistency between the two generations (controllers vs creators — duplication?).

**Reference evidence (Phase 3).** Origin: [[araya-2013-mhm]] (SINUM 2013) — coarse-skeleton multipliers, independent local problems, locally conservative dual, subdomain constant/rigid-body kernels as coarse unknowns (explains `IsRigidBodySpaces()=true` in the slice). Group variant for elasticity on polygonal meshes: [[devloo-mhm-elasticity-polygonal]] (displacement & stress-divergence superconvergence).

Related: [[hybridization]] · [[static-condensation]] · [[hdiv-space]] · [[flow-mhm-hdivconstant]]
