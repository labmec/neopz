---
type: app-survey
status: reviewed
updated: 2026-07-06
confidence: medium
evidence-commit: "app SymTensor @ 0a135dc (2025-05-19); embedded neopz develop @ 2937b5a90 (2025-05-26)"
tags:
  - neopz
  - downstream
  - elasticity
  - mixed-methods
  - mhm
---

# App survey: MixedElasticity (~/GitHub/MixedElasticityResearch)

> Downstream-usage evidence (Session 2). Claims [agent], load-bearing ones spot-verified [✓]. Citations refer to the app repo.

**What it is.** Research code for **mixed (stress-based, Hellinger–Reissner) elasticity**: stress tensor as the primary H(div) unknown, L2 displacement, and symmetry enforced either **weakly** (rotation/skew multiplier — PEERS/AFW style) or **strongly** (a Johnson–Mercier symmetric-tensor element built at app level). Three subprojects: `Mixed2D` (~11 targets: square/oscillatory/Girkmann/Yotov/HPC4E benchmarks, 3- and 5-field, MHM variants, H1 reference), `MHM-Elas-3D` (`voronoi_mixed_elas` via library creator), `SymTensor` (newest thread, strongly symmetric tensors).

**Library surface exercised:**
- **Tensor-valued H(div)**: stress rows as H(div) vectors with `NStateVariables = dim` — the same element family carrying matrix-valued fields, something scalar Darcy never shows. Displacement = discontinuous L2 via **`TPZCompElDiscScaled : TPZCompElDisc`** [✓ `Mixed2D/TPZCompelDiscScaled.h:17`] — shape functions scaled by element size for conditioning (a downstream subclass of the *discontinuous* element family).
- **Weak symmetry**: rotation/skew multiplier space (nstate 3 in 2D / 6 in 3D) + rigid-body constant multiplier spaces (distributed force, average displacement) — multiphysics meshes of **3, 5, and 7 fields** (`main.cpp:1249`, `main-five.cpp:676-713`, MHM 7-space drivers). Lagrange-level machinery orders their condensation.
- **Strong symmetry at app level (SymTensor)**: Johnson–Mercier macro-elements built from **custom runtime `TPZRefPattern`s** (`TriangleRef`/`QuadRef` → `CreateJohnsonMercier` [✓ `SymTensor/MeshConditioning.cpp:13,30,50`]) + continuous-but-disconnected elements with hand-stitched center-node continuity + a custom tensor interface `TPZInterfaceSymTensor : TPZLagrangeMultiplierCS<STATE>` — a *new space family assembled downstream from library primitives* (refinement patterns used as a space-construction device, not just adaptivity). Includes rigid-body-mode verification via `SolveEigenProblem` on the assembled matrix (`main-sym.cpp:329-347`).
- **All three space-construction idioms in one repo**: manual `SetAllCreateFunctions*` builders (`TPZMixedElasticityCMeshCreator`), the library `TPZHDivApproxCreator` (`voronoi_mixed_elas.cpp:94-150` [✓ `ProblemType::EElastic` :96; rigid-body spaces off :97]), and `TPZHybridizeHDiv` procedural hybridization + `TPZCompMeshTools::GroupElements/CondenseElements` + manual `TPZElementGroup`/`TPZCondensedCompElT` (`MeshConditioning.cpp:539-578`).
- **MHM at scale**: `TPZMHMeshControl`/`TPZMHMixedMeshControl` with `BuildComputationalMesh(substruct=true)` [✓ `main-seven.cpp:1277`] — the *controller* generation of MHM (which the divfreebubbles slice bypassed in favor of app-side creators) driving `TPZSubCompMesh` substructuring for 2D/3D elasticity, incl. nearly-incompressible cases.
- Materials upstreamed: the app's `TPZMixedElasticityND` (combined-spaces, Voigt handling, axisymmetric option) has a library twin in `Elasticity/TPZMixedElasticityND` — a documented app→lib migration case (cf. the same-name-class risk in `CPP_TECHNICAL_REVIEW.md` §6). Note the divfreebubbles finding [[finding-voronoi-null-ganalytic]] concerns a same-named `voronoi_mixed_elas` driver — the lineage spans repos.
- Geometry: gmsh (incl. curved Girkmann dome via `TPZArc3D`+blend), `TPZGenGrid3D`, honeycomb meshes; struct-matrix breadth incl. `TPZParFrontStructMatrix` (frontal solver actually used downstream).

**What it teaches about the library.** (1) The multiphysics + Lagrange-level design scales to 5–7 coupled fields without library changes. (2) Refinement patterns double as macro-element space constructors. (3) The MHM controller generation is alive downstream even as the creator generation replaces it — both must be treated as supported API. (4) App→lib material migration (`TPZMixedElasticityND`) is the concrete example of how the library grows.

Related: [[mixed-methods]] · [[hybridization]] · [[mhm]] · [[condensation-groups-submeshes]] · [[element-families]] · [[refinement-hanging-nodes]] · [[apps-overview]]
