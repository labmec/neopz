---
type: code
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - mesh
  - multiphysics
---

# TPZMultiphysicsCompMesh — combining approximation spaces

Session-2 deep dive (agent trace, load-bearing lines re-verified [✓]). The machinery that lets NeoPZ compose independent "atomic" meshes (flux, pressure, rotation, multipliers, …) into one coupled discretization — used by every mixed/hybrid path and, downstream, for 3–7-field couplings ([[app-mixed-elasticity]]) and same-physics background+enrichment composition ([[app-gfem]]).

## Mesh level (`Mesh/TPZMultiphysicsCompMesh.{h,cpp}`)
- State: `m_active_approx_spaces` (0/1 flags) + `m_mesh_vector`, equal length enforced (`.h:22-25`, `.cpp:83-95`). Inactive spaces contribute data but no equations.
- `BuildMultiphysicsSpace` (`.cpp:75-101`): reset references → multiphysics create-functions → `AutoBuild` skeleton → `AddElements` → `AddConnects` → `LoadSolutionFromMeshes` → `ComputeNodElCon` → `CleanUpUnconnectedNodes`.
- `AddElements` (`.cpp:229-321`): per active space, attach the referred atomic element to each multiphysics element; if no same-level match it **walks the geometric ancestry** (:264-280) to accept a coarser atomic element — the hook MHM-style spaces rely on.
- `AddConnects` (`.cpp:323-421`): concatenates atomic connect vectors with per-space offsets `FirstConnect[i]`, and **re-offsets dependency master indices** [✓ :366-377] so hanging-node restraints survive the merge — an easy-to-break invariant for any new build path.

## Element level (`Mesh/pzmultiphysicscompel.cpp`)
- `CalcStiffT` (:811-907) builds one `TPZMaterialDataT` per space (vector sized to `fElementVec`, :836-838); `InitMaterialDataT` (:667-723) marks `fActiveApproxSpace` and calls the material's `FillDataRequirements(dataVec)`; `ComputeRequiredData` (:498-533) evaluates all spaces at the shared integration point, reusing space 0's Jacobian (:514-517); the combined material's `Contribute(datavec, weight, ek, ef)` closes the loop — answering the old open question in [[material-system]]: **datavec order = mesh-vector order**.
- `InitializeElementMatrix` (:544-607) stacks per-space connect blocks; total nstate = sum over spaces.
- Interfaces: `TPZMultiphysicsElement::CreateInterfaces/CreateInterface/RemoveInterfaces` (`pzmultiphysicselement.h:130-138`) create `TPZMultiphysicsInterfaceElement` glue; `TPZBuildMultiphysicsMesh::AddWrap` builds hybridization wrapper stacks.

## Solution transfer (two symmetric paths)
- Instance: `LoadSolutionFromMeshes` / `LoadSolutionFromMultiPhysics` (`.cpp:436-541`) — block copies via `FirstConnectIndex` offsets; the latter finishes with `LoadSolution` on each atomic mesh and skips `NElConnected()==0` connects.
- Static: `TPZBuildMultiphysicsMesh::TransferFromMeshes / TransferFromMultiPhysics` (`Pre/pzbuildmultiphysicsmesh.cpp:305-403/405+`) — map each multiphysics connect to its `(atomic mesh, connect)` pair; `TransferFromMeshes` **recurses into `TPZSubCompMesh` children** (:394-402), i.e. it understands substructured multiphysics meshes.

## Notes
- The multiphysics mesh is itself subclassable downstream (`TPZMultiPhysicsMeshWindow` in [[app-error-estimation]]).
- Composition is not limited to different physics: same-space composition (background + enrichment H1) works because activity flags and connect offsets are per-mesh, not per-space-type ([[app-gfem]]).

Related: [[TPZCompMesh]] · [[TPZConnect]] · [[mixed-methods]] · [[material-system]] · [[condensation-groups-submeshes]] · [[approx-space-creators]]
