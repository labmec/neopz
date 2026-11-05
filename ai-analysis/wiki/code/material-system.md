---
type: code
status: draft
updated: 2026-07-02
confidence: high
evidence-commit: 6ffd38b12
tags:
  - neopz
  - material
  - weak-form
---

# Material system — weak forms & constitutive models

## Responsibility
A "material" is NeoPZ's unit of physics: it evaluates the weak form (`Contribute`: element matrix + rhs at an integration point), boundary conditions (`ContributeBC`), post-processing variables (`Solution`), and optionally exact-solution errors. Materials are keyed by material-id and attached to [[TPZCompMesh]].

## Architecture (verified [repo])
Layered + variadic-mixin design:
- `Material/TPZMaterial.h` — type-agnostic root ("Actual materials should derive from TPZMatBase").
- `Material/TPZMaterialT.h` — `TPZMaterialT<TVar>` type-parametrized layer (STATE vs CSTATE).
- `Material/TPZMatBase.h:21-23` — `template<class TVar, class... Interfaces> class TPZMatBase : public TPZMaterialT<TVar>, public virtual Interfaces...`. One mandatory space interface: `TPZMatSingleSpaceT` (one approximation space) or `TPZMatCombinedSpacesT` (multiphysics). Optional capability mixins: error computation (`TPZMatError*`), load cases, integration-point memory (`TPZMatWithMem`, plasticity/history), eigen problems, interface (DG) contributions.
- `Material/TPZBndCond(Base,T).h` — boundary conditions are themselves materials created via `TPZMatBase::CreateBC` (TPZMatBase.h:59-68) referencing the volumetric material.
- `Material/TPZMaterialData(T).h` — per-integration-point data carrier (shape values, gradients, axes, solution) passed into `Contribute` → [[assembly]].

## Physics families (dirs under `Material/`) [repo dirs; class lists agent-cited]
`Poisson/` (`TPZMatPoisson`), `DarcyFlow/` (`TPZDarcyFlow` primal, `TPZMixedDarcyFlow` H(div)×L², hybrid + fracture variants), `Elasticity/` (`TPZElasticity2D/3D`, `TPZMixedElasticityND`, `TPZHybridElasticity2D/3D` — the 2D hybrid one is in the 5-file develop delta, `TPZHybridMixedElasticityUP`), `Projection/` (L²/H(div)/H(curl) projections), `Electromagnetics/` (waveguides + PML), `Plasticity/` (~75 headers, `BUILD_PLASTICITY_MATERIALS`-gated), `ConsLaw/` (Euler), `BlackOil/`.
Glue materials: `TPZNullMaterial(CS)` (space placeholder), `TPZLagrangeMultiplier(CS)` (interface coupling in [[hybridization]]).

## Breadth map (Session 2, agent-traced, key lines verified [✓])
- **Compiled physics dirs** [✓ `Material/CMakeLists.txt` add_subdirectory list]: Plasticity (`BUILD_PLASTICITY_MATERIALS`-gated, double-only), Elasticity, ConsLaw (Euler), BlackOil, Projection, Poisson, Electromagnetics, DarcyFlow.
- **Galerkin vs mixed is visible in the Contribute signature**: single-space `Contribute(TPZMaterialDataT&, …)` (`TPZDarcyFlow`, NEvalErrors 3) vs combined `Contribute(TPZVec<TPZMaterialDataT>&, …)` (`TPZMixedDarcyFlow`, NEvalErrors 5). **datavec order = mesh-vector order** (resolved; `pzmultiphysicscompel.cpp:836-838` + [[multiphysics-composition]]). H1-hybrid materials sit astride both bases: `TPZHybridDarcyFlow : TPZDarcyFlow + TPZMatCombinedSpacesT + TPZMatErrorCombinedSpaces` [✓ `TPZHybridDarcyFlow.h:25`]; same pattern for `TPZHybridElasticity2D/3D`.
- **Scalar types**: `TVar`-templated with STATE+CSTATE instantiations (Poisson, projections, null/Lagrange materials); STATE-only (Darcy, Elasticity, ConsLaw, BlackOil); **CSTATE-only: all of Electromagnetics** [✓ `TPZWgma.h:20-23`] — waveguide modal analysis via eigen mixins `TPZMatGeneralisedEigenVal`/`TPZMatQuadraticEigenVal`, scattering + PML decorator `TPZMatPML<TMAT>`. The complex path is exercised by materials + `TPZEigenAnalysis`, not by any surveyed downstream app ([[apps-overview]]).
- **Post-processing protocol** (consumed by [[post-processing-vtk]]): `VariableIndex(name)` / `NSolutionVariables(var)` on `TPZMaterial` + typed `Solution(data, var, sol)`; `TPZVTKGenerator::InitFields` resolves names → indices through exactly these seams (`TPZVTKGenerator.cpp:299-321`). Error seam: `TPZMatError<TVar>` holds `fExactSol` + `NEvalErrors/ErrorNames`; pure-virtual `Errors(...)` on the single/combined error interfaces; driven by `TPZInterpolationSpace::EvaluateErrorT` → `TPZAnalysis::PostProcessError`.
- **BC framework**: `TPZBndCond` (type int + material back-pointer) → `TPZBndCondT<TVar>` (`fBCVal1` matrix, `fBCVal2` vector, `fForcingFunctionBC`) → variadic `TPZBndCondBase<TVar, Interfaces::TInterfaceBC...>` stamped out by `TPZMatBase::CreateBC` (C++17 fold `SetMaterialImpl` fan-out). BC-type ints are **per-material conventions** (Darcy: 0=Dirichlet/1=Neumann/2=Robin; Elasticity2D adds 3=directional Dirichlet, 4=stress field).
- **TPZMatWithMem** (integration-point memory): `shared_ptr<TPZAdmChunkVector<TMem>>` + `fUpdateMem`; index flows through `TPZMaterialData::intGlobPtIndex` [✓ `TPZMaterialData.h:141`], assigned by element loops, consumed via `MemItem(i)`. Live users: elastoplasticity (`TPZMatElastoPlastic(2D)`), EM sources (`TPZScatteringSrc`, `TPZPlanarWgScattSrc`).

## Legacy layer — corrected (Session 2)
`Material/needrefactor/` = 19 top-level entries + `REAL/` with 108 files [repo count] — old-style pre-mixin materials (CFD/Euler+k-ε, multiphase/reservoir, visco/poro-elastic, plates/shells, biharmonic, Poisson variants…). **Correction C3: it is *not* compiled into `pz`** — no `add_subdirectory(needrefactor)` exists and `libpz.dylib` contains none of its symbols [✓ verified by nm]. Residual risk is include-path shadowing only (headers duplicate modern class names); out-of-library targets (`SubStruct`, PerfTests, Publications, one unit test) still include its headers.

## Related
[[TPZCompMesh]] · [[assembly]] · [[mixed-methods]] · [[hybridization]] · [[approx-space-creators]] · [[error-estimation-convergence]] · [[multiphysics-composition]] · [[apps-overview]]

## Open questions
- Virtual-inheritance diamond (`public virtual Interfaces...`) cost/complexity — Phase 5.
