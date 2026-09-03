---
type: concept
status: reviewed
updated: 2026-07-06
confidence: high
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - discontinuous
  - dg
---

# Discontinuous (L²) spaces & interface/DG machinery

**Idea.** Piecewise-polynomial spaces with no inter-element continuity: the natural home of mixed-method pressures, multipliers/constants, and discontinuous-Galerkin primal fields. Coupling, where needed, is imposed weakly — via saddle-point structure or interface (jump/penalty) terms.

**In NeoPZ — two distinct realizations** (Session 2, [[element-families]]):
1. **`TPZCompElDisc`** — a true discontinuous element: single connect for the whole element, modal/orthogonal basis about the element center, optional external (enrichment) shape functions, no restraint machinery. Used for order-0 pressure/constant spaces (always for `EHDivConstant`), rotation multipliers, distributed-flux/average spaces — the rigid-body enrichment meshes of [[static-condensation]] are built from it.
2. **Broken-H1** — the standard trick for p>0 L² fields: continuous factory + `CreateDisconnectedElements(true)`, disconnection achieved by resetting geometric references during build so connects are never shared. Local basis = H1 hierarchical basis; global space = L². This is what `TPZHDivApproxCreator::CreateL2Space` produces for p>0 pairs.

**Interface/DG layer.** `TPZInterfaceElement` (single-space) and `TPZMultiphysicsInterfaceElement` (combined-spaces) assemble jump terms through the `TPZMatInterfaceSingleSpace`/`TPZMatInterfaceCombinedSpaces` material interfaces. DG is compositional: discontinuous space + `CreateInterfaceElements` + interface-capable material — there is no monolithic "DG creator". The same interface machinery is what [[hybridization]] uses for Lagrange-multiplier transmission (`TPZLagrangeMultiplier(CS)`), so DG and hybrid methods share one code path.

**Downstream evidence.** Every surveyed mixed app builds discontinuous pressure/multiplier meshes ([[apps-overview]] §1); `TPZCompElDisc` is subclassed downstream for conditioning (`TPZCompElDiscScaled`, [[app-mixed-elasticity]]); interface elements couple 3D/2D/1D physics in [[app-wann]] and tangential tractions in [[app-iterative-saddle-point]].

**What is *not* in-tree.** No upwinding/numerical-flux library for hyperbolic DG in the modern layer (the old `ConsLaw`/`needrefactor` CFD materials carry their own); `TPZAgglomerateElement` (agglomerated DG coarsening) exists but was not on any analyzed path.

Related: [[element-families]] · [[mixed-methods]] · [[hybridization]] · [[hdiv-space]] · [[static-condensation]]
