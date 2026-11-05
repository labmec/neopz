---
type: concept
status: draft
updated: 2026-07-02
confidence: medium
evidence-commit: 6ffd38b12
tags:
  - fem
  - neopz
  - visualization
---

# VTK output model

**Idea.** Visualization files for ParaView: legacy `.vtk` (ASCII, simple) vs XML `.vtu` (binary/appended, richer). High-order FE fields must either be subdivided into linear cells or use VTK high-order (Lagrange) cells.

**In NeoPZ.** Modern path: `TPZVTKGenerator` writes legacy-format `.vtk`; each computational element subdivided per `vtkRes` into linear cells (`TPZVTK::CellType` map: point/line/tri/quad/tet/pyr/prism/hex [repo:Post/TPZVTKGenerator.h:30-56]); fields = material `Solution()` variables by name. Legacy path: graph meshes (`pzvtkmesh` and other formats). Geometric-mesh dumps: `TPZVTKGeoMesh` (used for debugging, e.g. the 71MB `GeoMeshHybrid.vtk` in divfreebubbles build dir [repo ls]). Output *meaning*: pointwise evaluations of the FE solution at subdivision nodes — no projection/smoothing (verify), duplicated points across elements allow discontinuous fields (verify representation).

**Invariants to check (Phase 4).** Node ordering per VTK cell type; sub-element node placement for vtkRes>0; scalar/vector/tensor component conventions (ParaView expects 3-vectors); whether discontinuities across element boundaries are representable (point duplication) — matters for L² pressure fields in mixed methods.

**Reference anchors.** VTK legacy file-format spec (Kitware); NGSolve vtkoutput (declared origin of the adaptation [repo header]).

Related: [[post-processing-vtk]] · [[material-system]] · [[flow-dfreebubbles-1el]] · [[flow-mhm-hdivconstant]]
