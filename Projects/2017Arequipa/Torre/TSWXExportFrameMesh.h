#ifndef TSWXExportFrameMeshH
#define TSWXExportFrameMeshH

#include "pzcmesh.h"
#include "TSWXGraphMesh.h"

/** Class to create TSWXGraphMesh from a TPZCompMesh based on TPZEulerBernoulliBeam elements
 * TSWXGraphMesh can then be exported to paraview.
 */
class TSWXExportFrameMesh{

  private:

  TSWXGraphMesh fGraphMesh;

  public:

  void GenerateGraphMesh(TPZCompMesh &cmesh, double time);

  const TSWXGraphMesh & Mesh() const{ return fGraphMesh; }

};


#endif
