// -*- c++ -*-

//$Id: meshes.h,v 1.1 2005-11-28 13:49:47 tiago Exp $

class TPZCompMesh;

/** 4 quadrilateral
 */
TPZCompMesh * CreateMesh(int h, int p);
TPZCompMesh * CreateMeshSingularity(int h, int SingH, int p);
TPZCompMesh * CreateMeshSingularityKeepingAspectRatio(int h, int SingH, int p);
TPZCompMesh * CreateMeshSingularityKeepingAspectRatio_TriangularDomain(int h, int SingH, int p);
TPZCompMesh * CreateMeshSingularityKeepingAspectRatio_TriangularDomain2(int h, int SingH, int p);
void SetPOrder(int p);
TPZCompMesh * CreateMesh_TriangularDomain_ComoPhilippeQuer(int h, int SingH, int p);