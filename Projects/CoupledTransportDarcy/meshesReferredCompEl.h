// -*- c++ -*-

//$Id: meshesReferredCompEl.h,v 1.1 2006-06-29 12:28:13 tiago Exp $

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;

#include "pzvec.h"

//Cria malha simples com solucao analitica polinomial igual a CreateSimpleMeshWithExactSolution em meshes.h.
void CreateSimpleMeshesWithExactSolutionToReferredCompEl(TPZVec< TPZCompMesh * > & CompMeshes, int h, int p); 


//Mesh from NeoPZ/Project/TestePercolacao which creates a Darcy flow based on parameters obtained from a Propag simulation
void CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(TPZVec< TPZCompMesh * > & CompMeshes, int h, int SingH, int p);

void KRebocoVal1(TPZVec<REAL> &loc, TPZFMatrix &result);

void KRebocoVal2(TPZVec<REAL> &loc, TPZVec<REAL> &result);

TPZCompMesh * TesteConvectivoPuro(int h, int SingH, int p);
TPZCompMesh * TesteConvectivoPuro2(int h, int p);
