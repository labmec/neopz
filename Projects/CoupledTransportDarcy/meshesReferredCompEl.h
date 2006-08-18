// -*- c++ -*-

//$Id: meshesReferredCompEl.h,v 1.2 2006-08-18 13:35:28 tiago Exp $

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
void KRebocoVal1andVal2(TPZVec<REAL> &loc, TPZFMatrix &val1, TPZVec<REAL> &val2, int &type);

TPZCompMesh * TesteConvectivoPuro(int h, int SingH, int p);
TPZCompMesh * TesteConvectivoPuro2(int h, int p);
