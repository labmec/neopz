// -*- c++ -*-

//$Id: meshes.h,v 1.1 2006-01-10 19:40:31 tiago Exp $

class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;

#include "pzvec.h"

void SetPOrder(int p);

//Cria dominio triangular com refinamento em torno da ponta da fratura. O dominio eh triangular porque admite simetria entre a fratura de injecao e a fratura de producao
TPZCompMesh * CreateMesh_TriangularDomain_ComoPhilippeQuer(int h, int SingH, int p);

//Cria malha igual a CreateMesh_TriangularDomain_ComoPhilippeQuer. A diferenca esta nos valores da condicao mista. Aqui o reboco varia ao longo do comprimento da fratura.
TPZCompMesh * CreateMesh_TriangularDomain_ComoPhilippeQuer_Adimensional(int h, int SingH, int p);

//Igual a CreateMesh_TriangularDomain_ComoPhilippeQuer_Adimensional mas com dominio quadrado, sem simetria entre injecao e producao
TPZCompMesh * CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(int h, int SingH, int p);

//Calcula valor Val1 do reboco (condicao mista)
void KRebocoVal1(TPZVec<REAL> &loc, TPZFMatrix &result);

//Calcula valor Val2 do reboco (condicao mista)
void KRebocoVal2(TPZVec<REAL> &loc, TPZVec<REAL> &result);

//Cria malha simple com solucao analitica para validar TPZCoupledTransportDarcy
TPZCompMesh * CreateSimpleMeshWithExactSolution(int h, int p);
void ExactSol_p(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void ExactSol_u(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv);
void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &force);
void Forcing2(TPZVec<REAL> &x, TPZVec<REAL> &force);
void Dirichlet1(TPZVec<REAL> &pt, TPZVec<REAL> &force);