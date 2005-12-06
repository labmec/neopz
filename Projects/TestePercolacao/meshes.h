// -*- c++ -*-

//$Id: meshes.h,v 1.2 2005-12-06 13:37:24 tiago Exp $

class TPZCompMesh;
class TPZFMatrix;

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
