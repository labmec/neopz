// -*- c++ -*-

//$Id: meshes.h,v 1.2 2006-03-04 15:36:23 tiago Exp $

#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include <set>

/* CONFIGURACOES DA MALHA PARA CONTINUO/DESCONTINUO */
//Apenas um elemento continuo
void OneContinuous(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int continuousindex = 0);

//Apenas um elemento descontinuo
void OneDiscontinuous(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int discontinuousindex = 0);

//Nao funciona!!!
void DiscontinuousOnBoundaryLayer(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h);

//Auxiliar para DiscontinuousOnBoundaryLayer: procura subelementos a um lado do gel.
void GetSubElements( TPZGeoEl * gel, int side, TPZStack<TPZGeoElSide> &stack );


/* MALHAS */
//4 quadrados perto da camada limite
TPZCompMesh * DiscontinuousOnBoundaryLayer(int h);

//4 quadrados iguais, podendo ser divididos. Podem ser continuos ou descontinuos
TPZCompMesh *CreateMesh(int h);

//4 quadrados iguais, podendo ser divididos. Podem ser continuos e/ou descontinuos
TPZCompMesh *CreateMeshContDisc(int h);

//8 triangulos: Acho que nao estao funcional!
TPZCompMesh *CreateMesh2();

/* DADOS */
//Forcing function e condicoes de contorno Dirichlet
void ForcingFunction(TPZVec<REAL> &pto, TPZVec<REAL> &force);
void Dirichlet_X_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u);
void Dirichlet_Y_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u);

//Solucao exata do problema
void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv);


/** Malha com refinamento na camada limite 14fev2006 **/
TPZCompMesh * RefinedOnBoundLayer(int h, int ref_uniforme );
/** Malha com refinamento na camada limite 14fev2006 **/
void OneContinuousRefinementOnBoundLayer(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int ref_uniforme, int continuousindex, int BLMaterialId);
