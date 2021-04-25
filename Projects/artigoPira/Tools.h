//
//  Tools.hpp
//  pz
//
//  Created by Denise De Siqueira on 04/05/2018.
//

#ifndef Tools_hpp
#define Tools_hpp

#include <stdio.h>

#endif /* Tools_hpp */


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "tpzchangeel.h"
#include "tpzarc3d.h"
#include "tpzquadraticquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZStructMatrix.h"
#include "pzstepsolver.h"

#include "TPZGenGrid2D.h"

#include "pzlog.h"

#include "pzl2projection.h"
#include "pzmultiphysicselement.h"
#include "pzintel.h"
#include "TPZVTKGeoMesh.h"

#include "TPZMatDualHybridPoisson.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "tpzgeoelmapped.h"

#include "TPZParFrontStructMatrix.h"
#include "TPZIntQuadQuarterPoint.h"
#include "tpzintpoints.h"
#include "pzquad.h"
#include "tpzhierarquicalgrid.h"
#include "TPZVecL2.h"

#include <iostream>
#include <math.h>


TPZGeoMesh *GMesh();
TPZCompMesh *CMeshFlux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CMeshPressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZMixedPoisson* &mymaterial, bool QuarterPointRule);
TPZCompMesh *CMeshH1(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *MalhaCompH1QP(TPZGeoMesh * gmesh,int ordem);
TPZCompMesh *CMeshFluxL2(TPZGeoMesh *gmesh, int pOrder, int nodeAtOriginId);
void ErroL2NoElemento(TPZCompMesh *hdivmesh, std::ostream &out,  int nodeAtOriginId);

void TransferMatrixFromMeshes(TPZCompMesh *cmesh, TPZCompMesh *MFMesh, TPZAutoPointer< TPZMatrix<STATE> > matF,TPZAutoPointer< TPZMatrix<STATE> > matMP, int nodeAtOriginId);

//with hybrid method
TPZGeoMesh * MalhaGeo(const int h);
TPZCompMesh *CreateHybridCompMesh(TPZGeoMesh &gmesh,int porder,bool ismultiplierH1);
void GroupElements(TPZCompMesh *cmesh);
void Dirichlet2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);


void DirectionalRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void Prefinamento(TPZCompMesh * cmesh, int ndiv,int porder);

void QuarterPointRef(TPZGeoMesh *gmesh, int nodeAtOriginId);

void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void Dirichlet(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

void SolExata3D(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void f_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &f);
void Dirichlet_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannEsquerda_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannDireita_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannAcima_3D(const TPZVec<REAL> &loc, TPZVec<STATE> &result);



void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, int numThreads=0);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out);
void ErroHDivNoElemento(TPZCompMesh *hdivmesh, std::ostream &out,  int nodeAtOriginId);

void ChangeSideConnectOrderConnects(TPZCompMesh *mesh, int reduction_value);
void ChangeInternalConnectOrder(TPZCompMesh *mesh);

void PrefinamentoRibsHybridMesh(TPZCompMesh *cmesh);

void IntegrationRuleConvergence(bool intQuarterPoint);
void NormMax(TPZFMatrix<STATE> A, REAL &val, int &indexi, int &indexj);
void ChangeIntegrationRule(TPZCompMesh *cmesh, int porder,bool IsQPrule);
void GlobalSubMatrix(TPZCompMesh *cmesh, TPZAutoPointer< TPZMatrix<STATE> > mat, int nodeAtOriginId, bool matInicial, std::ofstream &subMat);

void DirectUniRef(TPZGeoMesh *gmesh, int nodeAtOriginId, int divide);

void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, REAL & rho_min);
TPZGeoMesh *GMeshTeste();

void ParametricCircle(REAL &radius,REAL &theta, TPZManVector<REAL,3> & xcoor);
TPZGeoMesh *MeshCircle( int ndiv);
TPZGeoMesh *CurvedMesh();

TPZGeoMesh *CurvedMesh2( int ndiv);
TPZGeoMesh *VerticalExtrusion(TPZGeoMesh *malha2D );
void ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X);

long ComputeAverageBandWidth(TPZCompMesh *cmesh);
