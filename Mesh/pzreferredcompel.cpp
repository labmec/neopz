// -*- c++ -*-

// $Id: pzreferredcompel.cpp,v 1.15 2007-04-19 11:41:14 tiago Exp $


#include "pzreferredcompel.h"
#include "pzelctemp.h"
#include "pzintel.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"

#include "pzquad.h"
#include "TPZGeoElement.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h" 
#include "TPZGeoLinear.h"
#include "pzshapelinear.h"
#include "pzgeoquad.h"
#include "pzshapequad.h"
#include "pzgeotriangle.h"
#include "pzshapetriang.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "pzgeoprism.h"
#include "pzshapeprism.h"
#include "pzgeopyramid.h"
#include "pzshapepiram.h"
#include "pzgeotetrahedra.h"
#include "pzshapetetra.h"

#include "pzmaterial.h"
#include "pzelmat.h"
#include "pzgeoel.h"
#include "pzcmesh.h"

#include "pzlog.h"



#include "pztempmat.h"
#include "tpzcompmeshreferred.h"

#ifdef DEBUG
  #define DEBUG2
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
#endif

using namespace std;

template< class TCOMPEL>
TPZCompEl * TPZReferredCompEl<TCOMPEL>::ReferredElement(){
  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (!refmesh) return NULL;
  TPZCompEl * other = refmesh->ReferredEl( this->Index() );
  return other;
}

template< class TCOMPEL>
void TPZReferredCompEl<TCOMPEL>::Print(std::ostream & out){
  out << "\n" << __PRETTY_FUNCTION__ << "\n";
  TCOMPEL::Print(out);

  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (refmesh){
    TPZCompEl * other = refmesh->ReferredEl( this->Index() );
    out << "My ReferredEl = " << other << "\n";
  }
  else out << "My ReferredEl = " << 0 << "\n";
  out << "end of " << __PRETTY_FUNCTION__ << "\n";
}//void

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index):TCOMPEL(mesh, gel,index){

}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::~TPZReferredCompEl(){

}//method

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol,
                                                      TPZFMatrix &dsol, TPZFMatrix &axes){
  TPZCompEl * other = this->ReferredElement();
  if (!other) return;

  TPZManVector<REAL> ThisSol(sol);
  TPZFNMatrix<100> ThisDSol(dsol);

  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol,OtherDSol2;
  TPZFNMatrix<9> otheraxes(3,3,0.);
  other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
  if(sol.NElements()){
    AdjustSolutionDerivatives(OtherDSol,otheraxes,OtherDSol2,axes);
  }
  else if(OtherSol.NElements()){
    OtherDSol2 = OtherDSol;
    axes = otheraxes;
  }
  Append(ThisSol,OtherSol,sol);
  Append(ThisDSol,OtherDSol2,dsol);
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol,
                                                      TPZFMatrix &dsol, const TPZFMatrix &axes){
  TPZCompEl * other = this->ReferredElement();
  if (!other) return;

  TPZManVector<REAL> ThisSol(sol);
  TPZFNMatrix<100> ThisDSol(dsol);

  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol,OtherDSol2;
  TPZFNMatrix<9> otheraxes(3,3,0.);
  other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
  AdjustSolutionDerivatives(OtherDSol,otheraxes,OtherDSol2,axes);
  Append(ThisSol,OtherSol,sol);
  Append(ThisDSol,OtherDSol2,dsol);
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi,
                                                       TPZVec<REAL> &normal,
                                                       TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol, TPZFMatrix &leftaxes,
                                                       TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  TPZCompEl * other = this->ReferredElement();
  if (!other) return;

  TPZManVector<REAL> ThisLeftSol(leftsol), ThisRightSol(rightsol);
  TPZFNMatrix<100> ThisDLeftSol(dleftsol), ThisDRightSol(drightsol);

  TPZManVector<REAL> OtherLeftSol(0), OtherRightSol(0), OtherNormal(0);
  TPZFNMatrix<100> OtherDSol2(0), OtherDLeftSol(0), OtherDLeftSol2(0), OtherDRightSol(0), OtherDRightSol2(0);
  TPZFNMatrix<9> OtherLeftAxes(3,3,0.), OtherRightAxes(3,3,0.);
  other->ComputeSolution(qsi, OtherNormal,
                         OtherLeftSol,  OtherDLeftSol,  OtherLeftAxes,
                         OtherRightSol, OtherDRightSol, OtherRightAxes);

  if (OtherLeftSol.NElements() || OtherRightSol.NElements()){//it means other has solution left/right
    if (normal.NElements() && OtherNormal.NElements()){ //then both element must have same normal
      if ( !AreEqual(normal,OtherNormal) ){
        PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
      }
    }
    if (normal.NElements() == 0){//however, this may be a interpolationspace and other is an interface.
       normal = OtherNormal;//Then OtherNormal is the corret value
    }//if (normal.NElements() == 0)
  }//if other has solution

  if(leftsol.NElements()){
    AdjustSolutionDerivatives(OtherDLeftSol,OtherLeftAxes,OtherDLeftSol2,leftaxes);
  }
  else if(OtherLeftSol.NElements()){
    OtherDLeftSol2 = OtherDLeftSol;
    leftaxes = OtherLeftAxes;
  }

  if(rightsol.NElements()){
    AdjustSolutionDerivatives(OtherDRightSol,OtherRightAxes,OtherDRightSol2,rightaxes);
  }
  else if(OtherRightSol.NElements()){
    OtherDRightSol2 = OtherDRightSol;
    rightaxes = OtherRightAxes;
  }

  Append(ThisLeftSol, OtherLeftSol, leftsol);
  Append(ThisDLeftSol, OtherDLeftSol, dleftsol);
  Append(ThisRightSol, OtherRightSol, rightsol);
  Append(ThisDRightSol, OtherDRightSol, drightsol);
}

template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi,
                                                   TPZFMatrix &phi,
                                                   TPZFMatrix &dphix,
                                                   const TPZFMatrix &axes,
                                                   TPZVec<REAL> &sol,
                                                   TPZFMatrix &dsol){
  TCOMPEL::ComputeSolution(qsi, phi, dphix, axes, sol, dsol);
  this->AppendOtherSolution(qsi, sol, dsol, axes);
}//method


template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi,
                                                   TPZVec<REAL> &normal,
                                                   TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol, TPZFMatrix &leftaxes,
                                                   TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  TCOMPEL::ComputeSolution(qsi, normal, leftsol, dleftsol, leftaxes, rightsol, drightsol, rightaxes);
  this->AppendOtherSolution(qsi, normal, leftsol, dleftsol, leftaxes, rightsol, drightsol, rightaxes);
}

void AdjustSolutionDerivatives(TPZFMatrix &dsolfrom, TPZFMatrix &axesfrom,
                               TPZFMatrix &dsolto, const TPZFMatrix &axesto)
{
  TPZFNMatrix<9> axesinner, axesfromlocal;
  axesfrom.Transpose(&axesfromlocal);
  axesto.ConstMultiply(axesfromlocal,axesinner);
  int nderiv = dsolfrom.Rows();
  int nstate = dsolfrom.Cols();
  dsolto.Resize(nderiv,nstate);
  int id,jd,is;
  for(is=0; is<nstate; is++)
  {
    TPZManVector<REAL> dval(nderiv,0.);
    for(id=0; id<nderiv; id++)
    {
      for(jd=0; jd<nderiv; jd++)
      {
        dval[id] += dsolfrom(jd,is)*axesinner(id,jd);
      }
    }
    for(id=0; id<nderiv; id++)
    {
      dsolto(id,is) = dval[id];
    }
  }
}

void Append(TPZVec<REAL> &u1, TPZVec<REAL> &u2, TPZVec<REAL> &u12)
{
  int nu1 = u1.NElements(),nu2 = u2.NElements();
  u12.Resize(nu1+nu2);
  int i;
  for(i=0; i<nu1; i++) u12[i] = u1[i];
  for(i=0; i<nu2; i++) u12[i+nu1] = u2[i];
}
  
void Append(TPZFMatrix &u1, TPZFMatrix &u2, TPZFMatrix &u12)
{
  int ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
  int ru12 = ru1 < ru2 ? ru2 : ru1;
  int cu12 = cu1+cu2;
  u12.Redim(ru12,cu12);
  int i,j;
  for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
  for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);
}

bool AreEqual(const TPZVec<REAL> &A, const TPZVec<REAL> &B, REAL tol){
  if (A.NElements() != B.NElements()) return false;
  int i;
  const int n = A.NElements();
  for(i = 0; i < n; i++){
    if ( fabs(A[i] - B[i]) > tol ) return false;
  }
  return true;
}

using namespace pzshape;
using namespace pzgeom;

template class TPZReferredCompEl< TPZInterfaceElement >;
template class TPZReferredCompEl< TPZCompElDisc >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >;

TPZCompEl * CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >(mesh,gel,index);
}

TPZCompEl * CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >(mesh,gel,index);
}

TPZCompEl * CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >(mesh,gel,index);
}

TPZCompEl * CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >(mesh,gel,index);
}

TPZCompEl * CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZCompElDisc >(mesh,gel,index);
}
