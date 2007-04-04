// -*- c++ -*-

// $Id: pzreferredcompel.cpp,v 1.11 2007-04-04 19:37:17 tiago Exp $


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

template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix, TPZFMatrix &axes, TPZVec<REAL> &u, TPZFMatrix &du){

  TCOMPEL::ComputeSolution(qsi, phi, dphix, axes, u, du);

  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (!refmesh) return;
  TPZCompEl * other = refmesh->ReferredEl( this->Index() );
  if (!other) return;
  
  TPZManVector<REAL> ThisSol(u);

  TPZFNMatrix<100> ThisDSol(du);
  
  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol,OtherDSol2;
  TPZFNMatrix<9> otheraxes(3,3,0.);
  other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
  AdjustSolutionDerivatives(OtherDSol,otheraxes,OtherDSol2,axes);
  Append(ThisSol,OtherSol,u);
  Append(ThisDSol,OtherDSol2,du);

}//method

 /**
 * Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
  */
template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi, 
                               TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes)
{
  TCOMPEL::ComputeSolution(qsi, sol, dsol, axes);

  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (!refmesh) return;
  TPZCompEl * other = refmesh->ReferredEl( this->Index() );
  if (!other) return;
  
  TPZManVector<REAL> ThisSol(sol);
  TPZFNMatrix<100> ThisDSol(dsol);
  
  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol,OtherDSol2;
  TPZFNMatrix<9> otheraxes(3,3,0.);
  other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
  if(sol.NElements()) 
  {
    AdjustSolutionDerivatives(OtherDSol,otheraxes,OtherDSol2,axes);
  }
  else if(OtherSol.NElements())
  {
    OtherDSol2 = OtherDSol;
    axes = otheraxes;
    
  }
  Append(ThisSol,OtherSol,sol);
  Append(ThisDSol,OtherDSol2,dsol);
}

template < >
void TPZReferredCompEl< TPZInterfaceElement >::ComputeSolution(TPZVec<REAL> &qsi,
                                                               TPZVec<REAL> &sol,
                                                               TPZFMatrix &dsol,
                                                               TPZFMatrix &axes)
{
  TPZManVector<REAL> LeftSol, RightSol;
  TPZFNMatrix<100> LeftDSol, RightDSol;
  this->NeighbourSolution(this->LeftElementSide(),  qsi, LeftSol,  LeftDSol,  axes);
  this->NeighbourSolution(this->RightElementSide(), qsi, RightSol, RightDSol, axes);
  Append(LeftSol, RightSol, sol);
  Append(LeftDSol,RightDSol,dsol);

  this->AppendOtherSolution(qsi, sol, dsol, axes);
}

template < class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::AppendOtherSolution(TPZVec<REAL> &qsi, TPZVec<REAL> &sol,
                                                       TPZFMatrix &dsol, const TPZFMatrix &axes){
  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (!refmesh) return;
  TPZCompEl * other = refmesh->ReferredEl( this->Index() );
  if (!other) return;

  TPZManVector<REAL> ThisSol(sol);
  TPZFNMatrix<100> ThisDSol(dsol);

  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol,OtherDSol2;
  TPZFNMatrix<9> otheraxes(3,3,0.), myAxes(axes);
  other->ComputeSolution(qsi, OtherSol, OtherDSol, otheraxes);
  AdjustSolutionDerivatives(OtherDSol,otheraxes,OtherDSol2,myAxes);
  Append(ThisSol,OtherSol,sol);
  Append(ThisDSol,OtherDSol2,dsol);
}

 /**
   * Computes solution and its derivatives in the local coordinate qsi.
   * @param qsi master element coordinate of the interface element
   * @param leftsol finite element solution
   * @param dleftsol solution derivatives
   * @param leftaxes axes associated with the left solution
   * @param rightsol finite element solution
   * @param drightsol solution derivatives
   * @param rightaxes axes associated with the right solution
  */
template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi, 
                            TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes,
                            TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                            TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes)
{
  ComputeSolution(qsi,sol,dsol,axes);
  ComputeSolution(qsi,leftsol,dleftsol,leftaxes,rightsol,drightsol,rightaxes);
}

template< >
void TPZReferredCompEl< TPZInterfaceElement >::ComputeSolution(TPZVec<REAL> &qsi, 
                                                               TPZVec<REAL> &sol, TPZFMatrix &dsol,TPZFMatrix &axes,
                                                               TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                                                               TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes){
  //ComputeSolution(qsi,sol,dsol,axes);Interface has no solution associated to
  this->ComputeSolution(qsi,leftsol,dleftsol,leftaxes,rightsol,drightsol,rightaxes);
  const int dim = this->Dimension();
  TPZFNMatrix<100> jac(dim,dim), jacinv(dim,dim), ThisAxes(3,3,0.);
  REAL detjac;
  this->Reference()->Jacobian(qsi, jac, ThisAxes, detjac, jacinv);
  axes = ThisAxes;
  this->AppendOtherSolution(qsi, sol, dsol, axes);
}

 /**
   * Computes solution and its derivatives in the local coordinate qsi.
   * This method will function for both volumetric and interface elements
   * @param qsi master element coordinate of the interface element
   * @param sol finite element solution
   * @param dsol solution derivatives
   * @param axes axes associated with the derivative of the solution
   * @param leftsol finite element solution
   * @param dleftsol solution derivatives
   * @param leftaxes axes associated with the left solution
   * @param rightsol finite element solution
   * @param drightsol solution derivatives
   * @param rightaxes axes associated with the right solution
  */
template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi, 
                            TPZVec<REAL> &leftsol, TPZFMatrix &dleftsol,TPZFMatrix &leftaxes,
                            TPZVec<REAL> &rightsol, TPZFMatrix &drightsol,TPZFMatrix &rightaxes)
{
  TCOMPEL::ComputeSolution(qsi, leftsol, dleftsol, leftaxes,rightsol, drightsol, rightaxes);

/*  TPZCompMesh * cmesh = this->Mesh();
  TPZCompMeshReferred * refmesh = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  if (!refmesh) return;
  TPZCompEl * other = refmesh->ReferredEl( this->Index() );
  if (!other) return;
  
  TPZManVector<REAL> LSol(leftsol),RSol(rightsol);
  TPZFNMatrix<100> LDSol(dleftsol),RDSol(drightsol);
  
  TPZManVector<REAL> OLSol,ORSol;
  TPZFNMatrix<100> OLDSol,OLDSol2,ORDSol,ORDSol2;
  TPZFNMatrix<9> OLaxes(3,3,0.),ORaxes(3,3,0.);
  other->ComputeSolution(qsi, OLSol,OLDSol,OLaxes,ORSol,ORDSol,ORaxes);
  if(leftsol.NElements()) 
  {
    AdjustSolutionDerivatives(OLDSol,OLaxes,OLDSol2,leftaxes);
  }
  else if(OLSol.NElements())
  {
    OLDSol2 = OLDSol;
    leftaxes = OLaxes;
    
  }
  Append(LSol,OLSol,leftsol);
  Append(LDSol,OLDSol2,dleftsol);
  if(rightsol.NElements()) 
  {
    AdjustSolutionDerivatives(ORDSol,ORaxes,ORDSol2,rightaxes);
  }
  else if(ORSol.NElements())
  {
    ORDSol2 = ORDSol;
    rightaxes = ORaxes;
    
  }
  Append(RSol,ORSol,rightsol);
  Append(RDSol,ORDSol2,drightsol);*/
  
}


/**
 * Adjust the derivatives from one system of axes to the other
 */
void AdjustSolutionDerivatives(TPZFMatrix &dsolfrom, TPZFMatrix &axesfrom,
                               TPZFMatrix &dsolto, TPZFMatrix &axesto)
{
  TPZFNMatrix<9> axesinner, axesfromlocal;
  axesfrom.Transpose(&axesfromlocal);
  axesto.Multiply(axesfromlocal,axesinner);
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
