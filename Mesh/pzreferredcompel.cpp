// -*- c++ -*-

// $Id: pzreferredcompel.cpp,v 1.1 2006-04-18 20:39:54 tiago Exp $

#include "pzreferredcompel.h"
#include "pzelctemp.h"
#include "pzintel.h"
#include "TPZCompElDisc.h"

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

#ifdef LOG4CXX
static LoggerPtr loggerIntel(Logger::getLogger("pz.mesh.tpzinterpolatedelement"));
static LoggerPtr loggerDisc(Logger::getLogger("pz.mesh.tpzcompeldisc"));
#endif

using namespace std;

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index):TCOMPEL(mesh, gel,index){

}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::~TPZReferredCompEl(){

}//method

template<class TCOMPEL>
void TPZReferredCompEl< TCOMPEL >::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef){
  int i;

  if(this->fMaterial == NULL){
    std::cout << "Exiting CalcStiff: no material for this element\n";
    LOGPZ_ERROR(loggerIntel,"Exiting CalcStiff: no material for this element");
    return;
  }
  int numdof = this->fMaterial->NStateVariables();
  int ncon = this->NConnects();
  int dim = this->Dimension();
  int nshape = this->NShapeF();

  int numeq = nshape*numdof;
  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  TPZManVector<REAL> sol(numdof,0.);
  for (i = 0; i < ncon ; i++) {
    ek.fBlock.Set(i,this->NConnectShapeF(i)*numdof);
    ef.fBlock.Set(i,this->NConnectShapeF(i)*numdof);
  }

  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = this->ConnectIndex(i);
    (ek.fConnect)[i] = this->ConnectIndex(i);
  }
  //suficiente para ordem 5 do cubo
  //  REAL phistore[220],dphistore[660],dphixstore[660];
  TPZFNMatrix<220> phi(nshape,1);
  TPZFNMatrix<660> dphi(dim,nshape),dphix(dim,nshape);
  TPZFNMatrix<9> axes(3,3,0.);
  TPZFNMatrix<9> jacobian(dim,dim);
  TPZFNMatrix<9> jacinv(dim,dim);
  REAL detjac;
  TPZManVector<REAL,3> x(3,0.);
  TPZManVector<REAL,3> intpoint(dim,0.);
  REAL weight = 0.;

  //  REAL dsolstore[90];
  TPZFNMatrix<90> dsol(dim,numdof);

  TPZIntPoints &intrule = this->GetIntegrationRule();
  if(this->fMaterial->HasForcingFunction()) {
    TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
    intrule.SetOrder(order);
  }

  int intrulepoints = intrule.NPoints();
  TPZGeoEl *ref = this->Reference();
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){

    intrule.Point(int_ind,intpoint,weight);

    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);

    ref->X(intpoint, x);

    weight *= fabs(detjac);

    this->Shape(intpoint,phi,dphi);

    int ieq;
    switch(dim) {
    case 0:
      break;
    case 1:
      dphix = dphi;
      dphix *= (1./detjac);
      break;
    case 2:
      for(ieq = 0; ieq < nshape; ieq++) {
        dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
        dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
      }
      break;
    case 3:
      for(ieq = 0; ieq < nshape; ieq++) {
        dphix(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
        dphix(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
        dphix(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
      }
      break;
    default:
      stringstream sout;
      sout << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
      LOGPZ_ERROR(loggerIntel,sout.str());
      //PZError.flush();
    }

    this->ComputeSolutionInOtherMesh(phi, dphix, sol, dsol);
    
    this->fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphix,ek.fMat,ef.fMat);

  }
}//method

template<>
void TPZReferredCompEl< TPZCompElDisc >::CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef){

  if(fMaterial == NULL){
    cout << "TPZCompElDisc::CalcStiff : no material for this element\n";
    LOGPZ_ERROR(loggerDisc,"Exiting CalcStiff: no material for this element");
    return;
  }
  TPZGeoEl *ref = Reference();
  int ncon = NConnects();
  int dim = Dimension();
  int nstate = fMaterial->NStateVariables();
  int nshape = NShapeF();
  int numeq = nshape * nstate;

  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  if(ncon){//pode serr no m�imo ncon = 1
    ek.fBlock.SetNBlocks(ncon);
    ef.fBlock.SetNBlocks(ncon); 
    ek.fBlock.Set(0,NShapeF()*nstate);
    ef.fBlock.Set(0,NShapeF()*nstate);
  }
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(int i=0;i<ncon;i++){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  if(ncon==0) return;//elemento CC no passa
  TPZFMatrix phix(nshape,1),dphix(dim,nshape);
  TPZFMatrix axes(3,3,0.);
  TPZFMatrix jacobian(dim,dim);
  TPZFMatrix jacinv(dim,dim);
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL detjac,weight;
  int integ = max( 2 * Degree(), 0);
  TPZIntPoints *intrule = 0;
  intrule = Reference()->CreateSideIntegrationRule(Reference()->NSides()-1,integ);
  if(fMaterial->HasForcingFunction())
  {
     int maxint = intrule->GetMaxOrder();
     TPZManVector<int> order(Reference()->Dimension());
     intrule->GetOrder(order);
     order.Fill(maxint);
     intrule->SetOrder(order);
  }
  int npoints = intrule->NPoints(),ip;                                              //integra fi*fj
  TPZVec<REAL> sol(nstate,0.);
  TPZFMatrix dsol(dim,nstate,0.);

  for(ip=0;ip<npoints;ip++){
    intrule->Point(ip,intpoint,weight);
    ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
    ref->X(intpoint, x);
    weight *= fabs(detjac);
    Shape(x,phix,dphix);
    
    this->ComputeSolutionInOtherMesh(phix, dphix, sol, dsol);
        
    axes.Identity();    
    fMaterial->Contribute(x,jacinv,sol,dsol,weight,axes,phix,dphix,ek.fMat,ef.fMat);
  }
}//method

template<class TCOMPEL>
void TPZReferredCompEl< TCOMPEL >::ComputeSolutionInOtherMesh(TPZFMatrix & phix, TPZFMatrix & dphix, TPZVec<REAL> &sol, TPZFMatrix &dsol){

    //solucao da iteracao anterior
    // In TPZReferredCompEl is possible to connect two different simulations. It is done by using two TPZCompMesh
    // associated to the same TPZGeoMesh. 
    // ref points to the geometric element
    TPZGeoEl * ref = this->Reference();
    if (!ref){
      std::cout << "Error!: " << __PRETTY_FUNCTION__ << std::endl;
      LOGPZ_ERROR(loggerDisc, __PRETTY_FUNCTION__);
    }//if
    // othercompel may be different of this. In that case the we will be able to get the solution of another simullation.
    TPZCompEl * othercompel = ref->Reference();
    if (!othercompel){
      std::cout << "Error!: " << __PRETTY_FUNCTION__ << std::endl;
      LOGPZ_ERROR(loggerDisc, __PRETTY_FUNCTION__);
    }//if    
    
    if (!othercompel->Material()){
      std::cout << "Error!: " << __PRETTY_FUNCTION__ << std::endl;
      LOGPZ_ERROR(loggerDisc, __PRETTY_FUNCTION__);
    }//if    

    const int othernstate = othercompel->Material()->NStateVariables();
    
    int othernshape = -1;
    {
      TPZCompElDisc * otherdisc = dynamic_cast<TPZCompElDisc*>(othercompel);
      if (otherdisc) othernshape = otherdisc->NShapeF();
      else{
        TPZInterpolatedElement * otherintel = dynamic_cast<TPZInterpolatedElement*>(othercompel);
        if (otherintel) othernshape = otherintel->NShapeF();
        else{
          std::cout << "Error!: " << __PRETTY_FUNCTION__ << std::endl;
          LOGPZ_ERROR(loggerDisc, __PRETTY_FUNCTION__);        
        }//else
      }//else
    }//escopo
     
    TPZBlock &otherblock = othercompel->Mesh()->Block();
    TPZFMatrix &otherMeshSol = othercompel->Mesh()->Solution();

    int otherncon = othercompel->NConnects();
    
    sol.Fill(0.);
    dsol.Redim(dphix.Rows(),othernstate);
    dsol.Zero();
    for(int in=0; in<otherncon; in++) {
      TPZConnect *df = &(othercompel->Connect(in));
      int dfseq = df->SequenceNumber();
      
      
      int dfvar = otherblock.Size(dfseq);
      int pos = otherblock.Position(dfseq);
      int iv = 0,d;
      for(int jn=0; jn<dfvar; jn++) {
        sol[iv%othernstate] += phix(iv/othernstate,0)*otherMeshSol(pos+jn,0);
        for(d=0; d<dphix.Rows(); d++){
          dsol(d,iv%othernstate) += dphix(d,iv/othernstate)*otherMeshSol(pos+jn,0);
        }//for d
        iv++;
      }//for jn
    }//for in
}//method


using namespace pzshape;
using namespace pzgeom;

template class TPZReferredCompEl< TPZCompElDisc >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >;

//Create functions
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >(mesh,gel,index);
// }
// 
// template<class TCOMPEL>
// TPZCompEl *TPZReferredCompEl<TCOMPEL>::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >(mesh,gel,index);
// }
// 
// template<>
// TPZCompEl *TPZReferredCompEl<TPZCompElDisc>::CreateDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl< TPZCompElDisc >(mesh,gel,index);
// }

///////////////////////////////////////////****************************/////////////////////////////////////

// template < class T >
// TPZCompEl *TPZReferredCompEl<T>::CreateFunction(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
//   return new TPZReferredCompEl<T>(mesh,gel,index);
// }

///////////////////////////////////////////****************************/////////////////////////////////////

TPZCompEl * CreateSpecialPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >(mesh,gel,index);
}

TPZCompEl * CreateSpecialDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZCompElDisc >(mesh,gel,index);
}
