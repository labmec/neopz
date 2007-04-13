//$Id: pzinterpolationspace.cpp,v 1.1 2007-04-13 13:14:44 tiago Exp $

#include "pzinterpolationspace.h"
#include "pzmaterialdata.h"
#include "pzelmat.h"
#include "pzquad.h"

TPZInterpolationSpace::TPZInterpolationSpace()
 : TPZCompEl(){}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy) 
 : TPZCompEl(mesh, copy){}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int,int> &gl2lcElMap)
  : TPZCompEl(mesh, copy, gl2lcElMap){}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int &index)
  : TPZCompEl(mesh, copy, index){}

TPZInterpolationSpace::TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index)
 : TPZCompEl(mesh,gel,index){}

TPZInterpolationSpace::~TPZInterpolationSpace(){}

int TPZInterpolationSpace::MaxOrder(){
  const int n = this->NConnects();
  int result = 0;
  int side;
  for(int i = 0; i < n; i++){
    side = this->Connect(i).Order();
    if (side > result) result = side;
  }//i
  return result;
}

void TPZInterpolationSpace::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                                         TPZFMatrix &jacobian, TPZFMatrix &axes,
                                         REAL &detjac, TPZFMatrix &jacinv,
                                         TPZFMatrix &phi, TPZFMatrix &dphix){
  TPZGeoEl * ref = this->Reference();
  if (!ref){
    PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
    return;
  }//if
  TPZFNMatrix<660> dphi;
  int dim = this->Dimension();
  int nshape = NShapeF();

  ref->Jacobian( intpoint, jacobian, axes, detjac , jacinv);
  this->Shape(intpoint,phi,dphi);
  int ieq;
  switch(dim){
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
    PZError << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n";
  } //switch
  ref->X(intpoint, X);
}

REAL TPZInterpolationSpace::InnerRadius(){
  if (!this->Reference()){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Reference() == NULL\n";
    return 0.;
  }
  return this->Reference()->ElementRadius();
}

void TPZInterpolationSpace::InitMaterialData(TPZMaterialData &data){
  this->Material()->FillDataRequirementsInterface(data);
  const int dim = this->Dimension();
  const int nshape = this->NShapeF();
  const int nstate = this->Material()->NStateVariables();
  data.phi.Redim(nshape,1);
  data.dphix.Redim(dim,nshape);
  data.axes.Redim(3,3);
  data.jacobian.Redim(dim,dim);
  data.jacinv.Redim(dim,dim);
  data.x.Resize(3);
  if (data.fNeedsSol){
    data.sol.Resize(nstate);
    data.dsol.Redim(dim,nstate);
  }
}//void

void TPZInterpolationSpace::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
  if (data.fNeedsNeighborSol){
    this->ComputeSolution(qsi, data.normal, data.soll, data.dsoll, data.axesleft, data.solr, data.dsolr, data.axesright);
  }//fNeedsNeighborSol

  if (data.fNeedsSol){
    if (data.phi.Rows()){//if shape functions are available
      this->ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    }
    else{//if shape functions are not available
      this->ComputeSolution(qsi, data.sol, data.dsol, data.axes);
    }
  }//fNeedsSol

  if (data.fNeedsHSize){
    data.HSize = 2.*this->InnerRadius();
  }//fNeedHSize
}//void

void TPZInterpolationSpace::CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) {
  int i;

  TPZAutoPointer<TPZMaterial> material = Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    ek.Reset();
    ef.Reset();
    return;
  }
  int numdof = material->NStateVariables();
  int ncon = NConnects();
  int dim = Dimension();
  int nshape = NShapeF();

  int numeq = nshape*numdof;
  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);
  TPZManVector<REAL> sol(numdof,0.);
  for (i = 0; i < ncon ; i++)	{
    ek.fBlock.Set(i,NConnectShapeF(i)*numdof);
    ef.fBlock.Set(i,NConnectShapeF(i)*numdof);
  }

  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);

  for(i=0; i<ncon; ++i){
    (ef.fConnect)[i] = ConnectIndex(i);
    (ek.fConnect)[i] = ConnectIndex(i);
  }
  //suficiente para ordem 5 do cubo
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

  TPZIntPoints &intrule = GetIntegrationRule();
  if(material->HasForcingFunction()) {
    TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
    intrule.SetOrder(order);
  }

  int intrulepoints = intrule.NPoints();
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){

    intrule.Point(int_ind,intpoint,weight);
    this->ComputeShape(intpoint, x, jacobian, axes, detjac, jacinv, phi, dphix);
    weight *= fabs(detjac);
    if (material->NeedsSolutionToContribute()){
      this->ComputeSolution(intpoint, phi, dphix, axes, sol, dsol);
    }
    material->Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphix,ek.fMat,ef.fMat);

  }//loop over integratin points
}

