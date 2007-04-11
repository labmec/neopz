//$Id: pzmaterialdata.cpp,v 1.1 2007-04-11 14:26:23 tiago Exp $

#include "pzmaterialdata.h"
#include "pzmaterial.h"
#include "pzcompel.h"
#include "pzelmat.h"
#include <sstream>
#include "pzerror.h"
#include "TPZInterfaceEl.h"

TPZMaterialData::TPZMaterialData(TPZCompEl &cel, TPZMaterial &material){
  this->fEl = &cel;
  this->fMat = &material;
  this->SetAllRequirements(true);
}

void TPZMaterialData::InitializeAttributes(TPZElementMatrix & ek, TPZElementMatrix & ef){
  const int numdof = this->fMat->NStateVariables();
  const int ncon = this->fEl->NConnects();
  const int dim = this->fEl->Dimension();
  const int nshape = this->fEl->NShapeF();
  const int numeq = nshape*numdof;
  ek.fMat.Redim(numeq,numeq);
  ef.fMat.Redim(numeq,1);
  ek.fBlock.SetNBlocks(ncon);
  ef.fBlock.SetNBlocks(ncon);

  for(int i = 0; i < ncon ; i++){
    ek.fBlock.Set(i,this->fEl->NConnectShapeF(i)*numdof);
    ef.fBlock.Set(i,this->fEl->NConnectShapeF(i)*numdof);
  }//i
  ek.fConnect.Resize(ncon);
  ef.fConnect.Resize(ncon);
  for(int i = 0; i < ncon; i++){
    (ef.fConnect)[i] = this->fEl->ConnectIndex(i);
    (ek.fConnect)[i] = this->fEl->ConnectIndex(i);
  }
  this->phi.Redim(nshape,1);
  this->dphi.Redim(dim,nshape);
  this->dphix.Redim(dim,nshape);
  this->axes.Resize(3,3);
  this->axes.Zero();
  this->jacobian.Redim(dim,dim);
  this->jacinv.Redim(dim,dim);
  this->X.Resize(3,0.);

  this->sol.Resize(numdof);
  this->sol.Fill(0.);
  this->dsol.Redim(dim,numdof);
}//void

TPZMaterialData::~TPZMaterialData(){
//NOTHING TO BE DONE!
}

int TPZMaterialData::AveragePOrder(TPZCompEl * cel){
  const int n = cel->NConnects();
  int result = 0;
  for(int i = 0; i < n; i++){
    result += cel->Connect(i).Order();
  }
  return static_cast<int>(result/n);
}//int

void TPZMaterialData::FillData(TPZVec<REAL> &qsi, REAL &weight){
    TPZGeoEl * ref = this->fEl->Reference();
    if (!ref){
      PZError << "\nError at " << __PRETTY_FUNCTION__ << " - ref == NULL\n";
    }
    REAL detjac;
    ref->Jacobian( qsi, jacobian, axes, detjac , jacinv);
    weight *= fabs(detjac);
    this->fEl->Shape(qsi,phi,dphi);

    const int nshape = this->fEl->NShapeF();
    int ieq, dim;
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
      PZError << "\nError at " << __PRETTY_FUNCTION__ << " - please implement the " << dim << "d Jacobian and inverse\n";
    } //switch

    if ( this->fNeedsSol ){
      if ( this->fNeedsNeighborSol ){
        this->fEl->ComputeSolution(qsi, this->sol, this->dsol, this->axes,
                                   this->normal,
                                   this->leftsol, this->leftdsol, this->leftaxes,
                                   this->rightsol, this->rightdsol, this->rightaxes);
      }
      else{
        this->fEl->ComputeSolution(qsi, this->phi, this->dphix, this->axes, this->sol, this->dsol);
      }
    }//fNeedsSol

    if ( this->fNeedsX ){
      ref->X(qsi, this->X);
    }//fNeedsX

    if (this->fNeedsPOrder){
      this->POrder = this->AveragePOrder(this->fEl);
      TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(this->fEl);
      if (face){
        this->leftPOrder = this->AveragePOrder(face->LeftElement());
        this->rightPOrder = this->AveragePOrder(face->RightElement());
      }//if
      else{
        this->leftPOrder = 0;
        this->rightPOrder = 0;
      }//else
    }//fNeedsPOrder

    if (this->fNeedsHSize){
      this->Hsize = 2. * this->fEl->Reference()->ElementRadius();
      TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(this->fEl);
      if (face){
        this->leftHsize = 2. * face->LeftElement()->Reference()->ElementRadius();
        this->rightHsize = 2. * face->RightElement()->Reference()->ElementRadius();
      }//if
      else{
        this->leftHsize = 0.;
        this->rightHsize = 0.;
      }//else
    }//fNeedsHSize

}//void


void TPZMaterialData::Contribute(TPZVec<REAL> &qsi, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  this->fMat->FillMaterialDataRequirements(*this);
  this->FillData(qsi, weight);
  if ( this->fNeedsNeighborSol ){
    this->fMat->Contribute(this->X, this->jacinv, this->sol, this->dsol, this->axes,
                           this->normal,
                           this->leftsol, this->leftdsol, this->leftaxes,
                           this->rightsol, this->rightdsol, this->rightaxes,
                           weight,
                           this->phi, this->dphix,
                           ek, ef);
  }
  else{
    this->fMat->Contribute(this->X, this->jacinv, this->sol, this->dsol, weight, this->axes, this->phi, this->dphix, ek, ef);
  }
}//void

void TPZMaterialData::ContributeInterface(TPZVec<REAL> &qsi, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  this->fMat->FillMaterialDataRequirements(*this);
  this->FillData(qsi, weight);
  if ( this->fNeedsNeighborSol ){
    this->fMat->Contribute(this->X, this->jacinv, this->sol, this->dsol, this->axes,
                           this->normal,
                           this->leftsol, this->leftdsol, this->leftaxes,
                           this->rightsol, this->rightdsol, this->rightaxes,
                           weight,
                           this->phi, this->dphix,
                           ek, ef);
  }
  else{
    this->fMat->Contribute(this->X, this->jacinv, this->sol, this->dsol, weight, this->axes, this->phi, this->dphix, ek, ef);
  }
}//void

