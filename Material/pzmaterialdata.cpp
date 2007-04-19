//$Id: pzmaterialdata.cpp,v 1.6 2007-04-19 20:19:55 tiago Exp $

#include "pzmaterialdata.h"
#include "pzmaterial.h"
#include "pzcompel.h"
#include "pzelmat.h"
#include <sstream>
#include "pzerror.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"



TPZMaterialData::TPZMaterialData(){
  this->SetAllRequirements(false);
}

TPZMaterialData::TPZMaterialData( const TPZMaterialData &cp ){
  this->operator =(cp);
}

TPZMaterialData & TPZMaterialData::operator= (const TPZMaterialData &cp ){
  this->fNeedsSol = cp.fNeedsSol;
  this->fNeedsNeighborSol = cp.fNeedsNeighborSol;
  this->fNeedsHSize = cp.fNeedsHSize;
  this->phi = cp.phi;
  this-> phil = cp.phil;
  this->phir = cp.phir;
  this->dphix = cp.dphix;
  this->dphixl = cp.dphixl;
  this->dphixr = cp.dphixr;
  this->axes = cp.axes;
  this->axesleft = cp.axesleft;
  this->axesright = cp.axesright;
  this->jacobian = cp.jacobian;
  this->leftjac = cp.leftjac;
  this->rightjac = cp.rightjac;
  this->jacinv = cp.jacinv;
  this->leftjacinv = cp.leftjacinv;
  this->rightjacinv = cp.rightjacinv;
  this->normal = cp.normal;
  this->x = cp.x;
  this->p = cp.p;
  this->leftp = cp.leftp;
  this->rightp = cp.rightp;
  this->sol = cp.sol;
  this->soll = cp.soll;
  this->solr = cp.solr;
  this->dsol = cp.dsol;
  this->dsoll = cp.dsoll;
  this->dsolr = cp.dsolr;
  this->HSize = cp.HSize;
  this->detjac = cp.detjac;
  this->leftdetjac = cp.leftdetjac;
  this->rightdetjac = cp.rightdetjac;
  return *this;
}

TPZMaterialData::~TPZMaterialData(){
//NOTHING TO BE DONE!
}

void TPZMaterialData::SetAllRequirements(bool set){
  this->fNeedsSol = set;
  this->fNeedsNeighborSol = set;
  this->fNeedsHSize = set;
}

void TPZMaterialData::InvertLeftRightData(){
  TPZMaterialData cp(*this);
  this->leftdetjac = cp.rightdetjac;
  this->leftjac = cp.rightjac;
  this->leftjacinv = cp.rightjacinv;
  this->leftp = cp.rightp; 
  this->phil = cp.phir;
  this->dphixl = cp.dphixr;
  this->axesleft = cp.axesright;
  this->soll = cp.solr;
  this->dsoll = cp.dsolr;
  
  this->rightdetjac = cp.leftdetjac;
  this->rightjac = cp.leftjac;
  this->rightjacinv = cp.leftjacinv;
  this->rightp = cp.leftp;
  this->phir = cp.phil;
  this->dphixr = cp.dphixl;
  this->axesright = cp.axesleft;
  this->solr = cp.soll;
  this->dsolr = cp.dsoll;  
  
  const int n = this->normal.NElements();
  for(int i = 0; i < n; i++){
    this->normal[i] *= -1.;
  }
}








