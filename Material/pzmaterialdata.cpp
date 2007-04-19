//$Id: pzmaterialdata.cpp,v 1.3 2007-04-19 19:01:26 tiago Exp $

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

TPZMaterialData & TPZMaterialData::operator= (const TPZMaterialData &A ){
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
  TPZMaterialData cp(data);
  this->leftdetjac = cp.rightdetjac;
  this->leftjac = cp.rightjac;
  this->leftjacinv = cp.rightjacinv;
  this->leftp = cp.rightp;


  this->rightdetjac = cp.leftdetjac;
  this->rightjac = cp.leftjac;
  this->rightjacinv = cp.leftjacinv;
  this->rightp = cp.leftp;
}








