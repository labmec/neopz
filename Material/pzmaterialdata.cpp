//$Id: pzmaterialdata.cpp,v 1.2 2007-04-12 20:01:13 tiago Exp $

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

TPZMaterialData::~TPZMaterialData(){
//NOTHING TO BE DONE!
}

void TPZMaterialData::SetAllRequirements(bool set){
  this->fNeedsSol = set;
  this->fNeedsNeighborSol = set;
  this->fNeedsHSize = set;
}

