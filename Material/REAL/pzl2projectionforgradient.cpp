//
//  pznewl2projection.cpp
//  PZ
//
//  Created by Agnaldo Farias on 1/9/13.
//
//

#include "pzl2projectionforgradient.h"


TPZL2ProjectionForGradient::TPZL2ProjectionForGradient(int matid, int dim, int nstate): TPZDiscontinuousGalerkin(matid)
{
	fgradients.Redim(0,0);
    this->fDim = dim;
    this->fNStateVars = nstate;
    this->fmatId = matid;
}


TPZL2ProjectionForGradient::~TPZL2ProjectionForGradient()
{
}


void TPZL2ProjectionForGradient::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
}

void TPZL2ProjectionForGradient::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){

    //Falta preencher
    DebugStop();
}

void TPZL2ProjectionForGradient::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
    
    //Falta preencher
    DebugStop();
}

int TPZL2ProjectionForGradient::VariableIndex(const std::string &name){
    
    if(!strcmp("Solution",name.c_str()))        return  1;
    return TPZMaterial::VariableIndex(name);
}

int TPZL2ProjectionForGradient::NSolutionVariables(int var){
    if(var == 1) return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZL2ProjectionForGradient::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
    
    //Falta preencher
    DebugStop();
}

