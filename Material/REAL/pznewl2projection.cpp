//
//  pznewl2projection.cpp
//  PZ
//
//  Created by Agnaldo Farias on 1/9/13.
//
//

#include "pznewl2projection.h"


TPZNewL2Projection::TPZNewL2Projection(int id, int dim, int nstate, TPZFMatrix<REAL> gradients): TPZDiscontinuousGalerkin(id)
{
    this->fDim = dim;
    this->fNStateVars = nstate;
    this->fgradients = gradients;
}

TPZNewL2Projection::TPZNewL2Projection(const TPZNewL2Projection &cp): TPZDiscontinuousGalerkin(cp)
{
    this->fDim = cp.fDim;
    this->fNStateVars = cp.fNStateVars;
    this->fgradients = cp.fgradients;
}

TPZNewL2Projection::~TPZNewL2Projection()
{
}