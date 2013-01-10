//
//  pznewl2projection.cpp
//  PZ
//
//  Created by Agnaldo Farias on 1/9/13.
//
//

#include "pznewl2projection.h"


TPZL2ProjectionForGradient::TPZL2ProjectionForGradient(int id, int dim, int nstate): TPZDiscontinuousGalerkin(id)
{
	fgradients.Redim(0,0);
    this->fDim = dim;
    this->fNStateVars = nstate;
}


TPZL2ProjectionForGradient::~TPZL2ProjectionForGradient()
{
}