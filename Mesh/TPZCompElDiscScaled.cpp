//
//  TPZCompElDiscScaled.cpp
//  MixedElasticityGirkmann
//
//  Created by Philippe Devloo on 03/05/18.
//

#include "pzlog.h"
#include "TPZCompElDiscScaled.h"
#include <math.h>


/** @brief Constructor of the discontinuous element associated with geometric element */
TPZCompElDiscScaled::TPZCompElDiscScaled(TPZCompMesh &mesh,TPZGeoEl *ref) : TPZCompElDisc(mesh,ref){
    //ComputeScale();
}



void TPZCompElDiscScaled::ComputeScale() {
    TPZGeoEl * ref = this->Reference();
    fScale = ref->CharacteristicSize();
}

/** @brief Compute shape functions multiplied by scale value*/
void TPZCompElDiscScaled::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
                          TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
                          REAL &detjac, TPZFMatrix<REAL> &jacinv,
                          TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphix){

    TPZCompElDisc::ComputeShape(intpoint,X,jacobian,axes,detjac,jacinv,phi,dphi,dphix);
    phi*=fScale;
    dphi*=fScale;
    dphix*=fScale;
}

TPZCompElDiscScaled::~TPZCompElDiscScaled()
{
}



