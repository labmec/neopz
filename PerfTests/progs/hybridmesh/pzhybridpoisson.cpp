//
//  mixedpoisson.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzhybridpoisson.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef PZ_LOG
static TPZLogger logdata("pz.hybridpoisson");
#endif

TPZHybridPoisson::TPZHybridPoisson(): TPZMatPoisson3d(){
}

TPZHybridPoisson::TPZHybridPoisson(int matid, int dim): TPZMatPoisson3d(matid,dim){
  
}

TPZHybridPoisson::~TPZHybridPoisson(){
}


void TPZHybridPoisson::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

/**
 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZHybridPoisson:: ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if (dataleft.fShapeType != TPZMaterialData::EVecandShape || dataright.fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    // loop over the flux functions
    int nvec = dataleft.fVecShapeIndex.size();
    int npressure = dataleft.numberdualfunctions;
    int nlagrange = dataright.phi.Rows();
    for (int ivec=0; ivec<nvec; ivec++) {
        int vec = dataleft.fVecShapeIndex[ivec].first;
        int ish = dataleft.fVecShapeIndex[ivec].second;
        REAL normal = 0.;
        for (int in=0; in<3; in++) {
            normal += data.normal[in]*dataleft.fDeformedDirections(in,vec);
        }
        REAL normalishape = normal * dataleft.phi(ish,0);
        for (int jlagr=0; jlagr < nlagrange; jlagr++) {
            ek(ivec,nvec+npressure+jlagr) += weight*normalishape*dataright.phi(jlagr,0);
            ek(nvec+npressure+jlagr,ivec) += weight*normalishape*dataright.phi(jlagr,0);
        }
    }
    // compute their normal component
    // loop over the pressure functions
    // contribute the Langrange multiplier
    
}

/**
 * @brief It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZHybridPoisson::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    if (dataleft.fShapeType != TPZMaterialData::EVecandShape || dataright.fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    DebugStop();
    
}


