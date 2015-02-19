/*
 *  TPZAxiSymmetricDarcyFlow.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZAxiSymmetricDarcyFlow.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow() : TPZMaterial()
{
    fReservoirdata=NULL;
}

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(int matid) : TPZMaterial(matid)
{
	fReservoirdata=NULL;
}


TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat) : TPZMaterial(mat)
{
	fReservoirdata = mat.fReservoirdata;
}

TPZAxiSymmetricDarcyFlow::~TPZAxiSymmetricDarcyFlow()
{
	
}

void TPZAxiSymmetricDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZAxiSymmetricDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZAxiSymmetricDarcyFlow::Print(std::ostream &out) {
	out << "\t Base class print:\n";
	out << " name of material : " << this->Name() << "\n";
	TPZMaterial::Print(out);
}

int TPZAxiSymmetricDarcyFlow::VariableIndex(const std::string &name) {
	if (!strcmp("Pressure", name.c_str())) return 0;
	if (!strcmp("Velocity", name.c_str())) return 1;
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

int TPZAxiSymmetricDarcyFlow::NSolutionVariables(int var) {
	switch(var) {
		case 0:
			return 1; // Scalar
		case 1:
			return 3; // Vector
		default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
	}
    return 0;
}

void TPZAxiSymmetricDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
	
//	TPZFMatrix<STATE> &dQsol = datavec[0].dsol[0];
//	TPZFMatrix<STATE> &dPsol = datavec[1].dsol[0];

    TPZVec<REAL> Q = datavec[0].sol[0];
    TPZVec<REAL> P = datavec[1].sol[0];
	
	switch(var) {
		case 0:
		{
			Solout[0] = P[0];
		}
			break;
		case 1:
		{
            Solout[0] = Q[0];
            Solout[1] = Q[1];
		}
			break;
		default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
	}
}

void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	
    // Here implements Volumetric integrals
    
    

    

}

void TPZAxiSymmetricDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	
	
}

int TPZAxiSymmetricDarcyFlow::ClassId() const {
	return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Write(TPZStream &buf, int withclassid) {

    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fReservoirdata->fPref);
    buf.Write(&fReservoirdata->fKab(0,0));
	
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fReservoirdata->fPref);
    buf.Read(&fReservoirdata->fKab(0,0));
	
}
