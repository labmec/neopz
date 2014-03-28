/**
 * @file
 * @brief Contains implementations of the TPZCohesiveBC methods.
 */

#include "TPZCohesiveBC.h"


TPZCohesiveBC::TPZCohesiveBC() : TPZMaterial(), fSigmaT(0.), fDeltaC(0.)
{
	
}

TPZCohesiveBC::TPZCohesiveBC(int id) : TPZMaterial(id), fSigmaT(0.), fDeltaC(0.)
{
	
}

TPZCohesiveBC::~TPZCohesiveBC()
{
	
}

void TPZCohesiveBC::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	
}

void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	
}

void TPZCohesiveBC::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	
}

void TPZCohesiveBC::Print(std::ostream &out)
{
	
}

void TPZCohesiveBC::Read(TPZStream &buf, void *context)
{
	PZError << "METHOD NOT IMPLEMENTED!!!\n";
	DebugStop();
}

void TPZCohesiveBC::Write(TPZStream &buf, int withclassid)
{
	PZError << "METHOD NOT IMPLEMENTED!!!\n";
	DebugStop();	
}

int TPZCohesiveBC::VariableIndex(const std::string &name)
{
	
}

int TPZCohesiveBC::NSolutionVariables(int var)
{
	
}


void TPZCohesiveBC::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
	
}

