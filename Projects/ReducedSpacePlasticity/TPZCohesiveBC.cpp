/**
 * @file
 * @brief Contains implementations of the TPZCohesiveBC methods.
 */

#include "TPZCohesiveBC.h"


TPZCohesiveBC::TPZCohesiveBC() : TPZMatWithMem<TPZFMatrix<REAL> >(), fSigmaT(0.), fDeltaC(0.), fDeltaT(0.)
{
	 
}

TPZCohesiveBC::TPZCohesiveBC(int id) : TPZMatWithMem<TPZFMatrix<REAL> >(id), fSigmaT(0.), fDeltaC(0.), fDeltaT(0.)
{
	
}

TPZCohesiveBC::TPZCohesiveBC(const TPZCohesiveBC &cp) : TPZMatWithMem<TPZFMatrix<REAL> >(cp), fSigmaT(cp.fSigmaT), fDeltaC(cp.fDeltaC), fDeltaT(cp.fDeltaT)
{
	
}

TPZCohesiveBC::~TPZCohesiveBC()
{
	
}

void TPZCohesiveBC::SetCohesiveData(const REAL &SigmaT, const REAL &DeltaC, const REAL &DeltaT)
{
	fSigmaT = SigmaT;
	fDeltaC = DeltaC;
	fDeltaT = DeltaT;
	TPZFNMatrix<2,REAL> DTST(2,1,0.);
	DTST(0,0) = DeltaT;
	DTST(1,0) = SigmaT;
	this->SetDefaultMem(DTST);
}

void TPZCohesiveBC::CalculateSigma(REAL &w,REAL &DeltaT, REAL &SigmaT, REAL &sigma) const 
{

	
	// Calculating the function
	if (w<=0) { // simulates contact between fracture walls in case of closening
		sigma = -w*fSigmaT/fDeltaT;
	}
	else if (w <= DeltaT) { // until it reaches the sigmaT max on the first time. It alters the curve at each time step
		sigma = w * SigmaT/DeltaT;
	}
	else if (w <= fDeltaC){ // sigma folow the linear law specified for the cohesive tension
		sigma = SigmaT * (1. - (w - DeltaT)/(fDeltaC - DeltaT) );
	}
	else { // if it passes the critical oppening, the cohesive tension ceases to exist
		sigma = 0;
	}
}

void TPZCohesiveBC::UpdateCohesiveCurve(TPZMaterialData &data)
{
	TPZManVector<REAL,3> sol_u = data.sol[0];
	REAL w = 2.*sol_u[1];
	REAL sigma = -6378.;
	
	const int index = data.intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 2 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 2 and Cols 1\n";
		DebugStop();
	}
	REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);	
	CalculateSigma(w,DeltaT,SigmaT,sigma);
	
	if (IsZero(sigma)) {
		w = 0.; // in case of fracture oppened, i dont need to store the sigma anymore or it is zero
	}
	DTST(0,0) = w;
	DTST(1,0) = sigma;
	
	//int index = data.intGlobPtIndex;
	this->MemItem(index) = DTST;
}


void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	TPZFMatrix<REAL> &phi = data.phi;
	int phr = phi.Rows();
	
	REAL CohesiveStress;
	TPZManVector<REAL,3> sol_u = data.sol[0];
	const REAL uy = sol_u[1];
	REAL w = 2.*uy;
	
	const int index = data.intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 2 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 2 and Cols 1\n";
		DebugStop();
	}
	REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);	
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress);
	
	for (int in = 0; in < phr; in++) 
	{
		for (int il = 0; il < fNumLoadCases; il++) 
		{
			ef(in,il) += CohesiveStress * phi(0,in) * weight;        // force in x direction
		}
	}
	
}

void TPZCohesiveBC::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	// datavec[0] = elasticity
	// datavec[1] = pressure
	TPZFMatrix<REAL> &phi = datavec[0].phi;
	int phr = phi.Rows();
	REAL CohesiveStress;
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL uy = sol_u[1];
	REAL w = 2.*uy;
	
	const int index = datavec[0].intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 2 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 2 and Cols 1\n";
		DebugStop();
	}
	REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress);
	
	for (int in = 0; in < phr; in++) 
	{
		for (int il = 0; il < fNumLoadCases; il++) 
		{
			ef(in,il) += CohesiveStress * phi(0,in) * weight;        // force in x direction
		}
	}	
}

void TPZCohesiveBC::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
  
}

void TPZCohesiveBC::Print(std::ostream &out)
{
	out << "TPZMatWithMem::Print:" << std::endl;
	TPZMatWithMem<TPZFMatrix<REAL> >::Print(out);
	out << "\nfSigmaT = " << fSigmaT << std::endl;
	out << "fDeltaT = " << fDeltaT << std::endl;
	out << "fDeltaC = " << fDeltaC << std::endl;
	out << std::endl;
}


