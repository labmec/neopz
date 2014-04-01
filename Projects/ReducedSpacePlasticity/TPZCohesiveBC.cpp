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

void TPZCohesiveBC::CalculateSigma(TPZVec<TPZMaterialData> &datavec, REAL &sigma) const 
{
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL uy = sol_u[1];
	const REAL w = 2.*uy;
	
	const int index = datavec[0].intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 2 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 2 and Cols 1\n";
		DebugStop();
	}
	const REAL DeltaT = DTST(0,0);
	const REAL SigmaT = DTST(1,0);
	
	// Calculating the function
	if (w <= DeltaT) {
		sigma = w * SigmaT/DeltaT;
	}
	else {
		sigma = SigmaT * (1. - (w - DeltaT)/(fDeltaC - DeltaT) );
	}

}

void TPZCohesiveBC::UpdateCohesiveCurve(TPZVec<TPZMaterialData> &datavec)
{
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL w = 2.*sol_u[1];
	REAL sigma = -6378.;
	CalculateSigma(datavec,sigma);
	
	TPZFMatrix<REAL> DTST(2,1,0.);
	DTST(0,0) = w;
	DTST(1,0) = sigma;
	
	int index = datavec[0].intGlobPtIndex;
	this->MemItem(index) = DTST;
}

void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	/*
	TPZFMatrix<REAL> &phi = data[0].phi; // 0 is the elastic mesh
	TPZFMatrix<REAL> &sol = data[0].sol[0];

	
	for(in=0 ; in<numnod ; ++in){
		for(idf = 0;idf<r;idf++) {
			(ef)(in*r+idf,0) += (STATE)phi(in,0)*bc.Val2()(idf,0)*(STATE)weight;
		}
	}	
	 */
}

void TPZCohesiveBC::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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


