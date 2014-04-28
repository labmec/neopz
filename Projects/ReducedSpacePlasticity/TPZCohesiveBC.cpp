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

void TPZCohesiveBC::CalculateSigma(TPZMaterialData &data, REAL &sigma) const
{
	TPZManVector<REAL,3> sol_u = data.sol[0];
	const REAL uy = sol_u[1];
	const REAL w = 2.*uy;
	
	const int index = data.intGlobPtIndex;
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 2 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 2 and Cols 1\n";
		DebugStop();
	}
	const REAL DeltaT = DTST(0,0);
	const REAL SigmaT = DTST(1,0);
	
	// Calculating the function
  if (w < 0) {
    sigma =  - w * SigmaT/DeltaT;
  }
	else if (w <= DeltaT) {
		sigma = w * SigmaT/DeltaT;
	}
	else if (w <= fDeltaC) {
		sigma = SigmaT * (1. - (w - DeltaT)/(fDeltaC - DeltaT) );
	}
  else{
    sigma = 0.;
  }

}

void TPZCohesiveBC::UpdateCohesiveCurve(TPZMaterialData &data)
{
	TPZManVector<REAL,3> sol_u = data.sol[0];
	const REAL w = 2.*sol_u[1];
	REAL sigma = -6378.;
	CalculateSigma(data,sigma);
	
	TPZFMatrix<REAL> DTST(2,1,0.);
	DTST(0,0) = w;
	DTST(1,0) = sigma;
	
	int index = data.intGlobPtIndex;
	this->MemItem(index) = DTST;
}

void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	TPZFMatrix<REAL> &phi = data.phi; // 0 is the elastic mesh
  REAL sigma = -6378.;
	CalculateSigma(data,sigma);
  
  
  int in = 0, phr = phi.Rows();
  for (in = 0; in < phr; in++)
  {
    for (int il = 0; il <fNumLoadCases; il++)
    {
      TPZFNMatrix<2,STATE> v2(2,1);
      v2(1,0) = -sigma;
      ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
      ef(2*in+1,il) +=  v2(1,0) * phi(in,0) * weight;      // force in y direction
    }
  }
  
}

/*

void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef)
{
	
	TPZFMatrix<REAL> &phi = data.phi; // 0 is the elastic mesh
	
  int in = 0, phr = phi.Rows();
  for (in = 0; in < phr; in++)
  {
    for (int il = 0; il <fNumLoadCases; il++)
    {
      TPZFNMatrix<2,STATE> v2(2,1);
      v2(1,0) = -1;
      ef(2*in,il) += v2(0,0) * phi(in,0) * weight;        // force in x direction
      ef(2*in+1,il) +=  v2(1,0) * phi(in,0) * weight;      // force in y direction
    }
  }
}

*/

void TPZCohesiveBC::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	
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


