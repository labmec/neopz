/**
 * @file
 * @brief Contains implementations of the TPZCohesiveBC methods.
 */

#include "TPZCohesiveBC.h"


TPZCohesiveBC::TPZCohesiveBC() : TPZMatWithMem<TPZFMatrix<REAL> >(), fSigmaT(0.), fDeltaC(0.), fDeltaT(0.)
{
 	this->SetCurrentState();
}

TPZCohesiveBC::TPZCohesiveBC(int id) : TPZMatWithMem<TPZFMatrix<REAL> >(id), fSigmaT(0.), fDeltaC(0.), fDeltaT(0.)
{
	this->SetCurrentState();
}

TPZCohesiveBC::TPZCohesiveBC(const TPZCohesiveBC &cp) : TPZMatWithMem<TPZFMatrix<REAL> >(cp), fSigmaT(cp.fSigmaT), fDeltaC(cp.fDeltaC), fDeltaT(cp.fDeltaT)
{
	gState = cp.gState;
}


TPZCohesiveBC::~TPZCohesiveBC()
{
	
}

void TPZCohesiveBC::SetCohesiveData(const REAL &SigmaT, const REAL &DeltaC, const REAL &DeltaT)
{
	fSigmaT = SigmaT;
	fDeltaC = DeltaC;
	fDeltaT = DeltaT;
	TPZFNMatrix<2,REAL> DTST(3,1,0.);
	DTST(0,0) = DeltaT;
	DTST(1,0) = SigmaT;
  REAL notpropageted = 0.;
  DTST(2,0) = notpropageted;
	this->SetDefaultMem(DTST);
}

void TPZCohesiveBC::CalculateSigma(REAL &w,REAL &DeltaT, REAL &SigmaT, REAL &sigma, REAL &propageted) const
{
#ifdef DEBUG
  if (w>=0 && propageted > 1.) { // if it is already propageted sigma is always zero, but i already do this check in the contribute
    DebugStop(); // i hope it never gets here
  }
#endif
  
	// Calculating the function
	if (w < 0.) { // simulates contact between fracture walls in case of closening
		REAL big = 1.;
		sigma = big*w*fSigmaT/fDeltaT;
	}
	else if (w <= DeltaT) { // until it reaches the sigmaT max on the first time. It alters the curve at each time step
		sigma = w * SigmaT/DeltaT;
	}
	else if (w <= fDeltaC){ // sigma folow the linear law specified for the cohesive tension
		sigma = SigmaT * (1. - (w - DeltaT)/(fDeltaC - DeltaT) );
	}
	else { // if it passes the critical oppening, the cohesive tension ceases to exist
		sigma = 0.;
	}
  
	sigma *= -1.0; // the cohesive stress in oposite to the displacement

}

void TPZCohesiveBC::CalculateCohesiveDerivative(REAL &w,REAL &DeltaT, REAL &SigmaT, REAL &deriv, REAL &propageted) const
{
 
#ifdef DEBUG
  if (w>=0 && propageted > 1.) { // if it is already propageted deriv is always zero, but i already do this check in the contribute
    DebugStop();
  }
#endif
	
  // Calculating the deriv of the function
	if (w < 0.) { // simulates compression between fracture walls in case of closening
		REAL big = 1.;
		deriv = big*fSigmaT/fDeltaT;
	}
	else if (w <= DeltaT) { // until it reaches the sigmaT max on the first time. It alters the curve at each time step
		deriv = SigmaT/DeltaT;
	}
	else if (w <= fDeltaC){ // sigma folow the linear law specified for the cohesive tension
		deriv = -SigmaT/(fDeltaC-DeltaT);
	}
	else { // if it passes the critical oppening, the cohesive tension ceases to exist
		deriv = 0;
	}
		deriv *= -1.0;
}

void TPZCohesiveBC::UpdateCohesiveCurve(TPZMaterialData &data)
{
	TPZManVector<REAL,3> sol_u = data.sol[0];
	REAL w = 2.*sol_u[1];
	REAL sigma = -6378.;
	
	const int index = data.intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	
	if (DTST.Rows() != 3 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 3 and Cols 1\n";
		DebugStop();
	}
	REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);
  REAL propageted = DTST(2,0);
  
  if (propageted > 1.) { // if propageted, there is nothing to be done
    return;
  }
  
  if (w >= fDeltaC) {
    w = fDeltaC;
    sigma = 0.;
    propageted = 2.;
  }
  else{
    CalculateSigma(w,DeltaT,SigmaT,sigma, propageted);
  }
  
	DTST(0,0) = w;
	DTST(1,0) = -sigma; // I assume it isnt beatiful, but the value of the sigma in function must be positive according to the way that i tought
  DTST(2,0) = propageted;
	
	//int index = data.intGlobPtIndex;
	this->MemItem(index) = DTST;
}


void TPZCohesiveBC::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	
	if(gState == ELastState)
	{
		return;
	}
	
  if(fUpdateMem){
    UpdateCohesiveCurve(data);
  }
  
	TPZFMatrix<REAL> &phi = data.phi;
	int phr = phi.Rows();
	
	REAL CohesiveStress, DerivCohesive;
	TPZManVector<REAL,3> sol_u = data.sol[0];
	const REAL uy = sol_u[1];
	REAL w = 2.*uy;
	
	const int index = data.intGlobPtIndex; 
	if (index == -1) {
		return;
	}
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 3 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 3 and Cols 1\n";
		DebugStop();
	}
	
	
	REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);
  REAL propageted = DTST(2,0); // Fracture is considered oppened when propageted is greater than 1!
  if (propageted > 1.) { // if propageted, there is nothing to be done
    return;
  }
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress,propageted);
	this->CalculateCohesiveDerivative(w,DeltaT,SigmaT,DerivCohesive,propageted);
		
	
	for (int in = 0; in < phr; in++) 
	{
		for (int il = 0; il < fNumLoadCases; il++) 
		{
			ef(2*in+1,il) += CohesiveStress * phi(in,0) * weight; // negativo da formula do residuo com negativo do residuo por causa do nonlinanalysis da positivo
		}
		for (int jn = 0; jn < phr; jn++) {
			ek(2*in+1,2*jn+1) += (- 2. * DerivCohesive) * phi(jn,0) * phi(in,0) * weight; // eh o negativo da formula do residuo
		}
	}
	
}

void TPZCohesiveBC::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	// datavec[0] = elasticity
	// datavec[1] = pressure

	if(gState == ELastState)
	{
		return;
	}	
	
  if(fUpdateMem){
    UpdateCohesiveCurve(datavec[0]);
  }
	TPZFMatrix<REAL> &phi = datavec[0].phi;
	int phc = phi.Cols();
	REAL CohesiveStress, DerivCohesive;
	TPZManVector<REAL,3> sol_u = datavec[0].sol[0];
	const REAL uy = sol_u[1];
	REAL w = 2.*uy;
	
	const int index = datavec[0].intGlobPtIndex; 
	TPZFMatrix<REAL> DTST = this->MemItem(index);
	if (DTST.Rows() != 3 || DTST.Cols() != 1) {
		PZError << "The DTST memory number of Rows must be 3 and Cols 1\n";
		DebugStop();
	}

  REAL DeltaT = DTST(0,0);
	REAL SigmaT = DTST(1,0);
  REAL propageted = DTST(2,0); // Fracture is considered oppened when propageted is greater than 1!
  if (propageted > 1.) { // if propageted, there is nothing to be done
    return;
  }
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress,propageted);
	this->CalculateCohesiveDerivative(w,DeltaT,SigmaT,DerivCohesive,propageted);
	
	/*
	// PARA COLOCAR COMO SE FOSSE UM DIRICHLET (TESTE)
	const REAL big  = TPZMaterial::gBigNumber;
	
	for(int in = 0 ; in < phc; in++)
	{
		for(int il = 0; il < this->fNumLoadCases; il++)
		{
			//termo big*u*v do vetor de carga
			TPZFNMatrix<3,STATE> v2(2,1,0.);
			ef(in,il) += big * ( v2(0,il)*phi(0,in) + v2(1,il)*phi(1,in) ) * weight;
			
			//termo big*u*v da matriz
			ef(in,il) += (-1.)*big * ( phi(0,in)*sol_u[0] + phi(1,in)*sol_u[1] ) * weight;
		}
		
		for (int jn = 0; jn < phc; jn++)
		{
			ek(in,jn) += big * ( phi(0,in)*phi(0,jn) + phi(1,in)*phi(1,jn) ) * weight;
		}
	}
	*/
	for (int in = 0; in < phc; in++)
	{
		for (int il = 0; il < fNumLoadCases; il++)
		{
			ef(in,il) += CohesiveStress * phi(1,in) * weight; // negativo da formula do residuo com negativo do residuo por causa do nonlinanalysis da positivo
		}
		for (int jn = 0; jn < phc; jn++) {
			ek(in,jn) += (- 2. * DerivCohesive) * phi(1,jn) * phi(1,in) * weight; // eh o negativo da formula do residuo
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


