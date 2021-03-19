/**
 * @file
 * @brief Contains implementations of the TPZCohesiveBC methods.
 */

#include "TPZCohesiveBC.h"


#ifdef PZ_LOG
#include "pzlog.h"
static PZLogger logger("pz.reducedspace.data");
#endif

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
#ifdef PZDEBUG
  if (w>=0 && propageted > 1.) { // if it is already propageted sigma is always zero, but i already do this check in the contribute
    DebugStop(); // i hope it never gets here
  }
#endif
  REAL big = 1.;
	// Calculating the function
	if (w < 0.) { // simulates contact between fracture walls in case of closening
		sigma = big * w*fSigmaT/fDeltaT;
	}
	else if (w <= DeltaT) { // until it reaches the sigmaT max on the first time. It alters the curve at each time step
		sigma = big * w * SigmaT/DeltaT;
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
 
#ifdef PZDEBUG
  if (w>=0 && propageted > 1.) { // if it is already propageted deriv is always zero, but i already do this check in the contribute
    DebugStop();
  }
#endif
	
  REAL big = 1.;
  // Calculating the deriv of the function
	if (w < 0.) { // simulates compression between fracture walls in case of closening
		deriv = big * fSigmaT/fDeltaT;
	}
	else if (w <= DeltaT) { // until it reaches the sigmaT max on the first time. It alters the curve at each time step
		deriv = big * SigmaT/DeltaT;
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
    globFractInputData.GetnElPropag()++;
    globFractInputData.SetPropagated();
  }
  else if (w < DeltaT){
    return;
  }
  else if (w == DeltaT){
    w+=fDeltaC*0.001; //AQUINATHAN
    CalculateSigma(w,DeltaT,SigmaT,sigma, propageted);
  }
  else if (w > DeltaT && w <= fDeltaC){
    CalculateSigma(w,DeltaT,SigmaT,sigma, propageted);
  }
  else{
    DebugStop(); // should be one of the above cases
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
    return;
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

	REAL InfluenceArea = globFractInputData.Lmax_edge();
	
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress,propageted);
	CohesiveStress *= InfluenceArea;
	this->CalculateCohesiveDerivative(w,DeltaT,SigmaT,DerivCohesive,propageted);
	DerivCohesive *= InfluenceArea;
	
	REAL big = TPZMaterial::gBigNumber;
	REAL dif = sol_u[0];
	for (int in = 0; in < phr; in++) 
	{
		for (int il = 0; il < fNumLoadCases; il++) 
		{
			ef(2*in,il) += - big * dif * phi(in,0) * weight; // para impedir deslocamento em x
			ef(2*in+1,il) += CohesiveStress * phi(in,0) * weight; // negativo da formula do residuo com negativo do residuo por causa do nonlinanalysis da positivo
		}
		for (int jn = 0; jn < phr; jn++) {
		  ek(2*in,2*jn)     += big * phi(in,0) *phi(jn,0) * weight;
			ek(2*in+1,2*jn+1) += (- 2. * DerivCohesive) * phi(jn,0) * phi(in,0) * weight; // eh o negativo da formula do residuo
		}
	}
	/*
	 REAL dif0 = (data.sol[0][0]-v2[0]);
	 REAL dif1 = (data.sol[0][1]-v2[1]);
	 switch (bc.Type()) {
	 case 0 :			// Dirichlet condition
	 {
	 for(in = 0 ; in < phr; in++) {
	 ef(2*in,0)   += - BIGNUMBER * dif0 * phi(in,0) * weight;        // forced v2 displacement
	 ef(2*in+1,0) += - BIGNUMBER * dif1 * phi(in,0) * weight;        // forced v2 displacement
	 for (jn = 0 ; jn < phi.Rows(); jn++)
	 {
	 ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
	 ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
	 }
	 }
	 }	 
	 */
}

void TPZCohesiveBC::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	// datavec[0] = elasticity
	// datavec[1] = pressure

	if(gState == ELastState)
	{
		return;
	}	
	
  TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
	if(shapetype != datavec[0].EVecShape)
	{
    this->Contribute(datavec[0], weight, ek, ef);
    return;
  }
  
  if(fUpdateMem){
    UpdateCohesiveCurve(datavec[0]);
    return;
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
	
	REAL InfluenceArea = globFractInputData.Lmax_edge();
	
	this->CalculateSigma(w,DeltaT,SigmaT,CohesiveStress,propageted);
	CohesiveStress *= InfluenceArea;
	this->CalculateCohesiveDerivative(w,DeltaT,SigmaT,DerivCohesive,propageted);
	DerivCohesive *= InfluenceArea;

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
	
#ifdef PZ_LOG
	if (logger.isDebugEnabled()) {
		std::stringstream str;
		str << "\n------- Contribute do Cohesive -------" << std::endl;
		str << "GeoElId = " << datavec[0].gelElId << std::endl; 
		ek.Print("ek",str);
		ef.Print("ef",str);
		LOGPZ_DEBUG(logger,str.str())		
	}
#endif
		
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


