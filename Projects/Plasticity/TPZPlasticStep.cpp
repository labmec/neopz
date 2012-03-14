// $Id: TPZPlasticStep.cpp,v 1.60 2010-11-03 18:21:36 diogo Exp $

#include "TPZPlasticStep.h"
#include "TPZElasticResponse.h"
#include "TPZLadeNelsonElasticResponse.h"
#include "TPZThermoForceA.h"
#include "TPZYCVonMises.h"
#include "TPZYCTresca.h"
#include "tpzyctrescaregularized.h"
#include "TPZYCLadeKim.h"
#include "TPZLadeKimThermoForceA.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZSandlerDimaggio.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZSandlerDimaggioThermoForceA.h"
#include "tpzycvonmisescombtresca.h"
#include "pzfmatrix.h"
#include "pzstepsolver.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"
#include "TPZYCDruckerPrager.h"
#include "TPZYCRankine.h"
#include <fenv.h>//NAN DETECTOR

TPZPlasticIntegrMem<REAL, 4> teste;

using namespace std;

#include "tfad.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr plasticIntegrLogger(Logger::getLogger("plasticity.plasticIntegr"));
#endif

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr logger(Logger::getLogger("PLASTIC_STEP"));
static LoggerPtr loggerPlasticResidual(Logger::getLogger("PLASTIC_RESIDUAL"));
static LoggerPtr loggerPlasticLoop(Logger::getLogger("loggerPlasticLoop"));
#endif

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::IntegrationSteps() const
{
	return fPlasticMem.NElements() -2;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(
				const TPZPlasticState<REAL> & state)
{
	fN = state;
	
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout1, sout2;
    sout1 << ">>> SetState_Internal ***";
    sout2 << "\nfN.fEpsP << " << fN.fEpsP << "\nfN.fAlpha << " << fN.fAlpha;
    LOGPZ_INFO(logger,sout1.str().c_str());
    LOGPZ_DEBUG(logger,sout2.str().c_str());
  }
#endif
}

template <class YC_t, class TF_t, class ER_t>
const TPZPlasticState<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::GetState_Internal() const
{
	return fN;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetState(
				const TPZPlasticState<REAL> & state)
{
	int multipl = SignCorrection();
	fN = state;
	fN.fEpsP *= multipl;
	fN.fEpsT *= multipl;
	
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout1, sout2;
    sout1 << ">>> SetUp ***";
    sout2 << "\nfN.fEpsP << " << fN.fEpsP << "\nfN.fAlpha << " << fN.fAlpha;
    LOGPZ_INFO(logger,sout1.str().c_str());
    LOGPZ_DEBUG(logger,sout2.str().c_str());
  }
#endif
}

template <class YC_t, class TF_t, class ER_t>
const TPZPlasticState<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::GetState() const
{
	int multipl = SignCorrection();
	TPZPlasticState<REAL> N(fN);
	N.fEpsP *= multipl;
	N.fEpsT *= multipl;
	
	return N;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetUp(
				const TPZTensor<REAL> & epsTotal)
{
	fN.fEpsT = epsTotal;

#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout1, sout2;
    sout1 << ">>> SetUp ***";
    LOGPZ_INFO(logger,sout1.str().c_str());
  }
#endif
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{
	TPZTensor<REAL> epsTotal_Internal(epsTotal);
	epsTotal_Internal *= SignCorrection();
	ApplyStrain_Internal(epsTotal_Internal);
}
	
template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
	int multipl = SignCorrection();
	TPZTensor<REAL> epsTotal_Internal(epsTotal);
	epsTotal_Internal *= multipl;
	ApplyStrainComputeSigma_Internal(epsTotal_Internal, sigma);
	sigma *= multipl;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix &Dep)
{
	int multipl = SignCorrection();
	TPZTensor<REAL> epsTotal_Internal(epsTotal);
	epsTotal_Internal *= multipl;
	ApplyStrainComputeDep_Internal(epsTotal_Internal, sigma, Dep);
	sigma *= multipl;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal)
{
	TPZTensor<REAL> sigma_Internal(sigma);
	sigma_Internal *= SignCorrection();
	ApplyLoad_Internal(sigma_Internal, epsTotal);
	epsTotal *= SignCorrection();
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const
{
	TPZTensor<REAL> epsTotal_Internal(epsTotal);
	epsTotal_Internal *= SignCorrection();

	Phi_Internal(epsTotal_Internal, phi);
	
    return;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetTensionSign(int s)
{
	fInterfaceTensionSign = (s >= 0) ? 1 : -1;
}

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::SignCorrection() const
{
	return fMaterialTensionSign * fInterfaceTensionSign;
}

template <class YC_t, class TF_t, class ER_t>
template<class T>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ComputePlasticVars(const TPZPlasticState<T> & state_T,
														  TPZTensor<T> & sigma_T,
														  T & A_T)const
{
  // variable representing the stress and elastic deformation
  TPZTensor<T> epsE_T(state_T.EpsT() );
  // subtract the value of the current plastic deformation
  epsE_T.Add(state_T.EpsP(), T(-1.));
  // compute the stress of the elastic response
  fER.Compute(epsE_T, sigma_T);
  // compute the value of the thermo dynamical force for the given damage variable
  A_T = fTFA.Compute(state_T.Alpha() );	
}
														

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::IsStrainElastic(const TPZPlasticState<REAL> &state)const
{
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout1, sout2;
    sout1 << ">>> IsStrainElastic ***";
    sout2 << "\nstate.EpsT() << " << state.EpsT();
    LOGPZ_INFO(logger,sout1.str().c_str());
    LOGPZ_DEBUG(logger,sout2.str().c_str());
  }
#endif
	
  TPZTensor<REAL> sigma;
  REAL A;
	
  ComputePlasticVars<REAL>(state, sigma, A);
	
  // compute the value of the yield functions
  TPZManVector<REAL,10> phi(YC_t::NYield);

  fYC.Compute(sigma, A, phi, 0);
  // verify if any yield function indicates plastification
  int i;
  for(i=0; i< YC_t::NYield; i++)
  {
    if(phi[i] > 0.) 
	{
		break;
	}
  }

	
  // if we are in the elastic range
  if(i == YC_t::NYield)
  {
	#ifdef LOG4CXX_PLASTICITY
	{
	std::stringstream sout;
	sout << "*** IsStrainElastic *** Strain yet in the elastic range - no damage variable needs update.\nExiting method ApplyStrain."
		 << "\n Phi = ";
	for(int j = 0; j < YC_t::NYield; j++)sout << phi[j] << "  ";
	LOGPZ_INFO(logger,sout.str().c_str());
	}
	#endif
    return 1;
  }
#ifdef LOG4CXX_PLASTICITY
{
    std::stringstream sout;
    sout << "*** IsStrainElastic *** Strain exceeds the elastic range - damage variables need update"
		 << "\n Phi = ";
	for(int j = 0; j < YC_t::NYield; j++)sout << phi[j] << "  ";
    LOGPZ_INFO(logger,sout.str().c_str());
}
#endif	
   return 0;
	
}


template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrain_Internal(const TPZTensor<REAL> &epsTotal)
{

#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout  << ">>> ApplyStrain_Internal *** Imposed epsTotal << " << epsTotal;
		LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
    }
#endif
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout  << ">>> ApplyStrain_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif
	
    
    ProcessStrainNoSubIncrement(epsTotal);
//	ProcessStrain(epsTotal);
	
	int n = fPlasticMem.NElements();
	
    // load the integrated values as the current state
	TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(fPlasticMem[n-1].fPlasticState);
	
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout << "*** ProcessStrain *** Exiting Method.";
        LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif
	
}


template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessStrain(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep)
{

#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrain ***";
		sout1 << ">>> EElastoPlastic ***" << ep;
        LOGPZ_INFO(logger,sout1.str().c_str());
    }
#endif
	
//#define _XX_ 0
//#define _XY_ 1
//#define _XZ_ 2
//#define _YY_ 3
//#define _YZ_ 4
//#define _ZZ_ 5
	

	
//	
    REAL yieldMultipl = 0;
	
    fPlasticMem.Resize(0);
    PushPlasticMem(fN, 
                   0. /*k*/,
                   0. /*lambda*/,
                   TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,
                   TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,
				   0 /*forceYield*/);
	
    TPZPlasticState<REAL> stateAtYield(fN),
	                      Np1(fN); // Np1 state with fN guesses
    Np1.fEpsT = epsTotal;
	
	bool elastic = true;
	
	if( ep == EForceElastic)
	{
	#ifdef LOG4CXX_PLASTICITY
	    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrain *** behaviour imposed to be Elastic";
        LOGPZ_INFO(logger,sout1.str().c_str());
	    }
	#endif
		elastic = true;
	}
	
	if( ep == EForcePlastic)
	{
	#ifdef LOG4CXX_PLASTICITY
	    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrain *** behaviour imposed to be Plastic";
        LOGPZ_INFO(logger,sout1.str().c_str());
	    }
	#endif
		elastic = false;
	}
	
	if( ep == EAuto ) 
	{
		elastic = IsStrainElastic(Np1) == 1;
	}

	/*
	fMaterialElasticOrPlastic=0;
	if(fMaterialElasticOrPlastic==0)
	{
		elastic = true;
		//return;
	}
	*/
	
    if(elastic)
    {
	    PushPlasticMem(Np1,
                       1. /*k*/,
                       0. /*lambda unused - elastic state*/,
                       TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,
                       TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,
					   0 /*forceYield*/);
        return;
    }	
	
	// Plastic Integration needed
	
	yieldMultipl = FindPointAtYield(Np1.EpsT(), stateAtYield);
  
    PushPlasticMem(stateAtYield,
	               yieldMultipl /*k*/,
                   0. /*lambda unused - elastic state*/,
                   TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,
                   TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,
				   0 /*forceYield*/);
	
	REAL multipl = 0.99;
	TPZTensor<REAL> DeltaEpsP_guess = Np1.fEpsT;
	DeltaEpsP_guess.Add(stateAtYield.fEpsT,-1.);
	Np1.fEpsP.Add(DeltaEpsP_guess, multipl);
	
	// we will integrate from stateAtYield .fEpsT to Np1.fEpsT
	// the plastic strain will evoluate from stateAtYield.fAlpha stateAtYield.fEpsP
    PlasticIntegrate(stateAtYield, Np1, fIntegrTol);
	  
}



template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessStrainNoSubIncrement(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep)
{
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout1;
        sout1 << ">>> ProcessStrain ***";
        sout1 << ">>> EElastoPlastic ***" << ep;
        LOGPZ_INFO(logger,sout1.str().c_str());
    }
#endif
    
    //    
    REAL yieldMultipl = 0;
    
    fPlasticMem.Resize(0);
    PushPlasticMem(fN, 0. /*k*/,0. /*lambda*/,TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,0 /*forceYield*/);
    
    TPZPlasticState<REAL> stateAtYield(fN),Np1(fN); // Np1 state with fN guesses
    
    Np1.fEpsT = epsTotal;
    
    bool elastic = true;
    
    if( ep == EForceElastic)
    {
        elastic = true;
    }
    
    if( ep == EForcePlastic)
    {
        elastic = false;
    }
    
    if( ep == EAuto ) 
    {
        elastic = IsStrainElastic(Np1) == 1;
    }
    
    if(elastic)
    {
        PushPlasticMem(Np1,1. /*k*/,0. /*lambda unused - elastic state*/,TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,
                       TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,0 /*forceYield*/);
        return;
    }    
    
    // Plastic Integration needed
    //else
    //{
    yieldMultipl = FindPointAtYield(Np1.EpsT(), stateAtYield);
    
    PushPlasticMem(stateAtYield,yieldMultipl /*k*/,0. /*lambda unused - elastic state*/,TPZManVector<REAL, YC_t::NYield>(YC_t::NYield,0)/*delGamma*/,
                   TPZManVector<int, YC_t::NYield>(YC_t::NYield,0)/*validEqs*/,0 /*forceYield*/);
    
    REAL multipl = 0.99;
    TPZTensor<REAL> DeltaEpsP_guess = Np1.fEpsT;
    DeltaEpsP_guess.Add(stateAtYield.fEpsT,-1.);
    Np1.fEpsP.Add(DeltaEpsP_guess, multipl);
    
    ///*********///////**********//////////****************************************************************************************************************
    
    
    int succeeded;
    
    TPZManVector<REAL,YC_t::NYield> delGamma(YC_t::NYield, 0.);
    TPZManVector<int, YC_t::NYield> validEqs(YC_t::NYield, 0);
    
    REAL normEpsPErr = 0.;
    REAL lambda = 0.;
    
    succeeded = PlasticLoop(stateAtYield, Np1, delGamma, normEpsPErr, lambda, validEqs);
    
    REAL TolEpsPErr = 0.00000001;
    
    if(normEpsPErr < TolEpsPErr && succeeded)
    {
        //        cout << "PLASTIC LOOP CONVERGED " << succeeded << endl;    
        PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());
        return; 
    }
    cout << "PLASTIC LOOP NOT ACHIEVED THE REQUIRED TOL " << succeeded << endl;
    PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());
    
}



template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyStrainComputeDep_Internal(const TPZTensor<REAL> &epsTotal,
															 TPZTensor<REAL> &sigma,
															 TPZFMatrix &Dep)
{

	#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> ApplyStrainComputeDep_Internal *** Imposed epsTotal << " << epsTotal;
		LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
    }
    #endif
	#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout << ">>> ApplyStrainComputeDep_Internal *** Imposed epsTotal << " << epsTotal;
        LOGPZ_INFO(logger,sout.str().c_str());
    }
    #endif

	ApplyStrain_Internal(epsTotal);
	
	#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout << "*** ApplyStrainComputeDep *** \n Calling ComputeDep";
        LOGPZ_INFO(logger,sout.str().c_str());
    }
    #endif
	
	ComputeDep(sigma, Dep);
	
	#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout1, sout2;
        sout1 << "<<< ApplyStrainComputeDep *** Exiting Method.";
        sout2 << "\nImposed epsTotal << " << epsTotal
		      << "\nResulted in Sigma = " << sigma << "\n and Dep = \n" << Dep;
        LOGPZ_INFO(logger,sout1.str().c_str());
        LOGPZ_INFO(logger,sout2.str().c_str());
    }
    #endif

}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ComputeDep(TPZTensor<REAL> & sigma, TPZFMatrix &Dep)
{
	
	const int nyield = YC_t::NYield;
	const int nVarsResidual = 7+nyield;
	const int nVarsTensor = 6;
	
	int i,j;
	
	typedef TFad<nVarsTensor, REAL> TFAD;
	typedef TFad<nVarsResidual, TFAD> TFAD_FAD;
	
	REAL normEpsPErr, resnorm;
	
	EStatus status;
	
	TPZPlasticState<TFAD_FAD> Nk_FADFAD,
	                          Nkp1_FADFAD;
	
	TPZPlasticState<REAL > Nk, Nkp1;
	
	TPZManVector<TFAD_FAD, nyield>  delGamma_FADFAD(nyield);
	TPZManVector<REAL, nyield> 		delGamma(nyield);
	TPZManVector<TFAD_FAD, nVarsResidual> epsRes_FADFAD(nVarsResidual);
	TPZManVector<REAL, nVarsResidual> epsRes(nVarsResidual);
	TPZDiffMatrix< TFAD > tangent_FAD(nVarsResidual,nVarsResidual),
		                  residual_FAD(nVarsResidual,1),
		                  Sol_FAD(nVarsResidual,1);	
	
	TPZPlasticState<TFAD> Nk_FAD;
	TPZTensor<TFAD> sigma_FAD;
	TFAD A_FAD;
	
	TPZTensor<REAL> diffPlasticStrain;
	
	int n = fPlasticMem.NElements();

	if(n < 2)
	{
		#ifdef LOG4CXX
		  {
		    std::stringstream sout;
		    sout << ">>> ComputeDep *** Insufficient Plastic Mem Entries: " << n << ".";
			LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
		  }
		#endif
		#ifdef LOG4CXX_PLASTICITY
		  {
		    std::stringstream sout;
		    sout << ">>> ComputeDep *** Insufficient Plastic Mem Entries: " << n << ".";
		    LOGPZ_ERROR(logger,sout.str().c_str());
		  }
		#endif
		return;
	}
#ifdef LOG4CXX_PLASTICITY
	{
	    std::stringstream sout;
	    sout << ">>> ComputeDep *** Plastic Mem Entries: " << n << ".";
		if( n == 2 ) sout << "\nTwo entries indicate pure elastic step";
	    LOGPZ_INFO(logger,sout.str().c_str());
	}
#endif	
	
	if(n==2)
	{// pure elastic step - no plastic loop necessary
		fPlasticMem[1].fPlasticState.CopyTo(Nk_FAD);
		for(i = 0; i < nVarsTensor; i++)Nk_FAD.fEpsT.fData[i].diff(i,nVarsTensor);
	
	}else
	{// plastification occurs
	
		// plastic Loop analogue with derivative evaluation
	
		// Initializing last available elastic step
		fPlasticMem[1].fPlasticState.CopyTo(Nk_FADFAD);
		for(j = 0; j < nVarsTensor; j++)Nk_FADFAD.fEpsT.fData[j].val().fastAccessDx(j) = fPlasticMem[1].fK;
		
		REAL disturbFactor = sqrt(fResTol); // imposing an initial guess very close to the real
		// solution in order to ensure residual drop and convergence.
		// This is very important to accumulate good FAD derivatives
	
		for(i = 2; i < n; i++)
		{		                  
			InitializePlasticFAD(fPlasticMem[i].fPlasticState,
								 fPlasticMem[i].fDelGamma,
								 Nkp1_FADFAD,
								 delGamma_FADFAD,
								 nVarsResidual);
			fYC.SetForceYield(fPlasticMem[i].fForceYield); // imposing the same assumptions made in the plasticLoop
			
			Nkp1_FADFAD.fAlpha.val().val() -=
				( fPlasticMem[i].fPlasticState.fAlpha - fPlasticMem[i-1].fPlasticState.fAlpha )
				* disturbFactor;
		
			for(j = 0; j < nyield; j++)delGamma_FADFAD[j].val().val() *= (1.0 - disturbFactor);
		
			for(j = 0; j < nVarsTensor; j++)
			{
				Nkp1_FADFAD.fEpsP.fData[j].val().val() -=
					( fPlasticMem[i].fPlasticState.fEpsP.fData[j] - fPlasticMem[i-1].fPlasticState.fEpsP.fData[j] )
					* disturbFactor;
				Nkp1_FADFAD.fEpsT.fData[j].val().diff(j, nVarsTensor);
				Nkp1_FADFAD.fEpsT.fData[j].val().fastAccessDx(j) = fPlasticMem[i].fK; // setting the derivatives of EpsT with respect to deltaEpsT
			}
					
			#ifdef LOG4CXX_PLASTICITY
        	{
        	    std::stringstream sout;
        	    sout << "*** ComputeDep *** Before Matrix Invertion: Plastic step number " << i-1 << " of " << n-1
					 << "\n*Nk_FADFAD=\n" << setw(10) << Nk_FADFAD
					 << "\n*Nkp1_FADFAD=\n" << setw(10) << Nkp1_FADFAD
			    	 << "\n*delGamma_FADFAD=\n" << setw(10) << delGamma_FADFAD;
        		LOGPZ_DEBUG(logger,sout.str().c_str());
        	}
        	#endif	
		
			int NewtonCounter = 0;

			do{//(resnorm > fResTol && NewtonCounter < fMaxNewton)
				
				PlasticResidual<TFAD_FAD, TFAD_FAD>(Nk_FADFAD, Nkp1_FADFAD, 
													delGamma_FADFAD, epsRes_FADFAD, 
													normEpsPErr);	
				tangent_FAD.Reset(); // resets the LU Decomposition flag
		
				ExtractTangent(epsRes_FADFAD, residual_FAD,
								resnorm, tangent_FAD, // TPZFMatrix for T1=fad<real> type
    	            	        fPlasticMem[i].fValidEqs,
							    1 /*precond*/, 1 /*resetInvalidEqs*/);

				status = tangent_FAD.Decompose_LU();
				if(status == EZeroPivot)
				{
					std::stringstream sout;
		    		sout << "*** ComputeDep *** ### Decompose_LU error! - ZeroPivot ### No inversion will be performed";
					#ifdef LOG4CXX
						LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
					#endif
					#ifdef LOG4CXX_PLASTICITY
						LOGPZ_ERROR(logger,sout.str().c_str());
					#endif
					cout << endl << sout.str().c_str();
				//	DebugStop();
				}
			
				Sol_FAD = residual_FAD;
				status = tangent_FAD.Substitution(&Sol_FAD);
				if(status != EOk)
				{
					std::stringstream sout;
		    		if(status == EIncompDim)sout << "*** ComputeDep *** ### LU Substitution error! - IncompatibleDimensions ### No inversion will be performed";
					if(status == EZeroPivot)sout << "*** ComputeDep *** ### LU Substitution error! - ZeroPivot ### No inversion will be performed";
					#ifdef LOG4CXX
						LOGPZ_ERROR(plasticIntegrLogger,sout.str().c_str());
					#endif
					#ifdef LOG4CXX_PLASTICITY
						LOGPZ_ERROR(logger,sout.str().c_str());
					#endif
					cout << endl << sout.str().c_str();
			//		DebugStop();
				}
		
				/*
				REAL lambda = UpdatePlasticVars(Nk_FADFAD, Nkp1_FADFAD, 
											   delGamma_FADFAD, epsRes_FADFAD, 
											   Sol_FAD, fPlasticMem[i].fValidEqs,
											   0); //Do not update variables internally
			
				*/
				REAL lambda = 1.; // forcing unity because the guess is very close to the solution
			
				
				#ifdef LOG4CXX_PLASTICITY
	        	{
	        	    std::stringstream sout;
	        	    sout << "*** ComputeDep *** After " << NewtonCounter
						 << "-th Matrix Invertion: Plastic step number " << i-1 << " of " << n-1
						 << "\nSol_FAD=\n " << setw(10) << Sol_FAD
						 << "\nepsRes_FADFAD=\n " << setw(10) << epsRes_FADFAD
						 << "\nNkp1_FADFAD=\n" << setw(10) <<  Nkp1_FADFAD
						 << "\ndelGamma_FADFAD=\n" << delGamma_FADFAD
						 << "\nfLambda = " << setw(10) << lambda;
	           		LOGPZ_DEBUG(logger,sout.str().c_str());
	        	}
	        	#endif
								
				for(j=0; j<nVarsTensor; j++) Nkp1_FADFAD.fEpsP.fData[j].val() -= lambda*Sol_FAD(j);
    			Nkp1_FADFAD.fAlpha.val() -= lambda*Sol_FAD(j++);
    			for(j=0; j<YC_t::NYield; j++) delGamma_FADFAD[j].val() -= lambda*Sol_FAD(j+7);
			
				for(j=0;j<nVarsTensor;j++)
				{
					Nk  .fEpsP.fData[j] = Nk_FADFAD  .fEpsP.fData[j].val().val();
					Nk  .fEpsT.fData[j] = Nk_FADFAD  .fEpsT.fData[j].val().val();
					Nkp1.fEpsP.fData[j] = Nkp1_FADFAD.fEpsP.fData[j].val().val();
					Nkp1.fEpsT.fData[j] = Nkp1_FADFAD.fEpsT.fData[j].val().val();
				}
				Nk  .fAlpha = Nk_FADFAD.  fAlpha.val().val();
				Nkp1.fAlpha = Nkp1_FADFAD.fAlpha.val().val();
				for(j=0; j<YC_t::NYield; j++) delGamma[j] = delGamma_FADFAD[j].val().val();
								
				PlasticResidual<REAL, REAL>(Nk, Nkp1, delGamma, epsRes, normEpsPErr);	
				//PlasticResidual<TFAD_FAD, TFAD_FAD>(Nk_FADFAD, Nkp1_FADFAD, delGamma_FADFAD, epsRes_FADFAD, normEpsPErr);
			
				// updating the residual
				for(j=0; j < nyield; j++)
					if(fPlasticMem[i].fValidEqs[j] == 0)epsRes[j + 7] = 0.;
				resnorm = 0.;
				for(j=0; j < nVarsResidual; j++)
		            resnorm += pow(/*epsRes_FADFAD[j].val().val()*/epsRes[j] , 2.);
  				resnorm = sqrt(resnorm);

				NewtonCounter++;
			}while(resnorm > fResTol && NewtonCounter < fMaxNewton || NewtonCounter < 2);

			for(j = 0; j < 6; j++)diffPlasticStrain.fData[j] = Nkp1_FADFAD.fEpsP.fData[j].val().val()
				                                               - fPlasticMem[i].fPlasticState.fEpsP.fData[j];
		
			#ifdef LOG4CXX
	    	{
	    	    std::stringstream sout;
				sout << "*** ComputeDep *** substep " << i-1 << " of " << n-2 
			     	 << " solved with " << NewtonCounter << " iterations and residual = " << resnorm;
				if(resnorm > fResTol)
				{
					sout << "\n#### Truncated Newton ####. Results are unpredictable";
					sout << "\nDifferences in the plastic strain:" << diffPlasticStrain;
					sout << "\nDifferences in alpha:" << Nkp1_FADFAD.fAlpha.val().val() - fPlasticMem[i].fPlasticState.fAlpha;
					sout << "\nAlpha = " << fPlasticMem[i].fPlasticState.fAlpha;
					LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
				}else
				{
					LOGPZ_DEBUG(plasticIntegrLogger,sout.str().c_str());
				}
	    	}
	    	#endif
		
			#ifdef LOG4CXX_PLASTICITY
	    	{
	        	std::stringstream sout;
				sout << "*** ComputeDep *** substep " << i-1 << " of " << n-2 
			    	 << " solved with " << NewtonCounter << " iterations and residual = " << resnorm;
				if(resnorm > fResTol)
				{
					sout << "\n#### Truncated Newton ####. Results are unpredictable";
					sout << "\nDifferences in the plastic strain:" << diffPlasticStrain;
					sout << "\nDifferences in alpha:" << Nkp1_FADFAD.fAlpha.val().val() - fPlasticMem[i].fPlasticState.fAlpha;
					sout << "\nAlpha = " << fPlasticMem[i].fPlasticState.fAlpha;
					LOGPZ_WARN(logger,sout.str().c_str());
				}else
				{
					LOGPZ_INFO(logger,sout.str().c_str());
				}
	    	}
	    	#endif
			
			// substep marching
			for(j = 0; j < nVarsTensor; j++)
			{
				Nk_FADFAD.fEpsT.fData[j].val()  = Nkp1_FADFAD.fEpsT.fData[j].val(); // ignoring first derivatives
				Nk_FADFAD.fEpsP.fData[j].val()  = Nkp1_FADFAD.fEpsP.fData[j].val(); // ignoring first derivatives
			}
			Nk_FADFAD.fAlpha.val() = Nkp1_FADFAD.fAlpha.val();// ignoring first derivatives
		}
	
		// decreasing derivative order
		for(i = 0; i < nVarsTensor; i++)
		{
			Nk_FAD.fEpsT.fData[i]  = Nk_FADFAD.fEpsT.fData[i].val();
			Nk_FAD.fEpsP.fData[i]  = Nk_FADFAD.fEpsP.fData[i].val();
		}
		Nk_FAD.fAlpha = Nk_FADFAD.fAlpha.val();
	
	}// end of plasticLoop analogue
	
	//at this point, all the residual variables should contain their real derivatives to detaEpsT
	// in the case the loading is either elastic or elastoplastic.
	ComputePlasticVars(Nk_FAD, sigma_FAD, A_FAD);
	
	for(i = 0; i < nVarsTensor; i++)for(j = 0; j < nVarsTensor; j++)Dep(i,j) = sigma_FAD.fData[i].dx(j); 
	
	#ifdef LOG4CXX_PLASTICITY
    {
       std::stringstream sout;
       sout << "*** ComputeDep *** \nsigma_FAD= \n" << sigma_FAD
			<< "\nDep =\n" << Dep;
       LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif
	
//	NAN detector
//	int res = fetestexcept(FE_ALL_EXCEPT);
//	if(res)
//	{
//		std::cout << " \n " << __PRETTY_FUNCTION__ <<"\n NAN DETECTED \n";
//		DebugStop();
//	}
	
    sigma_FAD.CopyTo(sigma);
	
}

// este metodo usa fN para indicar o ponto "atual" a partir do qual iremos para epsTotalNp1
// o dado de entrada stateAtYield tambem eh inicializado com fN
// stateAtYield eh a saida de um ponto na superficia de plastificacao correspondente a um epstotal modficado

template <class YC_t, class TF_t, class ER_t>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::FindPointAtYield(
				 const TPZTensor<REAL> &epsTotalNp1,
				 TPZPlasticState<REAL> &stateAtYield)const
{
//
//	
//#ifdef LOG4CXX_PLASTICITY
//	
//	int plasticIntegrOutput;
//	
//	{
//		std::stringstream sout;
//		sout << ">>> FindPointAtYield ***";
//		LOGPZ_INFO(logger,sout.str().c_str());
//	}
//#endif
//	
//	const int nVars = 1;
//	typedef TFad<nVars, REAL> TFAD_One;
//	TPZTensor<REAL> deltaEps;
//	TPZTensor<TFAD_One>	deltaEps_FAD, sigma_FAD;
//	TFAD_One A_FAD, multipl_FAD;
//	TPZManVector<TFAD_One,10> phi_FAD(YC_t::NYield);
//	REAL multiplN, minMultipl = 1.;
//	TPZPlasticState<TFAD_One> stateAtYield_FAD;
//	
//	stateAtYield.CopyTo(stateAtYield_FAD);
//	
//	epsTotalNp1.CopyTo(deltaEps);
//	deltaEps.Add(fN.EpsT(), -1.);
//	deltaEps.CopyTo(deltaEps_FAD);
//	
//	int i, nyield = YC_t::NYield;
//	for(i = 0; i < nyield; i++)
//	{ // searching for the multiplier for each yield surface
//		
//    	multiplN = 0.; // the plastic multiplier is likely to be zero
//		// because of possible previous plastifications,
//		// although a null value may lead to no derivative
//		// or numerical instabilities when deltaEpsT=0
//		int count = 0;
//		do
//		{
//			if(count  > 0)
//			{
//				REAL derX = phi_FAD[i].dx(0);
//				if(fabs(derX)<1.e-8) derX += 1.e-8;
//				multiplN -= phi_FAD[i].val() / derX; // avoiding division by zero
//				if(multiplN > 1.0) 
//				{
//#ifdef LOG4CXX_PLASTICITY
//					{
//						std::stringstream sout;
//						sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
//						sout << "\nPlastification is known to occur within this load step and the guess for the nest Newton's step is clipped to multipl=1";
//						LOGPZ_DEBUG(logger,sout.str().c_str());
//					}
//#endif
//					multiplN = 1.0; // The step is known as plastic in advance
//					// (IsStrainElastic previously called)
//					// The check above avoids the multiplicator to be too large when the
//					// first derivative evaluation is very low. In very nonlinear models
//					// and specially in those where the stiffness matrix is very low at the
//					// null stress state this could happen and, if this statement weren't here,
//					// this newton loop could be extremely slow to converge since the
//					// initial guesses would be too far from the correct answer.
//				}
//				
//				if(multiplN < -1.0)
//				{
//#ifdef LOG4CXX_PLASTICITY
//					{
//						std::stringstream sout;
//						sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
//						sout << "\nPlastification is known to occur within this load step. Forcing restart with initial guess biased towards the highest range attempting to help the code reach physically meaningful results";
//						LOGPZ_DEBUG(logger,sout.str().c_str());
//					}
//#endif
//					multiplN = 1.0;  // This check attempts to ensure that the
//					// solution of this plastification step does not explode backwards to the
//					// desired loading path. This may happen when the yield function is too
//					// nonlinear and may present a second solution in the negative range of
//					// the desired loading path. By setting the initial guess equal to 1.0
//					// (that means the whole loading path is plastic - true only in perfect
//					// plastic materials) the code attempts to bias the solver towards the 
//					// physical meaningful solution.
//				}
//			}
//			count++;
//			
//			multipl_FAD = multiplN; 
//			multipl_FAD.diff(0,nVars);
//			fN.EpsT().CopyTo(stateAtYield_FAD.fEpsT);
//			stateAtYield_FAD.fEpsT.Add(deltaEps_FAD, multipl_FAD);
//			
//		    ComputePlasticVars(stateAtYield_FAD,
//							   sigma_FAD,
//							   A_FAD);
//			
//			// compute the value of the yield functions
//			fYC.Compute(sigma_FAD, A_FAD, phi_FAD, 0);
//			
//			
//		}while(fabs (phi_FAD[i].val()) > fResTol && count < fMaxNewton);
//		
//#ifdef LOG4CXX
//        {
//			if(count >= fMaxNewton )
//			{
//            	std::stringstream sout;
//            	sout << "*** FindPointAtYield *** multiplication factor= ";
//				sout << multiplN << " such that yield of function ";
//				sout << i << " = " << phi_FAD[i].val() << "; ";
//				sout << "\n#### Truncated Newton after " 
//				<< fMaxNewton << " steps with phi[" << i << "] = " 
//				<< shapeFAD::val(phi_FAD[i]) << "####. Results are unpredictable";
//				LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//				//plasticIntegrOutput = 1;
//			}
//        }
//#endif
//#ifdef LOG4CXX_PLASTICITY
//        {
//            std::stringstream sout;
//            sout << "*** FindPointAtYield *** multiplication factor= ";
//			sout << multiplN << " such that yield of function ";
//			sout << i << " = " << phi_FAD[i].val();
//			if(count >= fMaxNewton )
//			{
//				sout << "\n#### Truncated Newton after " 
//				<< fMaxNewton << " steps with phi[" << i << "] = " 
//				<< shapeFAD::val(phi_FAD[i]) 
//				<< "####. It appears in such load path the plastic yield solution was found in the opposite direction."
//				<< "\nPlease check the other(s) yield function(s) to guarantee at least one of them is plastifying within the proposed load path.";
//				LOGPZ_WARN(logger,sout.str().c_str());
//				LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//				
//			}else{
//				LOGPZ_DEBUG(logger,sout.str().c_str());
//			}
//        }
//#endif
//		
//		if(multiplN < minMultipl)minMultipl = multiplN;
//	}
//	
//	if(minMultipl < - fResTol)
//	{
//#ifdef LOG4CXX
//        {
//			std::stringstream sout;
//			sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
//			<< " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
//			LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
//        }
//#endif	
//#ifdef LOG4CXX_PLASTICITY
//        {
//			std::stringstream sout;
//			sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
//			<< " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
//			LOGPZ_WARN(logger,sout.str().c_str());
//        }
//#endif
//	}
//	
//	
//	if(minMultipl < 0.)minMultipl = 0.; //avoiding rounding errors
//	
//	stateAtYield = fN;
//	stateAtYield.fEpsT.Add(deltaEps, minMultipl);
//	
//#ifdef LOG4CXX_PLASTICITY
//	{
//		std::stringstream sout;
//		sout << "<<< FindPointAtYield *** Exiting Method with deltaStrainMultiplier = " << minMultipl;
//		LOGPZ_INFO(logger,sout.str().c_str());
//		if(plasticIntegrOutput)LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
//	}
//#endif
//	
//	return minMultipl;
	
	
	int plasticIntegrOutput;
	

	
	const int nVars = 1;
	typedef TFad<nVars, REAL> TFAD_One;
	TPZTensor<REAL> deltaEps;
	TPZTensor<TFAD_One>	deltaEps_FAD, sigma_FAD;
	TFAD_One A_FAD, multipl_FAD;
	TPZManVector<TFAD_One,10> phi_FAD(YC_t::NYield);
	REAL multiplN, minMultipl = 1.;
	TPZPlasticState<TFAD_One> stateAtYield_FAD;
	
	
	stateAtYield.CopyTo(stateAtYield_FAD);
	
	epsTotalNp1.CopyTo(deltaEps);
	deltaEps.Add(fN.EpsT(), -1.);
	deltaEps.CopyTo(deltaEps_FAD);
	
	int i, nyield = YC_t::NYield;
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << ">>> FindPointAtYield ***";
		sout << " \n epsTotalNp1 "<<epsTotalNp1;
		sout << " \n deltaEps "<<deltaEps;
		sout << " \n deltaEps_FAD "<<deltaEps_FAD;
		sout << " \n nyield "<<nyield;
		LOGPZ_INFO(logger,sout.str().c_str());
	}
#endif
	
	
	for(i = 0; i < nyield; i++)
	{ // searching for the multiplier for each yield surface

    	multiplN = 1.; // the plastic multiplier is likely to be zero
		               // because of possible previous plastifications,
		               // although a null value may lead to no derivative
					   // or numerical instabilities when deltaEpsT=0
		int count = 0;
		do
		{
			if(count  > 0)
			{
				REAL derX = phi_FAD[i].dx(0);
				if(fabs(derX)<1.e-8) derX += 1.e-8;
				multiplN -= phi_FAD[i].val() / derX; // avoiding division by zero
				if(multiplN > 1.0) 
				{
					#ifdef LOG4CXX_PLASTICITY
					{
						std::stringstream sout;
						sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
						sout << "\nPlastification is known to occur within this load step and the guess for the nest Newton's step is clipped to multipl=1";
						LOGPZ_DEBUG(logger,sout.str().c_str());
					}
					#endif
					multiplN = 1.0; // The step is known as plastic in advance
				     // (IsStrainElastic previously called)
				     // The check above avoids the multiplicator to be too large when the
				     // first derivative evaluation is very low. In very nonlinear models
				     // and specially in those where the stiffness matrix is very low at the
				     // null stress state this could happen and, if this statement weren't here,
				     // this newton loop could be extremely slow to converge since the
				     // initial guesses would be too far from the correct answer.
				}
				
				if(multiplN < -1.0)
				{
					#ifdef LOG4CXX_PLASTICITY
					{
						std::stringstream sout;
						sout << "*** FindPointAtYield *** multiplication factor = " << multiplN << " set to 1.0";
						sout << "\nPlastification is known to occur within this load step. Forcing restart with initial guess biased towards the highest range attempting to help the code reach physically meaningful results";
						LOGPZ_DEBUG(logger,sout.str().c_str());
					}
					#endif
					multiplN = 1.0;  // This check attempts to ensure that the
				     // solution of this plastification step does not explode backwards to the
				     // desired loading path. This may happen when the yield function is too
				     // nonlinear and may present a second solution in the negative range of
				     // the desired loading path. By setting the initial guess equal to 1.0
				     // (that means the whole loading path is plastic - true only in perfect
				     // plastic materials) the code attempts to bias the solver towards the 
				     // physical meaningful solution.
				}
			}
			count++;
			
			multipl_FAD = multiplN; 
			multipl_FAD.diff(0,nVars);
			fN.EpsT().CopyTo(stateAtYield_FAD.fEpsT);
			stateAtYield_FAD.fEpsT.Add(deltaEps_FAD, multipl_FAD);

		    ComputePlasticVars(stateAtYield_FAD,
							   sigma_FAD,
							   A_FAD);

			// compute the value of the yield functions
			fYC.Compute(sigma_FAD, A_FAD, phi_FAD, 0);
			
			// if at the first count phi<0 there is no need to compute a multiplier...
			if(count == 1 && phi_FAD[i].val() < 0.) 
			{
				break;
			}
			
			
		}while(fabs (phi_FAD[i].val()) > fResTol && count < fMaxNewton);

#ifdef LOG4CXX
        {
			if(count >= fMaxNewton )
			{
            	std::stringstream sout;
            	sout << "*** FindPointAtYield *** multiplication factor= ";
				sout << multiplN << " such that yield of function ";
				sout << i << " = " << phi_FAD[i].val() << "; ";
				sout << "\n#### Truncated Newton after " 
					 << fMaxNewton << " steps with phi[" << i << "] = " 
					 << shapeFAD::val(phi_FAD[i]) << "####. Results are unpredictable";
				LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
				//plasticIntegrOutput = 1;
			}
        }
#endif
#ifdef LOG4CXX_PLASTICITY
        {
            std::stringstream sout;
            sout << "*** FindPointAtYield *** multiplication factor= ";
			sout << multiplN << " such that yield of function ";
			sout << i << " = " << phi_FAD[i].val();
			if(count >= fMaxNewton )
			{
				sout << "\n#### Truncated Newton after " 
					 << fMaxNewton << " steps with phi[" << i << "] = " 
					 << shapeFAD::val(phi_FAD[i]) 
					 << "####. It appears in such load path the plastic yield solution was found in the opposite direction."
					 << "\nPlease check the other(s) yield function(s) to guarantee at least one of them is plastifying within the proposed load path.";
				LOGPZ_WARN(logger,sout.str().c_str());
				LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
				
			}else{
				LOGPZ_DEBUG(logger,sout.str().c_str());
			}
        }
#endif
		
		if(multiplN < minMultipl)minMultipl = multiplN;
	}
	
	if(minMultipl < - fResTol)
	{
        #ifdef LOG4CXX
        {
           std::stringstream sout;
           sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
			    << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
		   LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
        }
        #endif	
        #ifdef LOG4CXX_PLASTICITY
        {
           std::stringstream sout;
           sout << "*** FindPointAtYield *** Ignoring deltaStrainMultiplier = " << minMultipl
			    << " and setting it to zero.\n\t###### WARN: EpsTotalN isn't elastic!! ######";
           LOGPZ_WARN(logger,sout.str().c_str());
        }
        #endif
	}
	
	
	if(minMultipl < 0.)minMultipl = 0.; //avoiding rounding errors
	
	stateAtYield = fN;
	stateAtYield.fEpsT.Add(deltaEps, minMultipl);

#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << "<<< FindPointAtYield *** Exiting Method with deltaStrainMultiplier = " << minMultipl;
    LOGPZ_INFO(logger,sout.str().c_str());
	if(plasticIntegrOutput)LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
  }
#endif
	
	return minMultipl;
}


template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticLoop(const TPZPlasticState<REAL> &N,TPZPlasticState<REAL> &Np1,TPZVec<REAL> &delGamma,REAL &normEpsPErr,REAL &lambda,TPZVec<int> & validEqs)
{

  // 6 variables for tensor and others: alpha and deltaGamma
  const int nVars = 7+YC_t::NYield;
  int i;
  const REAL Tol = fResTol;
	
#ifdef LOG4CXX_PLASTICITY
  {
	  std::stringstream sout;//1, sout2;
    sout << ">>> PlasticLoop ***";
    sout << "\nNp1 << \n" << Np1
	     << "\n N << \n "<< N
          << "\ndelGamma << " << delGamma
	  << "\n labmbda << " << lambda
          << "\nNumber of plasticity variables: " << nVars;
    LOGPZ_INFO(logger,sout.str().c_str());
    LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif

  typedef TFad<nVars, REAL> TFAD;

  TPZPlasticState<TFAD>    Np1_FAD;
  TPZTensor<TFAD>          sigmaNp1_FAD;
  TPZManVector<TFAD,nVars> epsRes_FAD(nVars), 
	                       delGamma_FAD(YC_t::NYield);
  TFAD                     phiRes_FAD;
  TPZFNMatrix<nVars>       ResVal(nVars,1,0.), // ResVal to hold the residual vector
					       Sol(nVars,1,0.);    // Sol will contain the values of the unknown values
  TPZFNMatrix<nVars*nVars> tangent(nVars,nVars,0.); // Jacobian matrix
  REAL                     resnorm = 0.;

  int countReset = 0;	
  int countNewton = 0;
	
	ofstream arg1("PlasticLoopNewton.txt");
	
  do{//while(RemoveInvalidEqs(epsRes_FAD));
	  
	 // verifying if it is necessary to impose post-peak material behavior
	 TPZTensor<REAL> sigmaGuess;
	 REAL AGuess;
	 ComputePlasticVars<REAL>(Np1  , sigmaGuess  , AGuess);//Chute inicial
	 fYC.SetYieldStatusMode(sigmaGuess, AGuess);
	  
     InitializePlasticFAD(Np1, delGamma, Np1_FAD, delGamma_FAD);

     PlasticResidual<REAL, TFAD>(N, Np1_FAD, delGamma_FAD, epsRes_FAD, normEpsPErr);
	  
	 if(countReset == 0)InitializeValidEqs(epsRes_FAD, validEqs);
	
#ifdef LOG4CXX_PLASTICITY
        {
   	     std::stringstream sout;
  	     sout << "*** PlasticLoop *** PlasticLoop main loop with " 
              << countReset << "-th valid equation set: {" 
              << validEqs << "}.";
   	     LOGPZ_INFO(logger,sout.str().c_str());
		}
#endif	  
	  
     ExtractTangent(epsRes_FAD, ResVal, resnorm, tangent, validEqs, 1/*precond*/, 1/*ResetUnvalidEqs*/); 
	  
	 countNewton = 0;
	
	  
	 do{

		 countNewton++;
	
         TPZFNMatrix<nVars*nVars> *matc = new TPZFNMatrix<nVars*nVars>(nVars,nVars);
         *matc = tangent;
         TPZStepSolver st(matc);
	  // enum DecomposeType {ENoDecompose, ELU, ELUPivot, ECholesky, ELDLt};
         st.SetDirect(ELU);
	  // st.SetDirect(ECholesky);
	  // st.SetDirect(ENoDecompose); 
		 
//		 cout << " RESVAL " <<endl;
//		 cout << ResVal << endl;
	
		 
		 //QUEBRA AQUI!
		 // invert the tangent matrix and put the correction in the Sol variable
         st.Solve(ResVal,Sol,0);

		 // update the independent variables (their partial derivatives remain unchanged)
         lambda = UpdatePlasticVars(N, Np1_FAD, delGamma_FAD, epsRes_FAD, Sol, validEqs);

         // recompute the residual
         PlasticResidual<REAL, TFAD>(N, Np1_FAD, delGamma_FAD, epsRes_FAD, normEpsPErr);
#ifdef LOG4CXX_PLASTICITY
         {
             std::stringstream sout1;
			 sout1 << "*** PlasticLoop *** Newton's " << countNewton; 
			 sout1 << "\n PLASTICRESIDUAL \n";
			 sout1 << "\n N \n" << N << "\n";
			 sout1 << "\n Np1_FAD \n" << Np1_FAD << "\n";
			 sout1 << "\n delGamma_FAD \n" << delGamma_FAD << "\n";
			 sout1 << "\n epsRes_FAD \n" << epsRes_FAD << "\n";
			 sout1 << "\n normEpsPErr \n" << normEpsPErr << "\n";
			 LOGPZ_INFO(loggerPlasticLoop,sout1.str().c_str());
         }
#endif
		 // extract the values of the residual vector
         ExtractTangent(epsRes_FAD, ResVal, resnorm, tangent, validEqs, 1, 1);
		 
#ifdef LOG4CXX_PLASTICITY
         {
             std::stringstream sout1;
			 sout1 << "*** PlasticLoop *** Newton's " << countNewton; 
			 sout1 << "\n EXTRACT TANGENT \n";
			 sout1 << "\n epsRes_FAD \n" << epsRes_FAD << "\n";
			 sout1 << "\n ResVal \n" << ResVal << "\n";
			 sout1 << "\n resnorm \n" << resnorm << "\n";
			 sout1 << "\n tangent \n" << tangent << "\n";
			 sout1 << "\n validEqs \n" << validEqs << "\n";
          LOGPZ_INFO(loggerPlasticLoop,sout1.str().c_str());
          //LOGPZ_DEBUG(logger,sout2.str().c_str());
         }
#endif
  		  		  
	  }while(resnorm > Tol && countNewton < fMaxNewton);
	
#ifdef LOG4CXX
      {
		  
		  std::stringstream sout;
          sout << "*** PlasticLoop *** Exiting Newton's scheme after " << countNewton
		       << " steps and with residual norm = " << resnorm;
		  if(resnorm > Tol)
		  {
			  sout << "\n###### Max Newton steps achieved - Truncating Newton ######";
			  LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
		  }else
		  {
			  LOGPZ_DEBUG(plasticIntegrLogger,sout.str().c_str());
		  }
      }
#endif
#ifdef LOG4CXX_PLASTICITY
      {
      std::stringstream sout1, sout2;
      sout1 << "*** PlasticLoop *** Exiting Newton's scheme after " << countNewton
		    << " steps and with residual norm = " << resnorm;
	  LOGPZ_INFO(logger,sout1.str().c_str());
	  if(resnorm > Tol)
		  {
			  sout2 << "\n###### Max Newton steps achieved - Truncating Newton ######";
			  LOGPZ_WARN(logger,sout2.str().c_str());
		  }
      }
#endif

      countReset++;
			
  }while(RemoveInvalidEqs(delGamma_FAD, epsRes_FAD, validEqs));
	 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  // updating output variables
  //for(i=0; i<6; i++) Np1.fEpsP.fData[i] = Np1_FAD.EpsP().fData[i].val();
  //Np1.fAlpha = Np1_FAD.Alpha().val();
  Np1_FAD.CopyTo(Np1);
	
  for(i=0; i<YC_t::NYield; i++)
    delGamma[i] = delGamma_FAD[i].val();

#ifdef LOG4CXX_PLASTICITY
    {
     std::stringstream sout1, sout2;
     sout1 << "<<< PlasticLoop *** Exiting Method";
     sout2 << "\nNp1.Alpha() << " << Np1.Alpha()
           << "\ndelGamma << " << delGamma
		   << "\nValidEqs << {" << validEqs << "}" ;
     LOGPZ_INFO(logger,sout1.str().c_str());
     LOGPZ_INFO(logger,sout2.str().c_str());
    }
#endif
	
	return resnorm<=Tol;		
}

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticIntegrate(const TPZPlasticState<REAL> &N,TPZPlasticState<REAL> &Np1,const REAL &TolEpsPErr)
{
	int succeeded;
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> PlasticIntegrate *** Np1.fEpsT = "
	     << Np1.fEpsT;
	LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
  }
#endif
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << ">>> PlasticIntegrate *** Np1.fEpsT = "
	     << Np1.fEpsT;
	  sout << ">>> PlasticIntegrate *** N.fEpsT = "
	  << N.fEpsT << "\n";

	  sout << "PRINT Np1 \n";
	  Np1.Print(sout);
	  sout << "PRINT N \n";
	  N.Print(sout);
    LOGPZ_INFO(logger,sout.str().c_str());
  }
#endif

    TPZManVector<REAL,YC_t::NYield> delGamma(YC_t::NYield, 0.);
	TPZManVector<int, YC_t::NYield> validEqs(YC_t::NYield, 0);

	REAL normEpsPErr = 0.;
	int counter = 0;
	REAL lambda = 0.;

	succeeded = PlasticLoop(N, Np1, delGamma, normEpsPErr, lambda, validEqs);
	
	if(normEpsPErr < TolEpsPErr && succeeded)
	{
		#ifdef LOG4CXX
 		 {
  		  std::stringstream sout;
  		  sout << "*** PlasticIntegrate *** The integration succeeded with the desired tolerance without any substepping";
		  sout << "\nValidEqs = " << validEqs << " and forceYield = " << fYC.GetForceYield();
		  LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
         }
        #endif
		#ifdef LOG4CXX_PLASTICITY
 		 {
  		  std::stringstream sout;
  		  sout << "*** PlasticIntegrate *** The integration succeeded with the desired tolerance without any substepping";
		  sout << "\nValidEqs = " << validEqs << " and forceYield = " << fYC.GetForceYield();
		  LOGPZ_INFO(logger,sout.str().c_str());
         }
        #endif
		
		PushPlasticMem(Np1, 1., lambda, delGamma, validEqs, fYC.GetForceYield());
				
		return 1; 
	}
	
	int discarded = 1;
	
	// Substepping state variables
	TPZPlasticState<REAL> Nk(N),
		                  Nkp1;
	
	TPZTensor<REAL> deltaEpsTotal(Np1.EpsT());
	deltaEpsTotal.Add(N.EpsT(), -1.);
	
	REAL k = 0., q = 1., kp1 = 0., multipl;

	while(k < 1.)
	{
		
		multipl = 0.95 * pow(TolEpsPErr/normEpsPErr, 0.5);
		
		if(multipl < 0.1) multipl = 0.1;
		if(multipl > 10.) multipl = 10.;
		
		if(!succeeded)multipl = 0.5;
		
		q*= multipl;
		
	    if(q < fMinStepSize)q = fMinStepSize;
		
		kp1 = min(1.0, k + q);
		q = kp1 - k; // needed when k+q exceeds 1
		Nkp1.fEpsT = N.EpsT();
		Nkp1.fEpsT.Add(deltaEpsTotal, kp1);
		Nkp1.fAlpha = Nk.Alpha();
		Nkp1.fEpsP  = Nk.EpsP();
		for(int i = 0; i < YC_t::NYield; i++)delGamma[i] = 0.;
		
        succeeded = PlasticLoop(Nk, Nkp1, delGamma, normEpsPErr, lambda, validEqs);
		
		#ifdef LOG4CXX
 		 {
  		  std::stringstream sout;
  		  sout << "*** PlasticIntegrate *** " << counter 
               << "-th substep with k=" << k << " and kp1=" << kp1
		       << " (normEpsPErr = " << normEpsPErr << " )"
			   << "\nValidEqs = " << validEqs << " and forceYield = " << fYC.GetForceYield();
		  LOGPZ_DEBUG(plasticIntegrLogger,sout.str().c_str());
         }
        #endif
		#ifdef LOG4CXX_PLASTICITY
 		 {
  		  std::stringstream sout;
  		  sout << "*** PlasticIntegrate *** " << counter 
               << "-th substep with k=" << k << " and kp1=" << kp1
		       << " (normEpsPErr = " << normEpsPErr << " )"
			   << "\nValidEqs = " << validEqs << " and forceYield = " << fYC.GetForceYield();
          LOGPZ_INFO(logger,sout.str().c_str());
         }
        #endif
		
    	if((normEpsPErr < TolEpsPErr && succeeded) || kp1-k < fMinStepSize * 1.001) // 1.001 because of rounding errors
		{
			
           #ifdef LOG4CXX
 		   if(normEpsPErr >= TolEpsPErr)
		   {
  		      std::stringstream sout;
  		      sout << "*** PlasticIntegrate *** ###### SUBSTEPPING CUTOFF ###### "
			       << "\nAccepting substep with normEpsPErr = " << normEpsPErr 
                   << " because the step became " <<(kp1-k)*100. << "% of the original step size";
              LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
           }
           #endif
           #ifdef LOG4CXX_PLASTICITY
 		   if(normEpsPErr >= TolEpsPErr)
		   {
  		      std::stringstream sout;
  		      sout << "*** PlasticIntegrate *** ###### SUBSTEPPING CUTOFF ###### "
			       << "\nAccepting substep with normEpsPErr = " << normEpsPErr 
                   << " because the step became " <<(kp1-k)*100. << "% of the original step size";
              LOGPZ_WARN(logger,sout.str().c_str());
           }
           #endif
			
		   PushPlasticMem(Nkp1, kp1, lambda, delGamma, validEqs, fYC.GetForceYield());
           
		   counter++;
		   // the k-th integration step respects the estimated tolerance
		   // proceeding with time evolution...
		   Nk = Nkp1;			
		    k =  kp1;
		}
		else{
		   discarded++;
		}// otherwise the answer isn't accepted and the next q will be
		// recomputed in order to lie within the desired tolerance.
		// If the method works fine this situation should rarely happen.
	}
	
	#ifdef LOG4CXX
    {
       std::stringstream sout;
       sout << "*** PlasticIntegrate *** ###### Exiting with " << counter << " substepping(s) and "
			<< discarded << " plasticLoop(s) discarded. Residual = " << normEpsPErr << " ######";
	   LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
     }
#endif
	#ifdef LOG4CXX_PLASTICITY
    {
       std::stringstream sout;
       sout << "*** PlasticIntegrate *** ###### Exiting Method with " << counter << " substepping(s) and " << discarded << " plasticLoop(s) discarded ######";
       LOGPZ_INFO(logger,sout.str().c_str());
     }
#endif
	
	return counter;
}

template <class YC_t, class TF_t, class ER_t>
template <class T>
int TPZPlasticStep<YC_t, TF_t, ER_t>::InitializeValidEqs(TPZVec<T> &res_T, TPZVec<int> & validEqs)
{
#ifdef LOG4CXX_PLASTICITY
    {
		int i, n = res_T.NElements();
    	std::stringstream sout;
		sout << ">>> InitializeValidEqs *** Res = ";
		for(i = 0; i < n; i++)sout << shapeFAD::val(res_T[i]) << " ";
		LOGPZ_DEBUG(logger,sout.str().c_str());
    }
#endif
	const int nyield = YC_t::NYield;
	validEqs.Resize(nyield);
	int i, count = 0;
	for(i = 0; i < nyield; i++)
	{
		validEqs[i] = 0;
		if(shapeFAD::val(res_T[i+7]) > 0.)
		{
			validEqs[i] = 1;
			count++;
		}
	}
	
	if(nyield == 1) // if there exists only one eqn it must be active
	{
		validEqs[0] = 1;
		count++;
	}
	
	return count;
}

template <class YC_t, class TF_t, class ER_t>
template <class T>
int TPZPlasticStep<YC_t, TF_t, ER_t>::RemoveInvalidEqs( TPZVec<T> & delGamma_T, TPZVec<T> & res_T, TPZVec<int> & validEqs)
{	
	const int nyield = YC_t::NYield;
	validEqs.Resize(nyield);
	
	if(nyield == 1)return 0;// remove Invalid Eqs should not invalidate the unique equation
	
	int i, count = 0;
	
#ifdef LOG4CXX_PLASTICITY
    {
		int i;
		std::stringstream sout;
		sout << ">>> RemoveInvalidEqs *** delGamma = ";
		for(i = 0; i < nyield; i++)sout << shapeFAD::val(delGamma_T[i]) << " ";
		sout << " phi = ";
		for(i = 0; i < nyield; i++)sout << shapeFAD::val(res_T[i+7]) << " ";
		LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif
	
	int BoolDelGamma,
	    BoolValidEq,
	    BoolRes;
	
	for(i = 0; i < nyield; i++)
	{
		BoolDelGamma = shapeFAD::val(delGamma_T[i]) > 0.; // deltaGamma is valid
		BoolValidEq  = validEqs[i];                      // equation is already in use
		BoolRes      = shapeFAD::val(res_T[i+7]) > fResTol;   // equation indicates plastification
			
		if(BoolRes) validEqs[i] = 1;
		switch(BoolRes)
		{
	    	case (true):
				// the lines below unfortunately led to instabilities in the integration process
				//validEqs[i] = 1; // if the equation indicates plastification then it is necessary to involve it in the process
				//count++;
            break;			
		    case (false):
                 // if the equation does not indicate plastification but is
				 // already involved in the plastification process and 
				 // has a valid deltaGamma then keep it in the integration process - do nothing
				 
				 // if it is involved in the plastification process but has an
				 // invalid deltagamma - invalidate it
				 if(BoolValidEq && !BoolDelGamma)
                     {
						 validEqs[i] = 0;
						 count++;
					 }
			
			     // if it isn't invoved in the plastification process - delgamma is meaningless and no
				 // testing should be performed.
		    break;
		}
	}
	
	return count;
}


template <class YC_t, class TF_t, class ER_t>
template<class T1, class T_VECTOR, class T_MATRIX>//T1:input residual fad type (FAD FAD or FAD), T_MATRIX: output matrix &vector of type (FAD or REAL, respectvly)
void TPZPlasticStep<YC_t, TF_t, ER_t>::ExtractTangent(
						const TPZVec<T1> & epsRes_FAD,
						T_VECTOR & ResVal, //TPZFMatrix for the T1=fad<real> type
						REAL & resnorm, // REAL 
						T_MATRIX & tangent, // TPZFMatrix for T1=fad<real> type
                        TPZVec<int> & validEqs,
						const int precond,
						const int resetInvalidEqs)
{
  const int nyield = YC_t::NYield;
  const int nVars = 7+YC_t::NYield;
  int i,j;

  // extract the partial derivatives to form the tangent matrix
  for(i=0; i<nVars; i++)
  {
     for(j=0; j<nVars; j++) 
        tangent(i,j) = epsRes_FAD[i].dx(j);
     ResVal(i,0) = epsRes_FAD[i].val();
  }
	
//	cout << "epsRes_FAD[i] "<<endl;
//	cout << epsRes_FAD <<endl;
//	
//	cout << "ResVal"<<endl;
//	cout << ResVal<<endl;
//	
//	cout << " TANGENTE " <<endl;
//	cout << tangent <<endl;
	
  // reseting the equations related to the invalid yield surfaces
  if(resetInvalidEqs)
	for(i=0; i<nyield; i++)	
		{
			if(validEqs[i] == 0)
			{
		        for(j=0; j<nVars; j++) 
                   tangent(i+7,j) = 0.;
				tangent(i+7,i+7) = 1.;
				ResVal(i+7,0) = 0.;
			}
		}
  // updating residual norm
  // before preconditioning such that it keeps its physical meaning
	//cout << "ResVal"<<endl;
//	cout << ResVal<<endl;
  resnorm = 0.;
  for(i = 0; i < nVars; i++)
		resnorm += pow(shapeFAD::val( ResVal(i,0) ) , 2.);
  resnorm = sqrt(resnorm);
	//cout << "resnorm"<<endl;
//	cout << resnorm<<endl;
  // preconditioning the matrix/residual system
  if(precond)
     for(i=/*7*/0; i<nVars; i++)
     {
        REAL pivot = 0., elem;
        for(j=0; j<nVars; j++)
		 {
			 elem = shapeFAD::val(tangent(i,j) );
			 if(fabs(pivot) < fabs(elem))pivot = elem;
             //pivot = fabs(shapeFAD::val(pivot) ) < fabs(shapeFAD::val(tangent(i,j) ) ) ? tangent(i,j) : pivot;
		 }
        for(j=0; j<nVars; j++)
           tangent(i,j) = tangent(i,j) / pivot;
        ResVal(i,0) = ResVal(i,0) / pivot;
	 }	
	
#ifdef LOG4CXX_PLASTICITY
  {
    std::stringstream sout;
    sout << ">>> ExtractTangent *** Residual norm = " << resnorm 
		 << "; precond = " << precond
		 << "; resetInvalidEqs = " << resetInvalidEqs 
		 << "; validEqs = {" << validEqs << "}";
    LOGPZ_INFO(logger,sout.str().c_str());
  }
#endif
}

template <class YC_t, class TF_t, class ER_t>
template <class T1, class T2>
void TPZPlasticStep<YC_t, TF_t, ER_t>::PlasticResidual(
				const TPZPlasticState<T1> &N_T1,
				TPZPlasticState<T2> &Np1_T2,
				const TPZVec<T2> &delGamma_T2,
				TPZVec<T2> &res_T2,
				REAL &normEpsPErr,
				int silent)const
{

  const REAL a = 0.5;
	
  // This function will be either called with template parameter T being REAL or FAD type
  // nyield indicates the number of yield functions
  const int nyield = YC_t::NYield;

  // variable to hold the value of the stress, elastic strain and residual corresponding to the plastic strain
  TPZTensor<T2> sigmaNp1_T2, epsResMidPt_T2;
  TPZTensor<T1> sigmaN_T1;
  // variable containing the error or deviation between solutions from Modified and
  // explicit Euler schemes
  TPZTensor<REAL> epsPErr;
	
  // vector holding the tensors of the N directions (usually derivative of the yield functions)
  TPZManVector<TPZTensor<T2>,nyield> NdirNp1_T2(nyield), NdirMidPt_T2(nyield);
  TPZManVector<TPZTensor<T1>, nyield> NdirN_T1(nyield);

  // vector holding the residual of the yield function equations
  TPZManVector<T2, nyield> phiRes_T2(nyield),HNp1_T2(nyield), HMidPt_T2(nyield);
  TPZManVector<T1, nyield> HN_T1(nyield);

  // residual of the equation corresponding to the damage variable
  T2 alphaRes_T2, ANp1_T2;
  T1 AN_T1;
	
  ComputePlasticVars<T1>(N_T1  , sigmaN_T1  , AN_T1);
  ComputePlasticVars<T2>(Np1_T2, sigmaNp1_T2, ANp1_T2);	

  // Compute the values of the N directions (which will update the plastic strain)
  fYC.N(sigmaN_T1  , AN_T1  , NdirN_T1  , 0);
  fYC.N(sigmaNp1_T2, ANp1_T2, NdirNp1_T2, 1);
		
  int i, j;
  for(i=0; i< nyield; i++)for(j = 0; j < 6; j++)NdirMidPt_T2[i].fData[j] = ( ( NdirNp1_T2[i].fData[j] ) * T2(1.-a) 
									 + ( NdirN_T1[i].fData[j]   ) * T1(a) );
	
  // Compute the value of H
  fYC.H(sigmaN_T1  , AN_T1  , HN_T1  , 0);
  fYC.H(sigmaNp1_T2, ANp1_T2, HNp1_T2, 1);
	
  for(i=0; i< nyield; i++)HMidPt_T2[i] = (HNp1_T2[i] * T2(1.-a) + T2( HN_T1[i] * T1(a) ) );
	
  //EpsRes = EpsilonPNp1 - fEpsP - deltagamma * Ndir; // 6 eqs
  //for
  for(i=0; i< nyield; i++)
  {
    NdirMidPt_T2[i].Multiply(delGamma_T2[i], T2(1.) );
    epsResMidPt_T2.Add(NdirMidPt_T2[i], T2(-1.) );
	  
	NdirN_T1[i].Multiply(T1(shapeFAD::val(delGamma_T2[i])), T1(1.) );
	for(j = 0; j < 6; j++)epsPErr.fData[j] +=   shapeFAD::val( NdirN_T1[i].    fData[j] )
		                      				  - shapeFAD::val( NdirMidPt_T2[i].fData[j] );
  }
  epsResMidPt_T2.Add(  N_T1.EpsP(), T1(-1.) );
  epsResMidPt_T2.Add(Np1_T2.EpsP(), T1( 1.) );
  
  // Explicit scheme relative error: estimated by the difference between the
  // first and second order schemes results.
  normEpsPErr = epsPErr.Norm();
  // The second order scheme error is estimated based on the explicit Euler
  // scheme. Its relative error measure coincides with the Explicit Scheme
  // error estimation.

  // alphaRes = alpha(n+1)-alpha(n)-Sum delGamma * H
  alphaRes_T2 = Np1_T2.Alpha() - N_T1.Alpha();
  for(i=0; i<nyield; i++)
  {
    alphaRes_T2 -= delGamma_T2[i] * HMidPt_T2[i]; // 1 eq
  }

  // compute the values of the yield functions
  fYC.Compute(sigmaNp1_T2, ANp1_T2,phiRes_T2, 1); // nyield eq

  // transfer the values of the residuals to the residual vector
  // depending on the type T, the residual and its partial derivatives will have been computed
  for(i=0; i<6; i++)
    res_T2[i] = epsResMidPt_T2.fData[i];
  res_T2[i++] = alphaRes_T2;
  for(j=0; j<nyield; j++)
    res_T2[i++] = phiRes_T2[j];

#ifdef LOG4CXX
 if(!silent)
	{
	int i, j;
	std::stringstream sout1, sout2;
		
		
//		 
//	sout1 << ">>> PlasticResidual *** normEpsPErr = " << normEpsPErr;
//	
//	sout2 << "\n\tN_T1 <<\n";
//	N_T1.Print(sout2, 0 /*no Fad Derivatives*/);
//		
//	sout2 << "\n\tNp1_T2 <<\n";
//	Np1_T2.Print(sout2, 0 /*no Fad Derivatives*/);
//		
//	sout2 << "\n\tdelGamma_T2 <<";
//	for(i = 0; i < YC_t::NYield; i++)sout2 << shapeFAD::val(delGamma_T2[i]) << " ";
//		
//	sout2 << "\tepsPErr <<" << epsPErr;
//		
//	sout2 << "\n\tHMidPt_T2 << ";
//	for(i = 0; i < YC_t::NYield; i++)sout2 << shapeFAD::val(HMidPt_T2[i]) << " ";
//		
//	sout2 << "\tHN_T1 << ";		
//	for(i = 0; i < YC_t::NYield; i++)sout2 << shapeFAD::val(HN_T1[i]) << " ";
//	
//	sout2 << "\tHNp1_T2 << ";	
//	for(i = 0; i < YC_t::NYield; i++)sout2 << shapeFAD::val(HNp1_T2[i]) << " ";
//		
//	sout2 << "\n\tNdirMidPt_T2 <<";		
//	for(i = 0; i < YC_t::NYield; i++)
//		  for(j = 0; j < 6; j++)sout2 << shapeFAD::val(NdirMidPt_T2[i].fData[j]) << " ";
//   
//	sout2 << "\n\tNdirN_T1 <<";
//	for(i = 0; i < YC_t::NYield; i++)
//		  for(j = 0; j < 6; j++)sout2 << shapeFAD::val(NdirN_T1[i].fData[j]) << " ";
//		
//	sout2 << "\n\tNdirNp1_T2 <<";		
//	for(i = 0; i < YC_t::NYield; i++)
//		  for(j = 0; j < 6; j++)sout2 << shapeFAD::val(NdirNp1_T2[i].fData[j]) << " ";
//
//	sout2 << "\n\tres_T2 <<";		
//	for(i = 0; i < YC_t::NYield + 7; i++)
//			sout2 << /*shapeFAD::val*/(res_T2[i]) << " ";
//			
//	sout2 << "\n\tsigmaN_T1 <<\n";
//	sout2 << sigmaN_T1;
//		
//	sout2 << "\n\tsigmaNp1_T2 <<\n";
//	sout2 << sigmaNp1_T2;

		
	
		
////////////////////////////////MYLOG/////////////////////////
		
		std::stringstream sout;
		
		sout << " \n ****----------------MYLOG----------------**** \n";
		
		sout << ">>> PlasticResidual *** normEpsPErr = " << normEpsPErr;
		
		sout << "\n\tN_T1 <<\n";
		N_T1.Print(sout, 0 /*no Fad Derivatives*/);
		
		sout << "\n\tNp1_T2 <<\n";
		Np1_T2.Print(sout, 0 /*no Fad Derivatives*/);
		
		sout << "\n\tdelGamma_T2 <<";
		for(i = 0; i < YC_t::NYield; i++)sout << shapeFAD::val(delGamma_T2[i]) << " ";
		
		sout << "\tepsPErr <<" << epsPErr;
		
		sout << "\n\tHMidPt_T2 << ";
		for(i = 0; i < YC_t::NYield; i++)sout << shapeFAD::val(HMidPt_T2[i]) << " ";
		
		sout << "\tHN_T1 << ";		
		for(i = 0; i < YC_t::NYield; i++)sout << shapeFAD::val(HN_T1[i]) << " ";
		
		sout << "\tHNp1_T2 << ";	
		for(i = 0; i < YC_t::NYield; i++)sout << shapeFAD::val(HNp1_T2[i]) << " ";
		
		sout << "\n\tNdirMidPt_T2 <<";		
		for(i = 0; i < YC_t::NYield; i++)
			for(j = 0; j < 6; j++)sout << shapeFAD::val(NdirMidPt_T2[i].fData[j]) << " ";
		
		sout << "\n\tNdirN_T1 <<";
		for(i = 0; i < YC_t::NYield; i++)
			for(j = 0; j < 6; j++)sout << shapeFAD::val(NdirN_T1[i].fData[j]) << " ";
		
		sout << "\n\tNdirNp1_T2 <<";		
		for(i = 0; i < YC_t::NYield; i++)
			for(j = 0; j < 6; j++)sout << shapeFAD::val(NdirNp1_T2[i].fData[j]) << " ";
		
		sout << "\n\tres_T2 <<";		
		for(i = 0; i < YC_t::NYield + 7; i++)
			sout << /*shapeFAD::val*/(res_T2[i]) << " ";
		
		sout << "\n\tsigmaN_T1 <<\n";
		sout << sigmaN_T1;
		
		sout << "\n\tsigmaNp1_T2 <<\n";
		sout << sigmaNp1_T2;
		
		sout << " \n ****------------ENDMYLOG-----------------**** \n";
		
		
	//LOGPZ_DEBUG(loggerPlasticResidual,sout.str().c_str());	
	//LOGPZ_INFO(logger,sout1.str().c_str());
	//LOGPZ_DEBUG(logger,sout2.str().c_str());
  }
#endif
	
}

template <class YC_t, class TF_t, class ER_t>
template<class T1, class T2, class TVECTOR>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::UpdatePlasticVars(
				const TPZPlasticState<T1> &N_T1,
				TPZPlasticState<T2> &Np1_T2,
				TPZVec<T2> &delGamma_T2,
				TPZVec<T2> &res_T2,
				TVECTOR & Sol_TVECTOR,
				TPZVec<int> & validEqs,
				int updateVars)const
{
	const int nyield = YC_t::NYield;
	const int nVars = 7+nyield;
	TPZManVector<REAL,nyield> delGamma(nyield);
	REAL sqrNormResN = 0, sqrNormResNp1, normEpsPErr;
	TPZManVector<REAL, nVars> res(nVars);
	TPZPlasticState<REAL> Np1, N;
	int i, j;
	
	Np1_T2.CopyTo(Np1);
	N_T1.CopyTo(N);
		
	for(i=0; i<6; i++) Np1.fEpsP.fData[i] -= shapeFAD::val( Sol_TVECTOR(i) );
	Np1.fAlpha -= shapeFAD::val( Sol_TVECTOR(i++) );
	for(j=0; j<nyield; j++)delGamma[j] = shapeFAD::val(delGamma_T2[j]) - shapeFAD::val( Sol_TVECTOR(i++) );
		
	PlasticResidual<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);
	
	sqrNormResNp1 = pow(sdot(res,res),0.5);	
	
	const REAL k = 2.;
	REAL lambda = 1.; // ensuring that the lambda value will be ONE at the first step
	
	do{
		sqrNormResN = sqrNormResNp1;
		
		lambda /= k;
		
		Np1_T2.CopyTo(Np1);
		for(i=0; i<6; i++) Np1.fEpsP.fData[i] -= lambda * shapeFAD::val( Sol_TVECTOR(i) );
		Np1.fAlpha -= lambda * shapeFAD::val( Sol_TVECTOR(i++) );
		for(j=0; j<nyield; j++)delGamma[j] = shapeFAD::val(delGamma_T2[j]) - lambda * shapeFAD::val( Sol_TVECTOR(i++) );
		
		PlasticResidual<REAL, REAL>(N, Np1, delGamma, res, normEpsPErr, 1 /*silent*/);
		
        // resetting invalid equations
		for(i=0; i<nyield; i++)	
			if(validEqs[i] == 0)
				res[i+7] = 0.;
		
		sqrNormResNp1 = pow(sdot(res,res),0.5);		
		
	}while(sqrNormResNp1 < sqrNormResN && lambda >= fMinLambda); // ensuring that the step will be larger than fMinLambda
			                                                     // to avoid wandering within local minima.

	lambda *= k;
	
	if(lambda < 1.)
	{
#ifdef LOG4CXX_PLASTICITY
      {
        std::stringstream sout;
        sout << "*** UpdatePlasticVars *** Line Search indicates lambda = " << lambda << " to ensure residual drop." ;
        //LOGPZ_INFO(logger,sout.str().c_str());
      }
#endif	
	}

	if(!updateVars) return lambda;
	
	for(i=0; i<6; i++) Np1_T2.fEpsP.fData[i] -= T2( lambda*Sol_TVECTOR(i) );
	Np1_T2.fAlpha -= T2( lambda*Sol_TVECTOR(i++) );
	for(j=0; j<YC_t::NYield; j++) delGamma_T2[j] -= T2( lambda*Sol_TVECTOR(i++) );
	
	return lambda;

}

template <class YC_t, class TF_t, class ER_t>
template<class T>
void TPZPlasticStep<YC_t, TF_t, ER_t>::InitializePlasticFAD(
                const TPZPlasticState<REAL> &state,
				const TPZVec<REAL> &delGamma,
                TPZPlasticState<T> &state_T,
				TPZVec<T> &delGamma_T,
				const int nVars)const
{
    int i;
	int nVarsPlastic = 7+YC_t::NYield;
	
    //  copying values
	state.CopyTo(state_T); // also initializing derivatives to null
    for(i=0; i<YC_t::NYield; i++)delGamma_T[i] = delGamma[i];

    // Initialize the partial derivative values
    // the first 6 independent variables are the values of the plastic strains
    for(i = 0; i < 6; i++) state_T.fEpsP.fData[i].diff(i,nVars);
    // the damage variable is the seventh variable
    state_T.fAlpha.diff(6,nVars);
    // the remaining variables are the yield function multipliers
    for(i=7; i<nVarsPlastic; i++) delGamma_T[i-7].diff(i,nVars);
	
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::SetMaterialElasticOrPlastic(int choice)
{
	fMaterialElasticOrPlastic = choice;
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ProcessLoad(const TPZTensor<REAL> &sigma, const EElastoPlastic ep) 
{
	const int nVars = 6;
	REAL resnorm;
	int i, k = 0;
	TPZFNMatrix<nVars * nVars> Dep_mat(nVars, nVars);
	TPZFNMatrix<nVars>         residual_mat(nVars, 1),
                               sigma_mat(nVars, 1),
	                           sol_mat(nVars, 1);
	
	TPZTensor<REAL> epsTotal(fN.fEpsT), EEpsilon;
	
#ifdef LOG4CXX
    {
	std::stringstream sout;
	sout << ">>> ProcessLoad *** Evaluating Sigma to compute the resultant stress for the guess strain - Preparing starting tangent matrix for Newton's scheme";
    sout << "\n sigma << " << sigma;
	//LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
    }
#endif
#ifdef LOG4CXX_PLASTICITY
    {
	std::stringstream sout1, sout2;
	sout1 << ">>> ProcessLoad *** Evaluating Sigma to compute the resultant stress for the guess strain - Preparing starting tangent matrix for Newton's scheme";
    sout2 << "\n sigma << " << sigma;
	//LOGPZ_INFO(logger,sout1.str().c_str());
    //LOGPZ_DEBUG(logger,sout2.str().c_str());
    }
#endif
	
	
//	if(fMaterialElasticOrPlastic==0)
//	{
//		ProcessStrain(epsTotal, EForceElastic);
//		ApplyStrainComputeDep(epsTotal, EEpsilon, Dep_mat);
//		return;
//	}
	
	//cout << "\nstarting ProcessStrain/ComputeDep";
	//cout.flush();
	
	// evaluating the plastic integration, stress tensor and jacobian ATUAL
//	ProcessStrain(epsTotal, EAuto);//, ep);
//	ApplyStrainComputeDep(epsTotal, EEpsilon, Dep_mat);
	
	
    
    ProcessStrainNoSubIncrement(epsTotal, EAuto);
    
//	    ProcessStrain(epsTotal, EAuto);//, ep);TESTE
		//ProcessStrain(epsTotal, ep);
		ComputeDep(EEpsilon, Dep_mat);
	
	//cout << "\nended ProcessStrain/ComputeDep";
	//cout.flush();
	
	resnorm = 0.;
	for(i = 0; i < nVars; i++)residual_mat(i,0) = EEpsilon.fData[i] - sigma.fData[i];
	for(i = 0; i < nVars; i++)resnorm += pow(residual_mat(i,0),2.);
	resnorm = sqrt(resnorm);
	
	while(resnorm > fResTol && k < fMaxNewton)
	{
		k++;
		
        TPZFNMatrix<nVars*nVars> *matc = new TPZFNMatrix<nVars*nVars>(nVars,nVars);
        *matc = Dep_mat;
		

		
        TPZStepSolver st(matc);
        st.SetDirect(ELU);
//		st.SetDirect(ELDLt);
		
//		cout << "\nNewton Method Step " << k << "\nresidual_mat=" << residual_mat << endl << matc;
//		matc->Print("matc= ",cout,EMathematicaInput);
		cout.flush();
		  
		// invert the tangent matrix and put the correction in the Sol variable
        st.Solve(residual_mat,sol_mat,0);

		//cout << "\n solve ended:" << sol_mat;
		//cout.flush();
		
		for(i = 0; i < nVars; i ++)epsTotal.fData[i] -= sol_mat(i,0);
		
		//cout << "\nstarting ProcessStrain/ComputeDep";
		//cout.flush();
		
        // evaluating the plastic integration, stress tensor and jacobianATUAL
   //     ProcessStrain(epsTotal, ep);
	//	ApplyStrainComputeDep(epsTotal, EEpsilon, Dep_mat);
		//cout << "\nended ProcessStrain/ComputeDep";
		//cout.flush();
		
		
		//TESTE
        
        ProcessStrainNoSubIncrement(epsTotal, ep);
//		ProcessStrain(epsTotal,ep);
		ComputeDep(EEpsilon,Dep_mat);
		
		
	//	cout << "\nNewton Method Step " << k << "\nDep=" << Dep_mat << endl << residual_mat;
//		cout.flush();
		
		
        resnorm = 0.;
        for(i = 0; i < nVars; i++)residual_mat(i,0) = EEpsilon.fData[i] - sigma.fData[i];
		for(i = 0; i < nVars; i++)resnorm += pow(residual_mat(i,0),2.);
		resnorm = sqrt(resnorm);
		//cout << "\nresidual = " << resnorm;
		
		#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "*** ProcessLoad *** " << k <<"-th iteration of Newton's scheme with residual = " << resnorm;
			//LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
        }
        #endif
		#ifdef LOG4CXX_PLASTICITY
        {
            std::stringstream sout;
            sout << "*** ProcessLoad *** " << k <<"-th iteration of Newton's scheme with residual = " << resnorm;
            //LOGPZ_INFO(logger,sout.str().c_str());
        }
        #endif
		
		if(k > fMaxNewton)cout << "\n*** ProcessLoad step " << k << " with res= " << resnorm;
	}
	
	#ifdef LOG4CXX_PLASTICITY
    {	
		std::stringstream sout;
		
		if( k > fMaxNewton)
		{
        	sout << "<<< ProcessLoad *** Exiting Method with residual = " << resnorm
				 << " after " << k << " steps.";
			sout << "\n#### Truncated Newton ####. Results are unpredictable";
        	//LOGPZ_WARN(logger,sout.str().c_str());
			LOGPZ_WARN(plasticIntegrLogger,sout.str().c_str());
		}else{
			sout << "<<< ProcessLoad *** Exiting Method with residual = " << resnorm;
        	//LOGPZ_INFO(logger,sout.str().c_str());
			LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
		}
    }
    #endif
	
}

    /**
    return the value of the yield function for the given strain
    */
template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::Phi_Internal(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const
{
    TPZTensor<REAL> sigma;
    REAL A;
	
	TPZPlasticState<REAL> state(fN);
	state.fEpsT = epsTotal;
	
    ComputePlasticVars<REAL>(state, sigma, A);
	
    fYC.Compute(sigma, A, phi, 0);
	
    return;
}

template <class YC_t, class TF_t, class ER_t>
template <int N>
void TPZPlasticStep<YC_t, TF_t, ER_t>::PushPlasticMem(
						const TPZPlasticState<REAL> & state,
                        const REAL & k,
						const REAL & lambda,
						const TPZManVector<REAL,N> & delGamma,
						const TPZVec<int> & validEqs,
						const int forceYield)
{
	TPZPlasticIntegrMem<REAL, YC_t::NYield>  
			Mem(state, k, lambda, delGamma, validEqs, forceYield);
  	fPlasticMem.Push(Mem);
		
}

template <class YC_t, class TF_t, class ER_t>
void TPZPlasticStep<YC_t, TF_t, ER_t>::ApplyLoad_Internal(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal)
{
	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
			 << " Imposed sigma << " << sigma;
		//LOGPZ_INFO(plasticIntegrLogger,sout.str().c_str());
    }
#endif
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
			 << " Imposed sigma << " << sigma;
        LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif
/*
	ProcessLoad(sigma, EAuto);
	int n = fPlasticMem.NElements();
*/
	
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout;
        sout << ">>> ApplyLoad_Internal ***"
			 << " Forcing Elastic behaviour";
        LOGPZ_INFO(logger,sout.str().c_str());
    }
#endif


	ProcessLoad(sigma, EForceElastic);
	int n = fPlasticMem.NElements();

	if(!IsStrainElastic(fPlasticMem[n-1].fPlasticState) )
	{
	#ifdef LOG4CXX_PLASTICITY
	    {
    	    std::stringstream sout;
        	sout << ">>> ApplyLoad_Internal ***"
				 << " Forcing Plastic behaviour - Elastic attempt led to a final plastic state";
    	    LOGPZ_INFO(logger,sout.str().c_str());
	    }
	#endif
		ProcessLoad(sigma, EForcePlastic);
		n = fPlasticMem.NElements();
	}
	
	TPZPlasticStep<YC_t, TF_t, ER_t>::SetState_Internal(fPlasticMem[n-1].fPlasticState);
	
    epsTotal = fN.fEpsT;
	
#ifdef LOG4CXX_PLASTICITY
    {
        std::stringstream sout1, sout2;
        sout1 << "<<< ApplyLoad_Internal ***";
        sout2 << "\nOutput epsTotal << " << epsTotal;
        LOGPZ_INFO(logger,sout1.str().c_str());
        LOGPZ_INFO(logger,sout2.str().c_str());
    }
#endif
	
}

template <class YC_t, class TF_t, class ER_t>
REAL TPZPlasticStep<YC_t, TF_t, ER_t>::IntegrationOverview(TPZVec<REAL> & plastifLen)
{
	int i, j, n = fPlasticMem.NElements();
	
	if(n <= 2) return 0;
	
	plastifLen.Fill(0.);
	REAL plasticLen = 0., deltaK;
	
	for(i = 2; i < n; i++)
	{
		deltaK = fPlasticMem[i].fK - fPlasticMem[i-1].fK;
		for(j = 0; j < YC_t::NYield; j++)if(fPlasticMem[i].fValidEqs[j])plastifLen[j] += deltaK;
		plasticLen += deltaK;
	}
	
	return plasticLen;
}

template <class YC_t, class TF_t, class ER_t>
TPZTensor<REAL> TPZPlasticStep<YC_t, TF_t, ER_t>::gRefDeform;

template class TPZPlasticStep<TPZYCVonMises, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCTresca, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCTrescaRegularized, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCVonMisesCombTresca, TPZThermoForceA, TPZElasticResponse>;

template class TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>;

template class TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCRankine< TPZYCDruckerPrager > , TPZThermoForceA, TPZElasticResponse>;

#include "TPZYCMohrCoulomb.h"
#include  "TPZYCWillamWarnke.h"
#include  "TPZYCModifiedMohrCoulomb.h"

template class TPZPlasticStep<TPZYCMohrCoulomb, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCWillamWarnke, TPZThermoForceA, TPZElasticResponse>;
template class TPZPlasticStep<TPZYCModifiedMohrCoulomb, TPZThermoForceA, TPZElasticResponse>;

template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
				PlasticResidual<REAL, REAL>(TPZPlasticState<REAL> const &, 
											TPZPlasticState<REAL> &, 
											TPZVec<REAL> const&, 
											TPZVec<REAL> &, 
											REAL &, int)const;

template void TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>::
				PlasticResidual<REAL, TFad<14,REAL> >(TPZPlasticState<REAL> const &, 
													  TPZPlasticState< TFad<14,REAL> > &, 
													  TPZVec<TFad<14,REAL> > const &, 
													  TPZVec<TFad<14,REAL> > &, 
													  REAL &, int)const;   

template class TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>;


template void TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>::ComputePlasticVars<double>(TPZPlasticState<double> const&, TPZTensor<double>&, double&) const;

