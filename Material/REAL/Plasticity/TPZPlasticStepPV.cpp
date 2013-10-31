/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "pzsandlerextPV.h"
#include "TPZElasticResponse.h"

#include "pzlog.h"

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
{
	TPZTensor<REAL>::TPZDecomposed DecompSig; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
	TPZTensor<REAL> sigtr;
	
	//
	TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
	epsPN = fN.fEpsP;
	epsTr = epsTotal;
	epsTr -= epsPN; // Porque soh tem implementado o operator -=
	
	// Compute and Decomposition of SigTrial
	fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
	sigtr.EigenSystem(DecompSig);
	TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
	
	// ReturMap in the principal values
	STATE nextalpha = -6378.;
	fYC.ProjectSigma(sigtrvec, fN.fAlpha, sigprvec, nextalpha);
	fN.fAlpha = nextalpha;
	
	
	// Reconstruction of sigmaprTensor
	DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
	sigma = TPZTensor<REAL>(DecompSig);
	
	fER.ComputeDeformation(sigma,epsElaNp1);
	fN.fEpsT = epsTotal;
	epsPN = epsTotal;
	epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
	fN.fEpsP = epsPN; 
	
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
{
	
}


template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;



