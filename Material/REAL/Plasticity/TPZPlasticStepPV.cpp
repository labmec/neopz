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

// EM FASE DE DEBUG!!!!!!!
template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)
{
    TPZTensor<REAL>::TPZDecomposed DecompSig,DecompEps; // It may be SigTr or SigPr Decomposition, dependes on the part of this method
	TPZTensor<REAL> sigtr;
	
	//
	TPZTensor<REAL> epsTr,epsPN,epsElaNp1,depsTensor;
	epsPN = fN.fEpsP;
	epsTr = epsTotal;
	epsTr -= epsPN; // Porque soh tem implementado o operator -=
	
	// Compute and Decomposition of SigTrial
	fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
    epsTotal.EigenSystem(DecompEps);
	sigtr.EigenSystem(DecompSig);
	TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
    // Pegando os autovetores
    TPZManVector<TPZFMatrix<REAL>,3> epsegveFromProj(3);
    TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);
    for (int i = 0; i < 3; i++)
	{
		EigenvecMat[i] = DecompSig.fEigenvectors[i];
		epsegveFromProj[i].Resize(3,1);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj[i](j,0) = EigenvecMat[i](j,0);
		}
	}
	for (int i = 0; i < 3; i++) {
		REAL normvec = 0.;
		normvec = NormVecOfMat(epsegveFromProj[i]);
		for (int j = 0; j < 3; j++) {
			epsegveFromProj[i](j,0) /= normvec;
		}
	}
	
	// ReturMap in the principal values
	STATE nextalpha = -6378.;
    TPZFNMatrix<9> GradSigma(3,3,0.);
	fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha,GradSigma);
	fN.fAlpha = nextalpha;
    
    // Aqui calculo minha matriz tangente ------------------------------------
    // Criando matriz tangente
	TPZFNMatrix<36> dSigDe(6,6,0.);
	
	//Montando a matriz tangente
	int kival[] = {0,0,0,1,1,2};
	int kjval[] = {0,1,2,1,2,2};
	REAL G = fER.G();
	REAL lambda = fER.Lambda();
	

	for (int k = 0; k < 6; k++)
	{
		int ki, kj;
		ki = kival[k];
		kj = kjval[k];
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int l = 0; l < 6; l++)
				{
					REAL temp = 2 * G * EigenvecMat[j](kj,ki); // * EigenvecMat[j](j,ki);
					
					if (ki == kj)
					{
						temp += lambda;
					}
					else {
						temp *= 2.;
					}
                    
					temp *= GradSigma(i,j);
					dSigDe(l,k) += temp * DecompSig.fEigenvectors[i][l];
				}///l
			}///j
		}///i		
	}///k
    


    REAL deigensig = 0., deigeneps = 0.;
	
	TPZFNMatrix<6,REAL> factorMat(3,3,0.);
    TPZFNMatrix<9> depsMat(3,3,0.);
    depsTensor = epsTotal;
    depsTensor -= fN.fEpsT;
    depsMat = depsTensor;
    
	for (int i = 0; i < 2; i++) {
		for (int j = i+1; j<3 ; j++) {
			deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
			deigensig = sigprvec[i] - sigprvec[j];
			TPZFNMatrix<9,REAL> tempMat(3,1,0.);
			depsMat.Multiply(epsegveFromProj[i], tempMat);
			REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
			REAL factor = deigensig * deij / deigeneps;
			tempMat.Redim(3, 3);
			tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
			factorMat += tempMat * factor;
		}
	}
	
	// Reconstruction of sigmaprTensor
	DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
	sigma = TPZTensor<REAL>(DecompSig);
	
	fER.ComputeDeformation(sigma,epsElaNp1);
	fN.fEpsT = epsTotal;
	epsPN = epsTotal;
	epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
	fN.fEpsP = epsPN;
}

REAL NormVecOfMat(TPZFNMatrix <9> mat)
{
	REAL norm = 0.;
	for (int i = 0; i < mat.Rows(); i++) {
		norm += mat(i,0) * mat(i,0);
	}
	norm = sqrt(norm);
	return norm;
}

REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	REAL dot = 0.;
	for (int i = 0; i < m1.Rows(); i++) {
		dot += m1(i,0) * m2(i,0);
	}
	return dot;
}

TPZFMatrix<REAL> ProdT(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2)
{
	TPZFMatrix<REAL> mat(3,3,0.);
	for (int i = 0; i < 3; i++) {
		for (int j = 0 ; j < 3; j++) {
			mat(i,j) = m1(i,0) * m2(j,0);
		}
	}
	return mat;
}

template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;



