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
		epsegveFromProj[i].Print("eg");
	}
	
	// ReturMap in the principal values
	STATE nextalpha = -6378.;
	TPZFNMatrix<9> GradSigma(3,3,0.);
	fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);
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
	
	// Correcao do giro rigido
	for (int i = 0; i < 2; i++) {
		for (int j = i+1; j<3 ; j++) {
			deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
			deigensig = sigprvec[i] - sigprvec[j];
			TPZFNMatrix<9,REAL> tempMat(3,1,0.);
			depsMat.Multiply(epsegveFromProj[i], tempMat);
			REAL deij = InnerVecOfMat(tempMat,epsegveFromProj[j]);
			REAL factor = 0.;
			if (!IsZero(deigeneps)) {
				factor = deigensig * deij / deigeneps;
			}
			else {
				factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) ) * deij;
			}
			std::cout << "factor = " << factor << std::endl;
			std::cout << "G = " << fER.G() << std::endl;
			GradSigma.Print("GradSigma");
			tempMat.Redim(3, 3);
			tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
			factorMat += tempMat * factor;
			factorMat.Print("factorMat");
		}
	}
	TPZFNMatrix<6> factorV = FromMatToVoight(factorMat);
	
	
	// Achando o primeiro termo nao nulo de dEps (Tenho medo disso estar feio e ter jeito melhor - Nathan)
	long pos = -1;
	REAL depspos = 0.;
	for (int i = 0 ; i < 6 ; i++){
		if(!IsZero(depsTensor[i])) {
			pos = i;
			depspos = depsTensor[i]; 
			break;
		}
	}
	
	if (pos == -1) {
		PZError << "dEps Totalmente Nulo!!" << std::endl;
		DebugStop(); // acho que devo tratar isso diferente. Somente nao fazer nada na matriz. Lembrar de perguntar pro Phil
	}
	
	for (int i = 0; i < 6; i++) {
		dSigDe(i,pos) += factorV(i,0)/depspos; //AQUI NATHAN
	}	
	
	// Reconstruction of sigmaprTensor
	DecompSig.fEigenvalues = sigprvec; // CHANGING THE EIGENVALUES FOR THE ONES OF SIGMAPR
	sigma = TPZTensor<REAL>(DecompSig);
	
	fER.ComputeDeformation(sigma,epsElaNp1);
	fN.fEpsT = epsTotal;
	epsPN = epsTotal;
	epsPN -= epsElaNp1; // Transforma epsPN em epsPNp1
	fN.fEpsP = epsPN;
	Dep = dSigDe;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev)
{
	TPZTensor<REAL> eps1,eps2, SigmaTemp,Sigma1,Sigma2;
	TPZFNMatrix <36> dSigDe1(6,6,0.),dSigDe2(6,6,0.);
	TPZStack<REAL> coef;
	
	this->ApplyStrainComputeSigma(EpsIni,SigmaTemp);
	fN.fEpsP.Scale(0.);
	fN.fEpsT.Scale(0.);
	fN.fAlpha = kprev;
	
	REAL scale = 0.01;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int i = 0; i < 6; i++) {
		alphatable[i] *= scale;
	}
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[ia];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = EpsIni;
		eps2 = EpsIni;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeDep(eps1,Sigma1,dSigDe1);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeDep(eps2,Sigma2,dSigDe2);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		TPZFNMatrix <6> deps1(6,1,0.),deps2(6,1,0.);
		TPZFNMatrix <9> depsMat(3,3,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		deps2 = FromMatToVoight(depsMat);
		
		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		dSigDe1.Multiply(deps1, tanmult1);
		dSigDe2.Multiply(deps2, tanmult2);
		
		for (int i = 0 ; i < 6; i++) {
			tanmult1(i,0) *= alpha1;
			tanmult2(i,0) *= alpha2;
		}
		
		TPZFNMatrix <9> SigMatTemp33(3,3,0.);
		TPZFNMatrix <6> sigprMat(6,1,0.),sigpr1Mat(6,1,0.),sigpr2Mat(6,1,0.);
		SigMatTemp33 = SigmaTemp;
		sigprMat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma1;
		sigpr1Mat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma2;
		sigpr2Mat = FromMatToVoight(SigMatTemp33);
		
		TPZFNMatrix<6> error1(6,1,0.), error2(6,1,0.);
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
		
		REAL n;
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) );
		coef.push_back(n);
	}
	std::cout << "coef = " << coef << std::endl;
}

template <class YC_t, class ER_t>
REAL TPZPlasticStepPV<YC_t, ER_t>::ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat)
{
	REAL norm1, norm2, n;
	norm1 = NormVecOfMat(error1Mat);
	norm2 = NormVecOfMat(error2Mat);
	n = log(norm1/norm2) / log(alpha1/alpha2);
	return n;
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

TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
	TPZFNMatrix <6> voi(6,1,0.);
	int k = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i ; j < 3; j++) {
			voi(k++,0) = mat(i,j);
		}
	}
	return voi;	
}

template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;



