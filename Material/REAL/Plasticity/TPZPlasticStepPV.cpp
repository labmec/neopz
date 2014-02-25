/**
 * @file
 */


#include "TPZPlasticStepPV.h"
#include "TPZElasticResponse.h"

#include "pzsandlerextPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElasticResponse.h"

#include "pzlog.h"

//#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("pz.material.TPZPlasticStepPV"));
//#endif
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic2"));
#endif

#ifdef LOG4CXX
static LoggerPtr logger2(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

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
//	fYC.ProjectSigma(sigtrvec, fN.fAlpha, sigprvec, nextalpha);
    TPZFNMatrix<9> GradSigma(3,3,0.);
    fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);
	fN.fAlpha = nextalpha;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

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
	TPZTensor<REAL> epsTr,epsPN,epsElaNp1;
	epsPN = fN.fEpsP;
	epsTr = epsTotal;
	epsTr -= epsPN; // Porque soh tem implementado o operator -=
	
	// Compute and Decomposition of SigTrial
	fER.Compute(epsTr, sigtr); // sigma = lambda Tr(E)I + 2 mu E
	epsTr.EigenSystem(DecompEps);

	sigtr.EigenSystem(DecompSig);
	TPZManVector<REAL,3> sigtrvec(DecompSig.fEigenvalues), sigprvec(3,0.);
	
	// Pegando os autovetores
	TPZManVector<TPZFNMatrix<3>,3> epsegveFromProj(3);
	TPZManVector<TPZFNMatrix<9,REAL>, 3 > EigenvecMat(3);     
	for (int i = 0; i < 3; i++)
	{
		EigenvecMat[i] = DecompSig.fEigenvectors[i];
		epsegveFromProj[i].Resize(3,1);
		for	(int k = 0 ; k < 3 ; k++){
			bool IsnoZero = false;
			for (int j = 0; j < 3; j++) {
				epsegveFromProj[i](j,0) = EigenvecMat[i](j,k);
				if (!IsZero(epsegveFromProj[i](j,0))) IsnoZero = true;
			}	
			if (IsnoZero) break;
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
	fYC.ProjectSigmaDep(sigtrvec, fN.fAlpha, sigprvec, nextalpha, GradSigma);
	//GradSigma.Print("Grad");
	fN.fAlpha = nextalpha;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Sig Trial " << sigtrvec << "\nSig Project " << sigprvec << std::endl;
        GradSigma.Print("GradSigma", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

	// Aqui calculo minha matriz tangente ------------------------------------
	// Criando matriz tangente
	TPZFNMatrix<36> dSigDe(6,6,0.);
	
	//Montando a matriz tangente
	int kival[] = {0,0,0,1,1,2};
	int kjval[] = {0,1,2,1,2,2};
	REAL G = fER.G();
	REAL lambda = fER.Lambda();
	int ki, kj;
	for (int k = 0; k < 6; k++)
	{
		
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
	TPZFNMatrix<36> RotCorrection(6,6,0.);
	// Correcao do giro rigido
	for (int i = 0; i < 2; i++) {
		for (int j = i+1; j<3 ; j++) {
			deigeneps = DecompEps.fEigenvalues[i]  - DecompEps.fEigenvalues[j];
			deigensig = sigprvec[i] - sigprvec[j];
			TPZFNMatrix<9,REAL> tempMat(3,3,0.);
			REAL factor = 0.;
			if (!IsZero(deigeneps)) {
				factor = deigensig / deigeneps;
			}
			else {
				factor = fER.G() * ( GradSigma(i,i) - GradSigma(i,j) - GradSigma(j,i) + GradSigma(j,j) );
			}
			tempMat = ProdT(epsegveFromProj[i],epsegveFromProj[j]) + ProdT(epsegveFromProj[j],epsegveFromProj[i]);
			for (int k = 0 ; k < 6 ; k++){
				ki = kival[k];
				kj = kjval[k];
				TPZFNMatrix<9> ColCorr(3,3,0.),ColCorrV(6,1,0.);
				if (ki == kj) {
					REAL tempval = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) );
					ColCorr = tempval * factor * tempMat;
				}
				else {
					ColCorr = (epsegveFromProj[j](ki,0) * epsegveFromProj[i](kj,0) + epsegveFromProj[j](kj,0) * epsegveFromProj[i](ki,0) ) * factor * tempMat;
				}
				ColCorrV = FromMatToVoight(ColCorr);
				for (int l = 0; l < 6; l++) {
					RotCorrection(l,k) += ColCorrV(l,0);
				}
			}
		}
	}

	dSigDe += RotCorrection;
	
#ifdef LOG4CXX
	{
		if(logger->isDebugEnabled())
    {
			std::stringstream str;
			str << "\n**********************MATRIZ TANGENTE**********************" << endl;
			dSigDe.Print("Matriz Tangente:",str);
			str << "\n**********************CORRECAO GIRO**********************" << endl;
			RotCorrection.Print("GiroCorrection",str);
			LOGPZ_DEBUG(logger,str.str())
		}
	}
#endif
	
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
void TPZPlasticStepPV<YC_t, ER_t>::TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv)
{
	TPZTensor<REAL> eps1,eps2, SigmaTemp,Sigma1,Sigma2;
	TPZFNMatrix <36> dSigDe(6,6,0.);
	TPZStack<REAL> coef;
	
	fN.fEpsP.Scale(0.);
	fN.fEpsT.Scale(0.);
	fN.fAlpha = kprev;
	this->ApplyStrainComputeDep(EpsIni,SigmaTemp,dSigDe);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "EpsIni " << EpsIni << "\nSigmaTemp " << SigmaTemp << "\ndSidDe " << dSigDe << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	fN.fEpsP.Scale(0.);
	fN.fEpsT.Scale(0.);
	fN.fAlpha = kprev;
	
	REAL scale = 1.;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int i = 0; i < 6; i++) {
		alphatable[i] *= scale;
	}
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[0];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = EpsIni;
		eps2 = EpsIni;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeSigma(eps1,Sigma1);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		fN.fEpsT = EpsIni;
		this->ApplyStrainComputeSigma(eps2,Sigma2);
		fN.fEpsP.Scale(0.);
		fN.fEpsT.Scale(0.);
		fN.fAlpha = kprev;
		
		TPZFNMatrix <6> deps1(6,1,0.),deps2(6,1,0.);
		TPZFNMatrix <9> depsMat(3,3,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		deps2 = FromMatToVoight(depsMat);
		
		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		dSigDe.Multiply(deps1, tanmult1);
		dSigDe.Multiply(deps2, tanmult2);
		
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
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sigprMat.Print("sigprMat",sout);
            sigpr1Mat.Print("sigpr1Mat",sout);
            tanmult1.Print("tanmult1",sout);
            sigpr2Mat.Print("sigpr2Mat",sout);
            tanmult2.Print("tanmult2",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            error1.Print("error1:",sout);
            error2.Print("error2:",sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		REAL n;
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) );
		coef.push_back(n);
	}
	conv = coef;
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




template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyStrain(const TPZTensor<REAL> &epsTotal)
{

    std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    DebugStop();
    
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &eps)
{

    
//    TPZFNMatrix<36> Dep(6,6,0.);
//    TPZFNMatrix <6> signcalc2(6,1,0.),epsguess(6,1,0.),epsprev(6,1,0.),res(6,1,0.),deltaeps(6,1,0.),sigimposed(6,1,0.);
//    TPZTensor<STATE> epsnt,signcalc,epsn1t;
//    //fER.ComputeDeformation(sigma,epsnt);
//    epsnt=fN.EpsT();
//   // epsnt.XX()=-0.007;
//    TPZPlasticState<STATE> tempstate(fN);
//     ApplyStrainComputeDep(epsnt,signcalc,Dep);
//    //fN=tempstate;
//    CopyFromTensorToFNMatrix(epsnt,epsguess);
//    CopyFromTensorToFNMatrix(signcalc,signcalc2);
//    CopyFromTensorToFNMatrix(sigma,sigimposed);
//    STATE resnormprev,scale;
//    
//    
//    
//
//    for(int i=0;i<6;i++)res(i)=signcalc2[i]-sigimposed[i];
//    //res*=-1;
//    
//    STATE normres=1.,tol=1.e-6;
//    int counter=0;
//    while (normres>tol &&counter<30)
//    {
//        scale=1.;
//        for(int i=0;i<6;i++)resnormprev+=res[i]*res[i];
//        resnormprev=sqrt(resnormprev);
//        for(int i=0;i<6;i++)deltaeps(i)=res[i];
//        Dep.Print("Dep");
//        deltaeps.Print("deltaSig");
//        Dep.Solve_LU(&deltaeps);
//        deltaeps.Print("deltaeps");
//        for(int i=0;i<6;i++)epsprev(i)=epsguess[i];
////        int counter2 =0;
//        do{
//            for(int i=0;i<6;i++)epsguess(i)=epsprev[i]-(deltaeps[i]*scale);
//            CopyFromFNMatrixToTensor(epsguess, epsnt);
//            fN=tempstate;
//            ApplyStrainComputeSigma(epsnt,signcalc);
//            Dep.Print("Dep");
//            CopyFromTensorToFNMatrix(signcalc,signcalc2);
//            for(int i=0;i<6;i++)res(i)=sigimposed[i]+signcalc2[i];
//            //res*=-1;
//            for(int i=0;i<6;i++)normres+=res[i]*res[i];
//            normres = sqrt(normres);
//            scale*=0.5;
////            counter2++;
////            if (counter2%30==0) {
////                scale=1;
////            }
//            
//        }while (normres>=resnormprev);
//        
//        counter++;
//
//        
//    }
    

  //  SetTensionSign(1);
    TPZTensor<REAL> sigma_Internal(sigma);
 //   sigma_Internal *= SignCorrection();
    
    const int nVars = 6;
    REAL resnorm;
    int i, k = 0;
    TPZFNMatrix<nVars*nVars> Dep_mat(nVars,nVars);
    TPZFNMatrix<nVars>  residual_mat(nVars,1),
    sigma_mat(nVars, 1),
    sol_mat(nVars, 1);
    
    TPZPlasticState<STATE> tempstate(fN);
    
    TPZTensor<REAL> epsTotal(fN.fEpsT), EEpsilon;
     
//    ProcessStrain(epsTotal, EAuto);
//    ComputeDep(EEpsilon, Dep_mat);
    ApplyStrainComputeDep(epsTotal, EEpsilon,Dep_mat);

    resnorm = 0.;
    for(i = 0; i < nVars; i++)residual_mat(i,0) = (fabs(EEpsilon.fData[i]) - fabs(sigma_Internal.fData[i]));
    for(i = 0; i < nVars; i++)resnorm += pow(residual_mat(i,0),2.);
    resnorm = sqrt(resnorm);
    STATE restol = 1.e-6;
    REAL resnormprev=0.;
    while(fabs((fabs(resnorm)-fabs(resnormprev))) > restol && k < 30)
    {
        k++;
        
        TPZFNMatrix<nVars*nVars> *matc = new TPZFNMatrix<nVars*nVars>(nVars,nVars);
        *matc = Dep_mat;

        
        TPZStepSolver<REAL> st(matc);
        st.SetDirect(ELU);

        st.Solve(residual_mat,sol_mat,0);

        TPZTensor<REAL> epsTotalPrev(epsTotal);
        REAL scalefactor = 1.;
         resnormprev= resnorm;
        fN=tempstate;
        do {
            for(i = 0; i < nVars; i ++)epsTotal.fData[i] =(fabs(epsTotalPrev.fData[i]) - fabs(scalefactor*sol_mat(i,0)));


//            ProcessStrain(epsTotal,/*o original e ep e nao EAuto */ ep);
//            ComputeDep(EEpsilon, Dep_mat);
              fN=tempstate;
              ApplyStrainComputeDep(epsTotal, EEpsilon,Dep_mat);
              
            resnorm = 0.;
            for(i = 0; i < nVars; i++)residual_mat(i,0) = (fabs(EEpsilon.fData[i]) - fabs(sigma_Internal.fData[i]));
            for(i = 0; i < nVars; i++)resnorm += pow(residual_mat(i,0),2.);
            resnorm = sqrt(resnorm);
            
            scalefactor *= 0.5;
            
        } while (resnorm >= resnormprev);

        //fN=tempstate;
    }
    
    eps=fN.EpsT();


}

template <class YC_t, class ER_t>
TPZPlasticState<STATE>  TPZPlasticStepPV<YC_t, ER_t>::GetState() const
{
    return fN;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const
{

    fYC.Phi(epsTotal,fN.Alpha(),phi);
}


template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetState(const TPZPlasticState<REAL> &state)
{
    fN=state;
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetIntegrTol(REAL integrTol)
{
    //std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    //DebugStop();
}


template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::SetTensionSign(int s)
{
    fInterfaceTensionSign = (s >= 0) ? 1 : -1;
   
}


template <class YC_t, class ER_t>
int TPZPlasticStepPV<YC_t, ER_t>::SignCorrection() const
{
    return fMaterialTensionSign * fInterfaceTensionSign;
}


template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Read(TPZStream &buf)
{
    //std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    //DebugStop();
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::Write(TPZStream &buf)
{
    //std::cout<< " \n this method is not implemented in PlasticStepPV. ";
    //DebugStop();
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromFNMatrixToTensor(TPZFNMatrix<6> FNM,TPZTensor<STATE> &copy)
{
    FNM.Resize(6,1);
    copy.XX()=FNM(0,0);
    copy.XY()=FNM(1,0);
    copy.XZ()=FNM(2,0);
    copy.YY()=FNM(3,0);
    copy.YZ()=FNM(4,0);
    copy.ZZ()=FNM(5,0);
}

template <class YC_t, class ER_t>
void TPZPlasticStepPV<YC_t, ER_t>::CopyFromTensorToFNMatrix(TPZTensor<STATE> tensor,TPZFNMatrix<6> &copy)
{
    
    copy(0,0)=tensor.XX();
    copy(1,0)=tensor.XY();
    copy(2,0)=tensor.XZ();
    copy(3,0)=tensor.YY();
    copy(4,0)=tensor.YZ();
    copy(5,0)=tensor.ZZ();
}

template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;
template class TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>;

/*
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
 
 }
 }
*/