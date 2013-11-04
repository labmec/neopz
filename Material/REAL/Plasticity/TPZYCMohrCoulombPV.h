/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Nathan Shauer on 5/4/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZYCMOHRCOULOMBPV_H
#define TPZYCMOHRCOULOMBPV_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticState.h"
#include "TPZElasticResponse.h"

#ifdef LOG4CXX
static LoggerPtr loggerMohrCoulombPV(Logger::getLogger("pz.plasticity.mohrcoulombpv"));
#endif

class TPZYCMohrCoulombPV
{
	REAL fPhi;
	REAL fPsi;
	REAL coesion;
	
public:
	
	/// Internal structure to represent the plastic memory (plastic deformation and damage)
	struct TPlasticState
	{
		TPlasticState() : fEpsPlastic(), fEpsPlasticBar(0.)
		{
			
		}
		
		TPlasticState(const TPlasticState &copy) : fEpsPlastic(copy.fEpsPlastic), fEpsPlasticBar(copy.fEpsPlasticBar)
		{
			
		}
		
		TPlasticState &operator=(const TPlasticState &copy)
		{
			fEpsPlastic = copy.fEpsPlastic;
			fEpsPlasticBar = copy.fEpsPlasticBar;
			return *this;
		}
		
		void Print(std::ostream &out) const
		{
			out << "Plastic Deformation tensor ";
			fEpsPlastic.Print(out);
			out << "Acumulated plastic deformation " << fEpsPlasticBar << std::endl;
		}
		
		/// plastic deformation tensor
		TPZTensor<REAL> fEpsPlastic;
		
		/// accumulated damage
		REAL fEpsPlasticBar;
	};
	
	
	/// structure which contains the decision tree of the return map
	// we can only expect a consistent tangent matrix if the decision tree remains the same
	struct TComputeSequence
	{
		TComputeSequence() : fWhichPlane(ENoPlane), fGamma(0)
		{
			
		}
		
		TComputeSequence(const TComputeSequence &copy) : fWhichPlane(copy.fWhichPlane), fGamma(copy.fGamma)
		{
			
		}
		
		TComputeSequence &operator=(const TComputeSequence &copy)
		{
			fWhichPlane = copy.fWhichPlane;
			fGamma = copy.fGamma;
			return *this;
		}
		
		enum MPlane {ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EHydroStatic };
		
		MPlane fWhichPlane;
		
		TPZManVector<REAL> fGamma;
		
		
	};
	
protected:
	
	/// information of the plastic state of the material point
	TPZYCMohrCoulombPV::TPlasticState fState;
	
	
public:
	
	TPZYCMohrCoulombPV() : fPhi(M_PI/9.),fPsi(M_PI/9.), coesion(9.35) //NATHAN: Depois terei que setar esse valores via construtor
	{
		
	}
	/*
	 REAL Lambda()
	 {
	 return fPoisson*fYoung/(1.+fPoisson)*(1.-2.*fPoisson); // NATHAN: acho que isso está errado
	 }
	 REAL Mu()
	 {
	 return fYoung/(2.*(1.+fPoisson));
	 }
	 REAL G()
	 {
	 return fYoung/(2.*(1+fPoisson));
	 }
	 REAL K()
	 {
	 return fYoung/(3.*(1.-2.*fPoisson));
	 }
	 */
	
	void Print(std::ostream &out) const
	{
		out << "TPZYCMohrCoulombPV\n";
		fState.Print(out);
	}
	
	/// the hardening function and its derivative
	// aqui eh a funcao da coesao em funcao de epsp e sua derivada em funcao de epsp
	// essas funcoes parecem ser quaisquer. Sera que devo passar com forcingfunction depois?????
	template<class T>
	void PlasticityFunction(T epsp, T &sigmay, T &H) const
	{
		//sigmay = T(15.)+(T(2571.43)-T(2.95238e6)*epsp)*(T(-0.0035)+epsp);
		//H = T(12904.8)-T(5.90476e6)*epsp;
		//sigmay=20.;
		// H=0.;
		sigmay = T(15.)+(T(2571.43))*(T(-0.0035)) + T(12904.8)*epsp*epsp; // c(epsp)
		H = T(12904.8*2.)*epsp; // dc(epsp)/depsp
	}
	
	/// a piecewise linear hardening function
	template<class T>
	void PieceWise(T epbar, T &m, T & fx)const // NATHAN: acho que não teremos isso, certo?
	{
		
		ifstream file("curvadehardening.txt");
		TPZFMatrix<REAL> mat;
		int sz;
		file >>sz;
		mat.Resize(sz,2);
		
		for(int i=0;i<sz-1;i++)
		{
			
			file >> mat(i,0);
			file >> mat(i,1);
		}
		
		//cout << "\n mat = "<< mat <<endl;
		for(int i=0;i<sz-1;i++)
		{
			if(epbar >= mat(i,0) && epbar <=  mat(i+1,0))
			{
				REAL x0 = mat(i,0);
				REAL y0 = mat(i,1);
				REAL x = mat(i+1,0);
				REAL y = mat(i+1,1);
				m = (y - y0)/(x - x0);
				fx=m*(epbar-x0)+y0;
				// return;
			}
			
		}
	}
	
	
	
	/// sigma = lambda Tr(E)I + 2 mu E
	template<class T>
	TPZVec<T> SigmaElastPV(const TPZVec<T> &deform, TPZElasticResponse &ER) 
	{
		T trace = deform[0]+deform[1]+deform[2];
		TPZVec<T> sigma(3,0.);
		sigma = trace * ER.Lambda() + 2 * ER.G() * deform[0];
		sigma = trace * ER.Lambda() + 2 * ER.G() * deform[1];
		sigma = trace * ER.Lambda() + 2 * ER.G() * deform[2];
		
		return sigma;
	}
	
	
	
	// Metodo que calcula SigmaTrial. Lembrando que que epsTr = epsTr + dEps
	// Considera como se fosse elastico, e caso nao seja, nos returnMaps ira subtrair a parcela correta
	// epstotal eh dado em valores principais!!
	template<class T>
	TPZVec<T> SigmaTrial(const TPZVec<T> &epstotal, TPZElasticResponse &ER)
	{
		TPZVec<T> epsela(epstotal);
		for(int i=0; i<3; i++)
		{
			epsela[i] -= T(fState.fEpsPlastic[i]);
		}
		TPZVec<T> sigma;
		sigma = SigmaElastPV(epsela,ER);
		
		
		//typename TPZTensor<T>::TPZDecomposed sigma_trial;
		//sigma.EigenSystem(sigma_trial);
		
		
#ifdef LOG4CXX
		if (loggerMohrCoulombPV->isDebugEnabled()) {
			std::stringstream sout;
			sout << "Input epst ";
			sout << epstotal << std::endl;
			sout << "Input sigmatrial in Principal Values form\n" << std::endl;
			sout << sigma << std::endl;
			LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())
		}
#endif
		
		return sigma;
	}
	
	
	/*
	 
	 void ComputeSigmaTangent(TPZTensor<REAL> &epstotal, TPZTensor<REAL> &sigma, TPZFNMatrix<36,REAL> &tangent, const TComputeSequence &memory) //NATHAN: Preciso?
	 {
	 typedef TFad<6,REAL> fadtype;
	 TPZTensor<fadtype> epstotalFAD, sigmaElastFAD, sigmaFAD;
	 for (int i=0; i<6; i++) {
	 epstotalFAD[i].val() = epstotal[i];
	 epstotalFAD[i].fastAccessDx(i) = 1.;
	 }
	 sigmaElastFAD = SigmaElast(epstotalFAD);
	 
	 switch (memory.fWhichPlane) {
	 case TComputeSequence::ENoPlane:
	 DebugStop();
	 break;
	 case TComputeSequence::EElastic:
	 sigmaFAD = sigmaElastFAD;
	 break;
	 case TComputeSequence::EMainPlane:
	 case TComputeSequence::ELeftEdge:
	 case TComputeSequence::ERightEdge:
	 {
	 TPZTensor<fadtype>::TPZDecomposed sigma_trial = SigmaTrial(epstotalFAD);
	 TPZTensor<fadtype>::TPZDecomposed sigma_projected;
	 TComputeSequence locmem(memory);
	 switch (memory.fWhichPlane) {
	 case TComputeSequence::EMainPlane:
	 //ReturnMapPlane<fadtype>(sigma_trial, sigma_projected, locmem);                        
	 break;
	 case TComputeSequence::ELeftEdge:
	 ReturnMapLeftEdge<fadtype>(sigma_trial, sigma_projected, locmem);
	 break;
	 case TComputeSequence::ERightEdge:
	 ReturnMapRightEdge<fadtype>(sigma_trial, sigma_projected, locmem);
	 break;
	 default:
	 DebugStop();
	 break;
	 }
	 sigmaFAD = TPZTensor<fadtype>(sigma_projected);
	 }   
	 default:
	 break;
	 }
	 for (int i=0; i<6; i++) {
	 sigma[i] = sigmaFAD[i].val();
	 for (int j=0; j<6; j++) {
	 tangent(i,j) = sigmaFAD[i].fastAccessDx(j);
	 }
	 }
	 }
	 
	 */
	/*
	 void CommitDeformation(TPZVec<REAL> &epstotal, TComputeSequence &memory) //NATHAN: Preciso?
	 {
	 
	 const REAL cosphi = cos(fPhi);
	 TPZVec<REAL> sigma;
	 switch (memory.fWhichPlane) {
	 case TComputeSequence::EElastic:
	 {
	 TPZVec<REAL> epsela(epstotal);
	 for(int i=0; i<3; i++)
	 {
	 epsela[i] -= fState.fEpsPlastic[i];
	 }
	 //epsela -= fState.fEpsPlastic;
	 sigma = SigmaElast(epsela);
	 //sigma = SigmaElast(epstotal);
	 }
	 break;
	 case TComputeSequence::EMainPlane:
	 case TComputeSequence::ELeftEdge:
	 case TComputeSequence::ERightEdge:
	 {
	 TPZVec<REAL> sigma_trial = SigmaTrial(epstotal);
	 TPZVec<REAL> sigma_projected;
	 TComputeSequence locmem(memory);
	 switch (memory.fWhichPlane) {
	 case TComputeSequence::EMainPlane:
	 ReturnMapPlane<REAL>(sigma_trial, sigma_projected, locmem);
	 fState.fEpsPlasticBar+=(locmem.fGamma[0]*2.*cosphi);
	 break;
	 case TComputeSequence::ELeftEdge:
	 ReturnMapLeftEdge<REAL>(sigma_trial, sigma_projected, locmem);
	 fState.fEpsPlasticBar+=(locmem.fGamma[0]+locmem.fGamma[1])*2.*cosphi;
	 break;
	 case TComputeSequence::ERightEdge:
	 ReturnMapRightEdge<REAL>(sigma_trial, sigma_projected, locmem);
	 fState.fEpsPlasticBar+=(locmem.fGamma[0]+locmem.fGamma[1])*2.*cosphi;
	 break;
	 default:
	 DebugStop();
	 break;
	 }
	 sigma = sigma_projected;
	 break;
	 }
	 default:
	 {
	 break;
	 }
	 }
	 
	 
	 
	 */        
	
	/*REAL tempval;
	 TPZTensor<REAL> StressDeviatoric,P,I,epsplastic(epstotal),epselastic;
	 
	 P.XX()=1; I.XX()=1;
	 P.YY()=1; I.YY()=1;
	 P.ZZ()=1; I.ZZ()=1;
	 
	 tempval=(sigma.I1()/3.)*(1./(3.*K()));
	 P.Multiply(tempval,1);
	 
	 cout << " \n P = "<< P <<endl;
	 sigma.S(StressDeviatoric);
	 
	 StressDeviatoric*=1./(2.*G());
	 cout << " \n S = "<< StressDeviatoric <<endl;
	 
	 StressDeviatoric.Add(P,1);
	 epselastic=StressDeviatoric;
	 
	 epsplastic-=epselastic;*/
	
	//Cálculo de e11ela, e22ela e e33ela
	
	
	/*
	 
	 TPZVec<REAL> epsplastic(3,0.), epselastic(3,0.);
	 REAL cte1, cte2, den;
	 den = 6 * Lambda() * Mu() + 4 * Mu() * Mu();
	 cte1 = Lambda()/den;
	 cte2 = Mu()/den;
	 
	 epselastic[0] = cte1 * (+ 2 * sigma[0] - sigma[1] - sigma[2] ) + cte2 * 2 * sigma[0];
	 epselastic[1] = cte1 * (- sigma[0] + 2 * sigma[1] - sigma[2] ) + cte2 * 2 * sigma[1];
	 epselastic[2] = cte1 * (- sigma[0] - sigma[1] + 2 * sigma[2] ) + cte2 * 2 * sigma[2];
	 
	 
	 for (int i = 0 ; i < 3; i++) 
	 {
	 epsplastic[i] = epstotal[i] - epselastic[i];
	 }
	 
	 fState.fEpsPlastic = epsplastic;
	 }
	 
	 template<class T>
	 TComputeSequence ComputeSigma(TPZVec<T> &epstotal, TPZVec<T> &sigma, bool commitdefor)
	 {
	 
	 TPZVec<T> sigma_trial = SigmaTrial(epstotal);
	 TComputeSequence memory;
	 T phi = PhiPlane<T>(sigma_trial);
	 if (shapeFAD::val(phi) <= 0.) {
	 memory.fWhichPlane = TComputeSequence::EElastic;
	 memory.fGamma.Resize(0);
	 sigma = sigma_trial;
	 //state.fEpsT = epstotal;
	 return memory;
	 }
	 TPZVec<T> sigma_projected;
	 memory.fGamma.Resize(1);
	 memory.fGamma[0] = 0.;
	 if (ReturnMapPlane<T>(sigma_trial, sigma_projected, memory)) {
	 sigma = sigma_projected;
	 memory.fWhichPlane = TComputeSequence::EMainPlane;
	 }
	 else {
	 memory.fGamma.Resize(2);
	 memory.fGamma[0] = 0.;
	 memory.fGamma[1] = 0.;
	 
	 const REAL sinpsi = sin(fPsi);
	 REAL val = (1-sinpsi)*shapeFAD::val(sigma_trial[0])-2.*shapeFAD::val(sigma_trial[2])+(1+sinpsi)*shapeFAD::val(sigma_trial[1]);
	 if (val > 0.) {
	 ReturnMapRightEdge<T>(sigma_trial, sigma_projected, memory);
	 memory.fWhichPlane = TComputeSequence::ERightEdge;
	 }
	 else {
	 ReturnMapLeftEdge<T>(sigma_trial, sigma_projected, memory);
	 memory.fWhichPlane = TComputeSequence::ELeftEdge;
	 }
	 #ifdef LOG4CXX
	 {
	 std::stringstream sout;
	 sout << "After the map to the edge, sigma_projected :\n";
	 // Não está imprimindo nada
	 LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())
	 }
	 #endif
	 sigma = sigma_projected;
	 }
	 if(commitdefor == true)
	 {
	 CommitDeformation(epstotal,memory);
	 }	
	 return memory;
	 }
	 
	 */ 
	
	/// Calcula o valor da funcao criteiro de plastificacao
	template<class T>
	T PhiPlane(TPZVec<T> &sigma) const
	{
		const REAL sinphi = sin(fPhi);
		const REAL cosphi = cos(fPhi);
		T sigmay,H;
		PlasticityFunction(T(fState.fEpsPlasticBar),sigmay, H);
		
		return sigma[0]-sigma[2]+(sigma[0]+sigma[2])*sinphi-2.*sigmay*cosphi;
	}
	
	
	template<class T>
	bool ReturnMapPlane(TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected, 
											TComputeSequence &memory, TPZElasticResponse &ER)
	{
		sigma_projected = sigma_trial;
		TPZManVector<T,3> eigenvalues = sigma_projected;
		//		TPZManVector<T,3> eigenvalues2 = sigma_projected;
		//        TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;
		const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL sinphi2 = sinphi*sinphi;
		const REAL cosphi2 = 1.-sinphi2;
		const REAL constA = 4.* ER.G() *(1.+ sinphi*sinpsi/3.) + 4.*ER.K() * sinphi*sinpsi;
		T sigmay,H;
		T epsbar = T(fState.fEpsPlasticBar+memory.fGamma[0]*2.*cosphi);//diogo aqui
		PlasticityFunction(epsbar,sigmay, H);
		T phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi;
		T gamma = memory.fGamma[0];
		REAL phival = shapeFAD::val(phi);
		REAL tolerance = 1.e-8;
		do {
			T denom = -constA- T(4.*cosphi2)*H;
			//            T d = T(-4.*G()*(1.+sinphi*sinpsi/3.)-4.*K()*sinphi*sinpsi)-T(4.*cosphi2)*H;
			T deriv_gamma = -phi/denom;
			gamma += deriv_gamma;
			epsbar = T(fState.fEpsPlasticBar)+gamma*T(2.*cosphi); //errado esta inicializando toda vez. diogo aqui
			PlasticityFunction(epsbar, sigmay, H);
			if (shapeFAD::val(H) < 0.) {
				DebugStop();
			}
			
			//			eigenvalues2 = eigenvalues; Testes de validacao: phi == phi2
			//			eigenvalues2[0] -= T(2.*G()*(1+sinpsi/3.)+2.*K()*sinpsi)*gamma;
			//			eigenvalues2[1] += T((4.*G()/3. - K()*2.)*sinpsi)*gamma;
			//			eigenvalues2[2] += T(2.*G()*(1-sinpsi/3.)-2.*K()*sinpsi)*gamma;
			
			phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi - constA*gamma;
			//			T phi2 = eigenvalues2[0]-eigenvalues2[2]+(eigenvalues2[0]+eigenvalues2[2])*sinphi-2.*sigmay*cosphi;
			phival = shapeFAD::val(phi);
			//			std::cout << "phi1 = " << phi << "\t phi2 = " << phi2 << std::endl;
			
		} while (abs(phival) > tolerance);
		
		TPZFNMatrix<9>tangAnal(3,3);
		REAL epsbarREAL = shapeFAD::val(epsbar);
		this->ComputePlaneTangent(ER,tangAnal,epsbarREAL);
		tangAnal.Print("tangAnal");
		
		memory.fGamma[0] = shapeFAD::val(gamma);
		eigenvalues[0] -= T(2.*ER.G()*(1+sinpsi/3.)+2.*ER.K()*sinpsi)*gamma;
		eigenvalues[1] += T((4.*ER.G()/3. - ER.K()*2.)*sinpsi)*gamma;
		eigenvalues[2] += T(2.*ER.G()*(1-sinpsi/3.)-2.*ER.K()*sinpsi)*gamma;
		sigma_projected = eigenvalues;
#ifdef DEBUG
		phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi;
#endif
		return (shapeFAD::val(eigenvalues[0])>shapeFAD::val(eigenvalues[1]) && shapeFAD::val(eigenvalues[1]) > shapeFAD::val(eigenvalues[2]));
	}
	
	void ComputePlaneTangent(TPZElasticResponse &ER, TPZMatrix<REAL> &tang, REAL &epsbarp)
	{
		const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL cosphi2 = cosphi*cosphi;
		const REAL G = ER.G(), K = ER.K();
		const REAL c1 = 2.*G*(1.+1./3.*sinpsi) + 2.*K*sinpsi;
		const REAL c2 = (4.*G/3.-2.*K)*sinpsi;
		const REAL c3 = 2.*G*(1.-1./3.*sinpsi) - 2.*K*sinpsi;
		const REAL constA = 4.* G *(1.+ sinphi*sinpsi/3.) + 4.* K * sinphi*sinpsi;
		REAL epsbar = epsbarp;
		REAL c, H;
		PlasticityFunction(epsbar, c, H);
		const REAL denom = constA + 4 * cosphi2*H;
		const REAL dGds1 = (1+sinphi)/denom; // Derivate of gamma with respect to Sigma1tr
		const REAL dGds2 = 0.; // Only created to remember that PHIfunc doesnt depend in it
		const REAL dGds3 = (-1+sinphi)/denom;
		tang.Redim(3, 3);
		
		// First column
		tang(0,0) = 1.-c1*dGds1;
		tang(1,0) = c2*dGds1;
		tang(2,0) = c3*dGds1;
		
		// Second column
		tang(1,1) = 1.; // The others are 0
		
		// Third column
		tang(0,2) = -c1*dGds3;
		tang(1,2) = c2*dGds3;
		tang(2,2) = 1.+c3*dGds3;
	}
	 
	
	template<class T>
	bool ReturnMapLeftEdge(TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
												 TComputeSequence &memory, TPZElasticResponse &ER)
	{
		
		sigma_projected = sigma_trial;
		TPZManVector<T,3> eigenvalues = sigma_projected;
		//        TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;        
		const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL sinphi2 = sinphi*sinphi;
		const REAL cosphi2 = 1.-sinphi2;
		TPZManVector<T,2> gamma(2,0.),phi(2,0.),sigma_bar(2,0.),ab(2,0.);
		gamma[0] = memory.fGamma[0];
		gamma[1] = memory.fGamma[1];
		TPZManVector<REAL,2> phival(2,0.);
		TPZFNMatrix<4,T> d(2,2,0.), dinverse(2,2,0.);
		sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
		sigma_bar[1] = eigenvalues[1]-eigenvalues[2]+(eigenvalues[1]+eigenvalues[2])*T(sinphi);
		T sigmay,H;
		T epsbar = T(fState.fEpsPlasticBar) + (gamma[0]+gamma[1])*T(2.*cosphi);//diogo aqui
		PlasticityFunction(epsbar,sigmay, H);
		phi[0] = sigma_bar[0] - T(2.*cosphi)*sigmay;
		phi[1] = sigma_bar[1] - T(2.*cosphi)*sigmay;
		ab[0] = T(4.*ER.G()*(1+sinphi*sinpsi/3.)+4.*ER.K()*sinphi*sinpsi);
		ab[1] = T(2.*ER.G()*(1.-sinphi-sinpsi-sinphi*sinpsi/3.)+4.*ER.K()*sinphi*sinpsi);
		T residual =1;
		REAL tolerance = 1.e-8;
		do {
			d(0,0) = -ab[0]-T(4.*cosphi2)*H;
			d(1,0) = -ab[1]-T(4.*cosphi2)*H;
			d(0,1) = -ab[1]-T(4.*cosphi2)*H;
			d(1,1) = -ab[0]-T(4.*cosphi2)*H;
			T detd = d(0,0)*d(1,1)-d(0,1)*d(1,0);
			dinverse(0,0) = d(1,1)/detd;
			dinverse(1,0) = -d(1,0)/detd;
			dinverse(0,1) = -d(0,1)/detd;
			dinverse(1,1) = d(0,0)/detd;
			gamma[0] -= (dinverse(0,0)*phi[0]+dinverse(0,1)*phi[1]);
			gamma[1] -= (dinverse(1,0)*phi[0]+dinverse(1,1)*phi[1]);
			//T epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi); diogo aqui
			epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
			PlasticityFunction(epsbar, sigmay, H);
			phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T(2.*cosphi)*sigmay;
			phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[1] - T(2.*cosphi)*sigmay;
			phival[0] = shapeFAD::val(phi[0]);
			phival[1] = shapeFAD::val(phi[1]);
			residual=(fabs(phival[0])+fabs(phival[1]));//aqui diogo
			//} while (abs(phival[0]) > tolerance || abs(phival[1]) > tolerance);//aqui diogo
		}while (residual>tolerance);//aqui diogo
		//        eigenvalues[0] -= T(2.*G()*(1+sinpsi/3.)+2.*K()*sinpsi)*gamma;
		//        eigenvalues[1] += T((4.*G()/3. - K()*2.)*sinpsi)*gamma;
		//        eigenvalues[2] += T(2.*G()*(1-sinpsi/3.)-2.*K()*sinpsi)*gamma;
		
		memory.fGamma[0] = shapeFAD::val(gamma[0]);
		memory.fGamma[1] = shapeFAD::val(gamma[1]);
		eigenvalues[0] -= T(2.*ER.G()*(1+sinpsi/3.)+2.*ER.K()*sinpsi)*gamma[0]+T((4.*ER.G()/3.-2.*ER.K())*sinpsi)*gamma[1];
		eigenvalues[1] += T((4.*ER.G()/3.- ER.K()*2.)*sinpsi)*gamma[0]-T(2.*ER.G()*(1.+sinpsi/3.)+2.*ER.K()*sinpsi)*gamma[1];
		eigenvalues[2] -= T(2.*ER.G()*(1-sinpsi/3.)-2.*ER.K()*sinpsi)*(gamma[0]+gamma[1]);
		sigma_projected = eigenvalues;
		
		return (shapeFAD::val(eigenvalues[0])>shapeFAD::val(eigenvalues[1]) && shapeFAD::val(eigenvalues[1]) > shapeFAD::val(eigenvalues[2]));
	}
	
	void ComputeLeftEdgeTangent(TPZElasticResponse &ER, TPZMatrix<REAL> &tang, REAL &epsbarp)
	{
		//IMPLEMENTAR!!!!
		/*const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL cosphi2 = cosphi*cosphi;
		const REAL G = ER.G(), K = ER.K();
		const REAL c1 = 2.*G*(1.+1./3.*sinpsi) + 2.*K*sinpsi;
		const REAL c2 = (4.*G/3.-2.*K)*sinpsi;
		const REAL c3 = 2.*G*(1.-1./3.*sinpsi) - 2.*K*sinpsi;
		const REAL constA = 4.* G *(1.+ sinphi*sinpsi/3.) + 4.* K * sinphi*sinpsi;
		REAL epsbar = epsbarp;
		REAL c, H;
		PlasticityFunction(epsbar, c, H);
		const REAL denom = constA + 4 * cosphi2*H;
		const REAL dGds1 = (1+sinphi)/denom; // Derivate of gamma with respect to Sigma1tr
		const REAL dGds2 = 0.; // Only created to remember that PHIfunc doesnt depend in it
		const REAL dGds3 = (-1+sinphi)/denom;
		tang.Redim(3, 3);
		
		// First column
		tang(0,0) = 1-c1*dGds1;
		tang(1,0) = c2*dGds1;
		tang(2,0) = c3*dGds1;
		
		// Second column
		tang(1,1) = 1; // The others are 0
		
		// Third column
		tang(0,2) = -c1*dGds3;
		tang(1,2) = c2*dGds3;
		tang(2,2) = 1+c3*dGds3;*/
	}
	
	
	template<class T>
	bool ReturnMapRightEdge(TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
													TComputeSequence &memory, TPZElasticResponse &ER)
	{
		sigma_projected = sigma_trial;
		TPZManVector<T,3> eigenvalues = sigma_projected;
		//      TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;
		const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL sinphi2 = sinphi*sinphi;
		const REAL cosphi2 = 1.-sinphi2;
		const REAL KV = ER.K();
		const REAL GV = ER.G();
		TPZManVector<T,2> gamma(2,0.),phi(2,0.),sigma_bar(2,0.),ab(2,0.);
		gamma[0] = memory.fGamma[0];
		gamma[1] = memory.fGamma[1];
		TPZManVector<REAL,2> phival(2,0.);
		TPZManVector<TPZManVector<T,2>,2> d(2),dinverse(2);
		for (int i = 0; i < 2; i++) {
			d[i].Resize(2,0.);
			dinverse[i].Resize(2,0.);
		}
		//TPZFNMatrix<4,T> d(2,2,0.), dinverse(2,2,0.);
		sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
		sigma_bar[1] = eigenvalues[0]-eigenvalues[1]+(eigenvalues[0]+eigenvalues[1])*T(sinphi);
		T sigmay,H;
		T epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
		PlasticityFunction(epsbar,sigmay, H);
		phi[0] = sigma_bar[0] - T(2.*cosphi)*sigmay;
		phi[1] = sigma_bar[1] - T(2.*cosphi)*sigmay;
		ab[0] = T(4.*GV*(1+sinphi*sinpsi/3.)+4.*KV*sinphi*sinpsi);
		ab[1] = T(2.*GV*(1.+sinphi+sinpsi-sinphi*sinpsi/3.)+4.*KV*sinphi*sinpsi);
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "phi = " << phi << std::endl;
			LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())
		}
#endif
		
		REAL tolerance = 1.e-8;
		int iter = 0;
		T residual =1;
		do {
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "epsbar = " << epsbar << std::endl;
				sout << "sigmay = " << sigmay << std::endl;
				sout << "H = " << H << std::endl;
				LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())
			}
#endif
			d[0][0] = -ab[0]-T(4.*cosphi2)*H;
			d[1][0] = -ab[1]-T(4.*cosphi2)*H;
			d[0][1] = -ab[1]-T(4.*cosphi2)*H;
			d[1][1] = -ab[0]-T(4.*cosphi2)*H;
			T detd = d[0][0]*d[1][1]-d[0][1]*d[1][0];
			dinverse[0][0] = d[1][1]/detd;
			dinverse[1][0] = -d[1][0]/detd;
			dinverse[0][1]= -d[0][1]/detd;
			dinverse[1][1] = d[0][0]/detd;
			gamma[0] -= (dinverse[0][0]*phi[0]+dinverse[0][1]*phi[1]);
			gamma[1] -= (dinverse[1][0]*phi[0]+dinverse[1][1]*phi[1]);
			epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
			PlasticityFunction(epsbar, sigmay, H);
			if (shapeFAD::val(H) < 0.) {
				DebugStop();
			}
			iter++;
			phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T(2.*cosphi)*sigmay;
			phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[1] - T(2.*cosphi)*sigmay;
			phival[0] = shapeFAD::val(phi[0]);
			phival[1] = shapeFAD::val(phi[1]);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "iter = " << iter << " phi = " << phival << std::endl;
				LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())
			}
#endif
			residual=(fabs(phival[0])+fabs(phival[1]));//aqui diogo
			cout << "\n residula = "<< residual << endl;
			//} while (abs(phival[0]) > tolerance || abs(phival[1]) > tolerance);//aqui diogo
		}while (residual>tolerance);//aqui diogo
		
		//        eigenvalues[0] -= T(2.*GV*(1+sinpsi/3.)+2.*KV*sinpsi)*gamma;
		//        eigenvalues[1] += T((4.*GV/3. - KV*2.)*sinpsi)*gamma;
		//        eigenvalues[2] += T(2.*GV*(1-sinpsi/3.)-2.*KV*sinpsi)*gamma;
		// epsbar = T(state.fAlpha)+(gamma[0]+gamma[1])*T(2.*cosphi);
		
		
		memory.fGamma[0] = shapeFAD::val(gamma[0]);
		memory.fGamma[1] = shapeFAD::val(gamma[1]);
#ifdef LOG4CXX
		{
			/*std::stringstream sout;
			sout << "gamma = " << gamma << std::endl;
			sout << "phival = " << phival << std::endl;
			sout << "ab = " << ab << std::endl;
			sout << "sigma_bar = " << sigma_bar << std::endl;
			d.Print("Jacobian",sout);
			dinverse.Print("Inverse Jacobian",sout);
			sout << "epsbar = " << epsbar << std::endl;
			LOGPZ_DEBUG(loggerMohrCoulombPV, sout.str())*/
		}
#endif
		eigenvalues[0] -= T(2.*GV*(1+sinpsi/3.)+2.*KV*sinpsi)*(gamma[0]+gamma[1]);
		eigenvalues[1] += T((4.*GV/3.- KV*2.)*sinpsi)*gamma[0]+T(2.*GV*(1.-sinpsi/3.)-2.*KV*sinpsi)*gamma[1];
		eigenvalues[2] += T(2.*GV*(1-sinpsi/3.)-2.*KV*sinpsi)*gamma[0]+T((4.*GV/3.-2.*KV)*sinpsi)*gamma[1];
		sigma_projected = eigenvalues;
		
		TPZFNMatrix<9>tangAnal(3,3);
		REAL epsbarREAL = shapeFAD::val(epsbar);
		this->ComputeRightEdgeTangent(ER,tangAnal,epsbarREAL);
		tangAnal.Print("tangAnal");
		
		return (shapeFAD::val(eigenvalues[0])>shapeFAD::val(eigenvalues[1]) && shapeFAD::val(eigenvalues[1]) > shapeFAD::val(eigenvalues[2]));    
	 
	}
	
	void ComputeRightEdgeTangent(TPZElasticResponse &ER, TPZMatrix<REAL> &tang, REAL &epsbarp)
	{
		
		const REAL sinphi = sin(fPhi);
		const REAL sinpsi = sin(fPsi);
		const REAL cosphi = cos(fPhi);
		const REAL cosphi2 = cosphi*cosphi;
		const REAL G = ER.G(), K = ER.K();
		const REAL c1 = 2.*G*(1.+1./3.*sinpsi) + 2.*K*sinpsi;
		const REAL c2 = (4.*G/3.-2.*K)*sinpsi;
		const REAL c3 = 2.*G*(1.-1./3.*sinpsi) - 2.*K*sinpsi;
		const REAL constA = 4.* G *(1.+ sinphi*sinpsi/3.) + 4.* K * sinphi*sinpsi;
		const REAL constB = 2.*G*(1+sinphi+sinpsi-1./3.*sinphi*sinpsi) + 4.*K*sinphi*sinpsi;
		REAL epsbar = epsbarp;
		REAL c, H;
		PlasticityFunction(epsbar, c, H);
		const REAL cos2H4 = 4.*cosphi2*H;
		const REAL denom = (constA-constB)*(constA+constB+8.*cosphi2*H);
		const REAL dGads1 = (-constB*(1+sinphi)+constA*(1+sinphi))/denom; // Derivative of DgammaA with respect to Sigma1tr
		const REAL dGads2 = (-cos2H4*(-1.+sinphi)-constB*(-1.+sinphi))/denom; // Derivative of DgammaA with respect to Sigma2tr 
		const REAL dGads3 = (cos2H4*(-1.+sinphi)+constA*(-1.+sinphi))/denom; // Derivative of DgammaA with respect to Sigma3tr
		const REAL dGbds1 = (constA*(1.+sinphi)-constB*(1.+sinphi))/denom; // Derivative of DgammaB with respect to Sigma1tr
		const REAL dGbds2 = (cos2H4*(-1.+sinphi)+constA*(-1.+sinphi))/denom; // Derivative of DgammaA with respect to Sigma2tr 
		const REAL dGbds3 = (-cos2H4*(-1.+sinphi)-constB*(-1.+sinphi))/denom; // Derivative of DgammaA with respect to Sigma3tr
				
		tang.Redim(3, 3);

		// First column
		tang(0,0) = 1.-c1*(dGads1+dGbds1);
		tang(1,0) = c2*dGads1+c3*dGbds1;
		tang(2,0) = c3*dGads1+c2*dGbds1;

		// Second column
		tang(0,1) = -c1*(dGads2+dGbds2);
		tang(1,1) = 1.+c2*dGads2+c3*dGbds2;
		tang(2,1) = c3*dGads2+c2*dGbds2;
		
		// Third column
		tang(0,2) = -c1*(dGads3+dGbds3);
		tang(1,2) = c2*dGads3+c3*dGbds3;
		tang(2,2) = 1.+c3*dGads3+c2*dGbds3;
	}	
	
	
};


#endif //TPZYCMOHRCOULOMBPV_H
