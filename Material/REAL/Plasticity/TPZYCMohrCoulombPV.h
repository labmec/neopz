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
    
public:
    enum{NYield=3};
    
private:
	REAL fPhi;
	REAL fPsi;
	REAL fc;
	TPZElasticResponse fER;

protected:
	REAL fEpsPlasticBar;
	
public:
    
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
		
		enum MPlane {ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EApex };
		
		MPlane fWhichPlane;
		
		TPZManVector<REAL> fGamma;
	};
		
public:
	
	/**
	 * @brief empty constructor
	 */
	TPZYCMohrCoulombPV();

	/**
	 * @brief Constructor seting yc parameters
	 */
	TPZYCMohrCoulombPV(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER);
	
	/**
	 * @brief Copy Constructor
	 */
	TPZYCMohrCoulombPV(const TPZYCMohrCoulombPV &cp);
	
	/**
	 * @brief Sets up the data
	 */
	void SetUp(REAL Phi, REAL Psi, REAL c, TPZElasticResponse &ER){
		fPhi = Phi;
		fPsi = Psi;
		fc = c;
		fER = ER;
	}
	
	/**
	 * @brief Operator =
	 */
	TPZYCMohrCoulombPV & operator=(const TPZYCMohrCoulombPV &cp);
	
    
    void Read(TPZStream &buf);
    
    void Write(TPZStream &buf) const;
    

	/**
	 * @brief Sets epsbar
	 */
	void SetEpsBar(REAL &epsbar)
	{
		fEpsPlasticBar = epsbar;
	}
	 
	/**
	 * @brief Print Method
	 */
	void Print(std::ostream &out) const
	{
		out << "TPZYCMohrCoulombPV\n";
		out << "Still have to implement the print" << std::endl;
	}
	
	/**
	 * @brief Calculates the value c(epsp) and its derivative
	 */
	template <class T>
	void PlasticityFunction(const T epsp, T &c, T &H) const;
	
	/**
	 * @brief sigma = lambda Tr(E)I + 2 mu E
	 */
	template<class T>
	TPZVec<T> SigmaElastPV(const TPZVec<T> &deform) const;
	
	/**
	 * @brief Calcula o valor da funcao criteiro de plastificacao
	 */
	template<class T>
	T PhiPlane(const TPZVec<T> &sigma) const;
	
	/**
	 * @brief Implements the return map in the plane of the surface
	 */
	template<class T>
	bool ReturnMapPlane(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected, 
											TComputeSequence &memory, REAL &epsbarnew) const;

	/**
	 * @brief Computes dsigmapr/dsigmatr for the ReturnMapPlane
	 */
	void ComputePlaneTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;
	 
	/**
	 * @brief Implements the return map in the left edge of the surface
	 */
	template<class T>
	bool ReturnMapLeftEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
												 TComputeSequence &memory, REAL &epsbarnew) const;

	/**
	 * @brief Computes dsigmapr/dsigmatr for the ReturnMapLeftEdge
	 */
	void ComputeLeftEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;
	
	/**
	 * @brief Implements the return map in the right edge of the surface
	 */
	template<class T>
	bool ReturnMapRightEdge(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
													TComputeSequence &memory, REAL &epsbarnew) const;
	
	/**
	 * @brief Computes dsigmapr/dsigmatr for the ReturnMapRightEdge
	 */
	void ComputeRightEdgeTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;
	
	/**
	 * @brief Implements the return map in the apex
	 */
	template<class T>
	bool ReturnMapApex(const TPZVec<T> &sigma_trial, TPZVec<T> &sigma_projected,
													TComputeSequence &memory, REAL &epsbarnew) const;
	
	/**
	 * @brief Computes dsigmapr/dsigmatr for the ReturnMapApex
	 */
	void ComputeApexTangent(TPZMatrix<REAL> &tang, REAL &epsbarp) const;
	
	/**
	 * @brief Choses the correct projection and returns projected sigma and new epspbar
	 */
	void ProjectSigma(const TPZVec<STATE> &sigma_trial, STATE eprev, TPZVec<STATE> &sigma, STATE &eproj);
	
	/**
	 * @brief Choses the correct projection and returns projected sigma, new epspbar and tangent matrix
	 */
	void ProjectSigmaDep(const TPZVec<STATE> &sigmatrial, STATE kprev, TPZVec<STATE> &sigmaproj, STATE &kproj, TPZFMatrix<STATE> &tang);
    
	/**
	 * @brief Calculates the value of phi based on eps
	 */
  void Phi(TPZVec<STATE> sigvec,STATE alpha,TPZVec<STATE> &phi)const;

  STATE Phi() { return fPhi; }
  STATE Psi() { return fPsi; }
  STATE Cohesion() { return fc; }
  STATE E() { return fER.E(); }
  STATE Poisson() { return fER.Poisson(); }
	
};


#endif //TPZYCMOHRCOULOMBPV_H
