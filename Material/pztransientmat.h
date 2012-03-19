/**
 * \file
 * @brief Contains the TPZTransientMaterial class which implements an implicit Euler time integrator.
 */

//$Id: pztransientmat.h,v 1.7 2009-05-06 20:13:37 fortiago Exp $


#ifndef TRANSIENTMATH
#define TRANSIENTMATH

#include "pzmaterial.h"
// #include "pzvec.h"
// #include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief Implements an implicit Euler time integrator. \ref material Material 
 */
/**
 * Weak statement is supposed to be Integral[(un+1 - un)/deltaT * v, Omega] + Bilinear Form = Linear Form
 * This class implements only Integral[(un+1 - un)/deltaT * v, Omega]. Bilinear and linear form must be implemented in base class TBASEMAT.
 */
template<class TBASEMAT>
class TPZTransientMaterial : public TBASEMAT {
	
public:
	
	/** @brief Class constructor */
	TPZTransientMaterial(int nummat, int dim, REAL TimeStep);
	
	/** @brief Default destructor */
	~TPZTransientMaterial();
	
	/** @brief Copy constructor */
	TPZTransientMaterial(const TPZTransientMaterial &cp);
	
	/** @brief Sets integral scheme as an explicit Euler */
	void SetExplicit();
	
	/** √Sets integral scheme as an implicit Euler */
	void SetImplicit();
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<REAL> &ek,
                            TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<REAL> &ek,
                                     TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                       REAL weight,
                                       TPZFMatrix<REAL> &ek,
                                       TPZFMatrix<REAL> &ef,
                                       TPZBndCond &bc);
	
	/** @brief Set material to compute only Integral[- un/deltaT * v, Omega] */
	void SetLastState();
	
	/** @brief Set material to compute Integral[un+1/deltaT * v, Omega] + Bilinear Form = Linear Form  */
	void SetCurrentState();
	
	/** @brief Set material to compute ek = Integral[phi_i phi_j, Omega]/deltaT */
	void SetMassMatrix();
	
	/** @brief Set material to compute ef = Linear Form - Bilinear Form(u) = F -ku */ 
	void SetFluxOnly();
	
	/** @brief Define time step DeltaT */
	void SetTimeStep(REAL TimeStep);
	
	/** @brief Returns time step value. */
	REAL TimeStep();
	
	/** @brief Indicates if the material requires the solution to compute Contribute */
	/** 
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsSolutionToContribute(){
		return true;
	}
	
	/** @brief Indicates if the material requires the global coordinate X to compute Contribute */
	/**
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsXCoord(){
		return (this->fStep != ELast);
	}
	
protected:
	
	enum ETemporalScheme{EImplicit = 1, EExplicit = 2};
	
	ETemporalScheme fTemporalIntegrator;
	
	enum STEPS{ENone = -1, ELast = 0, ECurrent = 1, EMassMatrix = 2, EFluxOnly = 3};
	
	STEPS fStep;
	
	REAL fTimeStep;
	
	virtual void ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix<REAL> &phi, REAL weight, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeTangent(TPZVec<REAL> &sol, TPZFMatrix<REAL> &phi, REAL weight, TPZFMatrix<REAL> &ek);
};

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetLastState(){
	this->fStep = ELast;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetCurrentState(){
	this->fStep = ECurrent;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetMassMatrix(){
	this->fStep = EMassMatrix;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetFluxOnly(){
	this->fStep = EFluxOnly;
}

template<class TBASEMAT>
inline void TPZTransientMaterial< TBASEMAT >::SetTimeStep(REAL TimeStep){
	this->fTimeStep = TimeStep;
}

template<class TBASEMAT>
inline REAL TPZTransientMaterial< TBASEMAT >::TimeStep(){
	return this->fTimeStep;
}

#endif
