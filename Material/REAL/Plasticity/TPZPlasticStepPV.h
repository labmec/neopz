/**
 * @file
 */

#ifndef TPZPlasticStepPV_H
#define TPZPlasticStepPV_H


#include "TPZTensor.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZPlasticState.h"
#include "TPZPlasticIntegrMem.h"
#include "TPZPlasticStep.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include <set>
#include <ostream>

// Metodos para deixar o prog mais "encapsulado"
REAL NormVecOfMat(TPZFNMatrix <9> mat);
REAL InnerVecOfMat(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2);
TPZFMatrix<REAL> ProdT(TPZFMatrix<REAL> &m1,TPZFMatrix<REAL> &m2);
TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat);

template <class TVar>
class TPZFMatrix;

/*
 
 enum EElastoPlastic
 {
 EAuto = 0,
 EForceElastic = 1,
 EForcePlastic = 2
 };
 */

/**
 * @brief Classe que efetua avanco de um passo de plastificacao utilizando o metodo de Newton
 */
template <class YC_t, class ER_t>
class TPZPlasticStepPV : public TPZPlasticBase
{
public:
	
	/**
	 * @brief Constructor which Initialize the plastic material damage variable only
	 *
	 * @param[in] alpha damage variable
	 */
	
  TPZPlasticStepPV(REAL alpha=0.):fYC(), fER(), fResTol(1.e-12), fMaxNewton(30), fN()
	{ 
		fN.fAlpha = alpha; 
	}
	
	/**
	 * @brief Copy Constructor
	 *
	 * @param[in] source of copy
	 */
	TPZPlasticStepPV(const TPZPlasticStepPV & source)
	{
		fYC			= source.fYC;
		fER			= source.fER;
		fResTol		= source.fResTol;
		fMaxNewton = source.fMaxNewton;
		fN			= source.fN;
	}
	
	/**
	 * @brief Operator =
	 *
	 * @param[in] source of copy
	 */
	TPZPlasticStepPV & operator=(const TPZPlasticStepPV & source)
	{
		fYC			= source.fYC;
		fER			= source.fER;
		fResTol		= source.fResTol;
		fMaxNewton = source.fMaxNewton;
		fN			= source.fN;

		return *this;
	}
	
	/**
	 * @brief Name of the class ina string
	 */
	virtual const char * Name() const
	{
		return "TPZPlasticStepPV";	
	}
	
	virtual void Print(std::ostream & out) const
	{
		out << "\n" << this->Name();
		out << "\n YC_t:";
		//fYC.Print(out); FAZER O PRINT
		out << "\n ER_t:";
		fER.Print(out);
		out << "\nTPZPlasticStepPV Internal members:";
		out << "\n fResTol = " << fResTol;
		out << "\n fMaxNewton = " << fMaxNewton;
		out << "\n fN = "; // PlasticState
		fN.Print(out);
	}
		
	typedef YC_t fNYields;
	
	/**
	 * Imposes the specified strain tensor, evaluating the plastic integration if necessary.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal);
    
	/**
	 * Imposes the specified strain tensor and returns the correspondent stress state.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 */	
	virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma);
    
    
    
    //virtual void ApplySigmaComputeStrain(const TPZTensor<REAL> &sigma, TPZTensor<REAL> &epsTotal);
	
	/**
	 * Imposes the specified strain tensor and returns the corresp. stress state and tangent
	 * stiffness matrix.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 * @param[out] Dep Incremental constitutive relation
	 */
	virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep);
    
    /**
	 * Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
	 * This method should be used only for test purposes because it isn't fully robust. Some
	 * materials, specially those perfectly plastic and with softening, may fail when applying
	 * the Newton Method on ProcessLoad.
	 *
	 * @param[in] sigma stress tensor
	 * @param[out] epsTotal deformation tensor
	 */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal);
    
    virtual TPZPlasticState<REAL> GetState() const;
    /**
     * @brief Return the value of the yield functions for the given strain
     * @param[in] epsTotal strain tensor (total strain)
     * @param[out] phi vector of yield functions
	 */
	virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const;

    
    
    /**
	 * @brief Update the damage values
	 * @param[in] state Plastic state proposed
	 */
    virtual void SetState(const TPZPlasticState<REAL> &state);
    
    
    //void CopyFromFNMatrixToTensor(TPZFNMatrix<6> FNM,TPZTensor<STATE> &copy);
    void CopyFromTensorToFMatrix(TPZTensor<STATE> tensor,TPZFMatrix<STATE> &copy);
    

    //void CopyFromTensorToFNMatrix(TPZTensor<STATE> tensor,TPZFNMatrix<6> &copy);
    void CopyFromFMatrixToTensor(TPZFMatrix<STATE> FNM,TPZTensor<STATE> &copy);
    

    
    
    virtual void Read(TPZStream &buf);
    
    virtual void Write(TPZStream &buf) const;
    
	
	/**
	 * Does the TaylorCheck of the tangent matrix
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 * @param[out] Dep Incremental constitutive relation
	 */
	void TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev, TPZVec<REAL> &conv);

	
	REAL ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat);
	
	

public:
	
	void SetResidualTolerance(STATE tol)
	{
		fResTol = tol;
	}
    
    void ResetPlasticMem()
    {
        //fPlasticMem.Resize(0);
    }
	
	//virtual void Write(TPZStream &buf) const;
	
	//virtual void Read(TPZStream &buf);
public:
	
	/** @brief Object which represents the yield criterium */
	YC_t fYC;

	/** @brief Object representing the elastic response */
	ER_t fER;
	
protected:
	
	/** @brief Residual tolerance accepted in the plastic loop processes */
	REAL fResTol;
	
	/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
	int fMaxNewton;	// COLOCAR = 30 (sugestao do erick!)
    
	
	
public:
	
	/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */	
	TPZPlasticState<STATE> fN;
    
    int fYield;
	
	
};

#endif //TPZPlasticStepPV_H
