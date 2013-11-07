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

#include "pzlog.h"

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

class TPZPlasticBase
{
public:
	
  virtual	~TPZPlasticBase(){}; 
	/*virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) = 0;
	virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma) = 0;
	virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) = 0;
  virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) = 0;
  virtual void SetState(const TPZPlasticState<REAL> &state) = 0;
	virtual const TPZPlasticState<REAL> GetState() const = 0;
	virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const = 0;
	virtual int IntegrationSteps()const = 0;
	virtual void SetIntegrTol(REAL integrTol)=0;
	virtual const char * Name()const = 0;
	virtual void Print(std::ostream & out)const = 0;
	virtual void SetTensionSign(int sign) = 0;
	virtual void Write(TPZStream &buf) const = 0;
	virtual void Read(TPZStream &buf) = 0;*/
	
};

/**
 * @brief Classe que efetua avanco de um passo de plastificacao utilizando o metodo de Newton
 */
template <class YC_t, class ER_t>
class TPZPlasticStepPV: public TPZPlasticBase 
{
public:
	
	/**
	 * @brief Constructor which Initialize the plastic material damage variable only
	 *
	 * @param[in] alpha damage variable
	 */
	
  TPZPlasticStepPV(REAL alpha=0.):fYC(), fER(), fResTol(1.e-12), fMaxNewton(30)
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
	 * Imposes the specified strain tensor and returns the correspondent stress state.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 */	
	virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma);
	
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
	 * Does the TaylorCheck of the tangent matrix
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 * @param[out] Dep Incremental constitutive relation
	 */
	void TaylorCheck(TPZTensor<REAL> &EpsIni, TPZTensor<REAL> &deps, REAL kprev);
	
	REAL ComputeNFromTaylorCheck(REAL alpha1, REAL alpha2, TPZFMatrix<REAL> &error1Mat, TPZFMatrix<REAL> &error2Mat);
	
	/**
	 * Defines whether the tension/extension sign is desired by the external user.
	 *
	 * @param s [in] sign (1 or -1)
	 */
	//void SetTensionSign(int s);
	
	
protected:
	

	
public:
	/**
	 * @brief Return the value of the yield functions for the given strain
	 * @param[in] epsTotal strain tensor (total strain)
	 * @param[out] phi vector of yield functions
	 */
	//virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const;
	
		
	/**
	 * @brief Update the damage values
	 * @param[in] state Plastic state proposed
	 */
	//virtual void SetState(const TPZPlasticState<REAL> &state);

	
	/** @brief Retrieve the plastic state variables */	
	//virtual const TPZPlasticState<REAL> GetState()const;
	
	/** @brief Return the number of plastic steps in the last load step. Zero indicates elastic loading. */
	virtual int IntegrationSteps() const
	{
		return 1;
	}
	

protected:	
	
	/**
	 * @brief Updates the damage values - makes no interface sign checks
	 * @param[in] state Plastic state proposed
	 */
	//virtual void SetState_Internal(const TPZPlasticState<REAL> &state);




public:
	
	void SetResidualTolerance(STATE tol)
	{
		fResTol = tol;
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
	
protected:
	
	/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */	
	TPZPlasticState<REAL> fN;
	
	/**
	 * @brief Stores the plastic evolution in the last evaluated PlasticIntegration call,
	 * It includes the N-1 data, the elastic step until yield when it exists,
	 * the plastic substeppings and the N step.
	 */
	//TPZStack< TPZPlasticIntegrMem<REAL, YC_t::NYield> > fPlasticMem;

	
	//ofstream fOutfile(string &str);
	
	
};

#endif //TPZPlasticStepPV_H
