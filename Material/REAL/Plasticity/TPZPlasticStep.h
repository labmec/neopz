/**
 * @file
 */

#ifndef TPZPLASTICSTEP_H
#define TPZPLASTICSTEP_H

#include "TPZTensor.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "TPZPlasticState.h"
#include "TPZPlasticIntegrMem.h"
#include "TPZElasticResponse.h"

#include "pzlog.h"
#include "TPZPlasticBase.h"

#include <set>
#include <ostream>

enum EElastoPlastic
{
	EAuto = 0,
	EForceElastic = 1,
	EForcePlastic = 2
};

/**
 * @brief Classe que efetua avanco de um passo de plastificacao utilizando o metodo de Newton
 */
template <class YC_t, class TF_t, class ER_t>
class TPZPlasticStep: public TPZPlasticBase 
{
public:
	
virtual int ClassId() const override;
	/**
	 * Initialize the plastic material damage variable only
	 *
	 * @param[in] alpha damage variable
	 */
	
    TPZPlasticStep(REAL alpha=0.):fYC(), fTFA(), fER(), fResTol(1.e-12),
	fIntegrTol(1.e-4), fMaxNewton(30), fMinLambda(1.),
	fMinStepSize(1.e-3), fN(), fPlasticMem(),
	fMaterialTensionSign(1), fInterfaceTensionSign(1)
	{ fN.m_hardening = alpha; }
	
    TPZPlasticStep(const TPZPlasticStep & source)
	{
		fYC			= source.fYC;
		fTFA		= source.fTFA;
		fER			= source.fER;
		fResTol		= source.fResTol;
		fIntegrTol  = source.fIntegrTol;
		fMaxNewton  = source.fMaxNewton;	
		fMinLambda  = source.fMinLambda;
		fMinStepSize= source.fMinStepSize;
		fN			= source.fN;
		fPlasticMem = source.fPlasticMem;
		fMaterialTensionSign  = source.fMaterialTensionSign;
		fInterfaceTensionSign = source.fInterfaceTensionSign;
	}
	
	TPZPlasticStep & operator=(const TPZPlasticStep & source)
	{
		fYC			= source.fYC;
		fTFA		= source.fTFA;
		fER			= source.fER;
		fResTol		= source.fResTol;
		fIntegrTol  = source.fIntegrTol;
		fMaxNewton  = source.fMaxNewton;	
		fMinLambda  = source.fMinLambda;
		fMinStepSize= source.fMinStepSize;
		fN			= source.fN;
		fPlasticMem = source.fPlasticMem;
		fMaterialTensionSign  = source.fMaterialTensionSign;
		fInterfaceTensionSign = source.fInterfaceTensionSign;
		return *this;
	}
	
	virtual const char *Name() const override
	{
		return "TPZPlasticStep generic class";	
	}
	
	virtual void Print(std::ostream & out) const override
	{
		out << "\n" << this->Name();
		out << "\n YC_t:";
		fYC.Print(out);
		out << "\n TF_t:";
		fTFA.Print(out);	
		out << "\n ER_t:";
		fER.Print(out);
		out << "\nTPZPlasticStep Internal members:";
		out << "\n fResTol = " << fResTol;
		out << "\n fMaxNewton = " << fMaxNewton;
		out << "\n fMinLambda = " << fMinLambda;
		out << "\n fPlasticMem.NElements() = " << fPlasticMem.NElements();
		out << "\n fN = ";
		fN.Print(out);
	}
	
	/**
	 * Imposes the specified strain tensor, evaluating the plastic integration if necessary.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 */
	virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override;
    
    
    void SetOutFile(std::string outfile);
	
	typedef YC_t fNYields;
	
	/**
	 * Imposes the specified strain tensor and returns the correspondent stress state.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 */	
	virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent = NULL) override;
	
	/**
	 * Imposes the specified strain tensor and returns the corresp. stress state and tangent
	 * stiffness matrix.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 * @param[out] Dep Incremental constitutive relation
	 */
	virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override;
	
	/**
	 * Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
	 * This method should be used only for test purposes because it isn't fully robust. Some
	 * materials, specially those perfectly plastic and with softening, may fail when applying
	 * the Newton Method on ProcessLoad.
	 *
	 * @param[in] sigma stress tensor
	 * @param[out] epsTotal deformation tensor
	 */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override;
	
    /**
     * @brief Reset the plastic memory
     */
    void ResetPlasticMem()
    {
        fPlasticMem.Resize(0);
    }
    
	/**
	 * Defines whether the tension/extension sign is desired by the external user.
	 *
	 * @param s [in] sign (1 or -1)
	 */
	void SetTensionSign(int s);
	
	/** @brief Indicates whether or not to correct Stress/Strain sign */
	int SignCorrection()const;
        
        TPZPlasticCriterion& GetYC() override{
            return fYC;
        }

	
protected:
	
	/**
	 * @brief Imposes the specified strain tensor, evaluating the plastic integration if necessary. \n
	 * Internal Method
	 * @param[in] epsTotal Imposed total strain tensor
	 */
    void ApplyStrain_Internal(const TPZTensor<REAL> &epsTotal);
	
	/**
	 * @brief Imposes the specified strain tensor and returns the corresp. stress state and tangent
	 * stiffness matrix.
	 * Internal Method
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 * @param[out] Dep Incremental constitutive relation
	 */
    void ApplyStrainComputeDep_Internal(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep);
	
	/**
	 * @brief Imposes the specified strain tensor and returns the correspondent stress state.
	 * Internal Method
	 *
	 *	@note The method 	ApplyStrainComputeSigma calls Apply Strain followed by an evaluation
	 *	of the Elastic Response. This method is of very low efficiency because it evaluates the
	 *	whole ApplyStrain, evaluating sigma internally and discarding it. It evaluates
	 *	ComputePlasticVars in order to retrieve the value of sigma. This method should be used
	 *	for test purposes only or pieces of code where efficiency is not important.
	 *	The methods ApplyStrain and ApplyStrainComputeDep shoud be the most useful
	 *	implementations.
	 *
	 * @param[in] epsTotal Imposed total strain tensor
	 * @param[out] sigma Resultant stress
	 *
	 */
	void ApplyStrainComputeSigma_Internal(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
	{
		ApplyStrain_Internal(epsTotal);
		
		REAL A;
		
		int n = fPlasticMem.NElements();
		
		ComputePlasticVars(fPlasticMem[n-1].m_elastoplastic_state, sigma, A);
	}
	
	/**
	 * @brief Attempts to compute an epsTotal value in order to reach an imposed stress state sigma.
	 * This methid should be used only for test purposes because it isn't fully robust. Some
	 * materials, specially those perfectly plastic and with softening, may fail when applying
	 * the Newton Method on ProcessLoad.
	 * Internal Method
	 *
	 * @param[in] sigma stress tensor
	 * @param[out] epsTotal deformation tensor
	 */
    void ApplyLoad_Internal(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal);
	
	
public:
    /**
     * @brief Return the value of the yield functions for the given strain
     * @param[in] epsTotal strain tensor (total strain)
     * @param[out] phi vector of yield functions
	 */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override;
	
	
public:
	
    /**
     * @brief Return the value of the yield functions for the given strain
	 * Internal Method
     * @param[in] epsTotal strain tensor (total strain)
     * @param[out] phi vector of yield functions
     */
    virtual void Phi_Internal(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const;
	
	/**
	 * @brief Verifies if the proposed epsTotalNp1 is still in the elastic range
	 * @param[in] state Plastic state proposed
	 */
    bool IsStrainElastic(const TPZPlasticState<REAL> &state)const;
	
    /// modify the elastic response. Needs to be reimplemented for each instantiation
    virtual void SetElasticResponse(TPZElasticResponse &ER) override;

    virtual TPZElasticResponse GetElasticResponse() const override;
    /**
	 * @brief Update the damage values
	 * @param[in] state Plastic state proposed
	 */
    virtual void SetState(const TPZPlasticState<REAL> &state) override;
    
    /**
	 * @brief Overwrite the current total strain only
	 * @param[in] epsTotal Strain
	 */
    virtual void SetUp(const TPZTensor<REAL> & epsTotal);
	
    /** @brief Retrieve the plastic state variables */	
    virtual TPZPlasticState<REAL> GetState()const override;
	
	/** @brief Return the number of plastic steps in the last load step. Zero indicates elastic loading. */
	virtual int IntegrationSteps()const override;
	
	/** @brief Sets the tolerance allowed in the pde integration */
	virtual void SetIntegrTol(REAL integrTol)
	{
		fIntegrTol = integrTol;
	}
	
	/** @brief Sets the minimum loading substep size in the plastic integration */
	virtual void SetMinStepSize(REAL minStepSize)
	{
		fMinStepSize = minStepSize;
	}
	
	/**
	 * @brief Similar to IntegrationSteps, it returns the plastic parcel of the last loading.
	 * 1.0 means that the whole step was plastic, 0.0 means the whole step was elastic.
	 * @param plastifLen [out] returns the plastified lenght relative to the whole load
	 * length in which each yield function took part. 1.0 means that that specific
	 * yield function was ruling the plastification process in the whole loading step and
	 * 0.0 that that yield function did not take part in the plastification process.
	 * Each result of plastifLen may vary in the range 0.0 ~ 1.0 and the summation of
	 * all elements may vary in the range 0.0 ~ number_of_yield_functions, because
	 * more than one yield function may plastify simultaneously in the case of a
	 * subdifferentiation in non-smooth yield function grouping.
	 */
	virtual REAL IntegrationOverview(TPZVec<REAL> & plastifLen);
	
protected:	
	
    /**
	 * @brief Updates the damage values - makes no interface sign checks
	 * @param[in] state Plastic state proposed
	 */
    virtual void SetState_Internal(const TPZPlasticState<REAL> &state);
	
public:
    /** @brief Retrieves the plastic state variables - makes no interface sign checks */	
    virtual TPZPlasticState<REAL> GetState_Internal()const;
	
protected:

    /**
	 * @brief Imposes the specified strain tensor and performs plastic integration when necessary.
	 * This function creates a new plastic integration history epserimenting the proposed
	 * epsTotal. It does not update the current plastic state.
	 * @param ep Type indicating if it is elastic or plastic 
	 * @param epsTotal Imposed total strain tensor
	 */
    virtual void ProcessStrain(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep = EAuto);

    /**
     * @brief Imposes the specified strain tensor and performs plastic integration when necessary.
     * This function DO NOT calls PlasticIntegrate
     * @param epsTotal Imposed total strain tensor
	 * @param ep Type indicating if it is elastic or plastic 
     */
    virtual void ProcessStrainNoSubIncrement(const TPZTensor<REAL> &epsTotal, const EElastoPlastic ep = EAuto);
	
    /**
	 * @brief Imposes the specified stress tensor and performs plastic integration when necessary.
	 * This function evaluates a newton's method on ProcessStrain until the stress state matches.
	 * It does not update the current plastic state.
	 * @param sigma Imposed total strain tensor
	 * @param ep Type indicating if it is elastic or plastic 
	 */
    virtual void ProcessLoad(const TPZTensor<REAL> &sigma, const EElastoPlastic ep = EAuto);

	/**
	 * @brief Evaluates the constitutive matrix (DSigma/DEpsT) based on the data from the plastic
	 * integration history without modifying it.
	 * @param[out] sigma resultant stress tensor
	 * @param[out] Dep Incremental constitutive relation
	 */
    virtual void ComputeDep(TPZTensor<REAL> & sigma, TPZFMatrix<REAL> &Dep);
//    virtual void ComputeDep2(TPZTensor<REAL> & sigma, TPZFMatrix<REAL> &Dep);
	
	/**
	 *  @brief Evaluates the sigma stress tensor and the thermoForceA for use in several
	 *  pieces of this code.
	 *  The output stress is elastic, admitting that no plastification occurs.
	 * @param[in] state_T Plastic variables state
	 * @param[out] sigma_T resultant predicted elastic strain tensor
	 * @param[out] A_T resultant thermoForce A
	 */
	template<class T>
    void ComputePlasticVars(const TPZPlasticState<T> & state_T,
							TPZTensor<T> & sigma_T,
							T& A_T)const;
	
	/**
	 * @brief Finds the strain point in the linear path from epsTotalN and
	 * towards epsTotalNp1 that matches the yield criterium at at least one
	 * plastic surface.
	 * Returns as its value the deltaEpsT multiplier
	 * @param[in] epsTotalNp1 Proposed total strain
	 * @param[out] stateAtYield Plastic state at yield of at least one plastic surface.
	 */
	REAL FindPointAtYield(const TPZTensor<REAL>& epsTotalNp1,
						  TPZPlasticState<REAL> & stateAtYield)const ;
	
	/**
	 Evaluates the residual vector for the plasticity problem
	 @param[in] N_T1 Plastic state variables at the time N
	 @param[in] Np1_T2 Plastic state variables at the time N+1
	 @param[in] delGamma_T2 plastic multipliers, one for each yield function
	 @param[out] res_T2 computed residual
	 @param[out] normEpsPErr norm of the estimation of the relative plastic strain error
	 @param[in] silent indicates whether or not to generate Log4cxx logs
     * The template parameter will be either REAL or a FAD parameter
	 */
	
    template <class T1, class T2>
    void PlasticResidual(const TPZPlasticState<T1> &N_T1,
						 TPZPlasticState<T2> &Np1_T2,
						 const TPZVec<T2> &delGamma_T2,
						       TPZVec<T2> &res_T2,
                               REAL &normEpsPErr,
                               int silent = 0)const; 
    

    /**
	 * @brief Evaluates the residual vector for the plasticity problem RK
     * @param [in] N_T1 Plastic state variables at the time N
     * @param [in] Np1_T2 Plastic state variables at the time N+1
     * @param [in] delGamma_T2 plastic multipliers, one for each yield function
     * @param [out] res_T2 computed residual
     * @param [out] normEpsPErr norm of the estimation of the relative plastic strain error
     * @param [in] silent indicates whether or not to generate Log4cxx logs
     * The template parameter will be either REAL or a FAD parameter
     */
    template <class T1, class T2>
    void PlasticResidualRK(
                                                             const TPZPlasticState<T1> &N_T1,
                                                             TPZPlasticState<T2> &Np1_T2,
                                                             const TPZVec<T2> &delGamma_T2,
                                                             TPZVec<T2> &res_T2,
                                                             REAL &normEpsPErr,
                                                             int silent=0)const;
		
	/**
	 * @brief Updates the N+1 plastic state variables based on the solution of a Newton's scheme.
	 * A very simple line search is performed in order to attempt to guarantee residual drop.
	 * @param[in] N_T1 Plastic state variables at time N
	 * @param[out] Np1_T2 Proposed plastic state variables at time Np1
	 * @param[in] delGamma_T2 plastic multipliers, one for each yield function
	 * @param[in] res_T2 computed residual
	 * @param[in] Sol_TVECTOR Opposite (in sign) direction of residual drop - results from Newton's scheme
	 * @param[in] validEqs Indicates the equations not to use in residual norm computations
	 * @param[in] updateVars when 1, the Np1_T2 variables return updated, otherwise untouched. It is important to set it to 0 when the sol vector contains derivatives
	 */
	
    template<class T1, class T2, class TVECTOR>
    REAL UpdatePlasticVars(const TPZPlasticState<T1> &N_T1,
						   TPZPlasticState<T2> &Np1_T2,
						   TPZVec<T2> &delGamma_T2,
						   TPZVec<T2> &res_T2,
						   TVECTOR & Sol_TVECTOR,
						   TPZVec<int> & validEqs,
						   int updateVars = 1)const;   
	
	/**
	 * @brief This method copies the REAL variables into FAD variables and
	 * intializes the derivatives
	 * The nVars variable is usually the number of plastic independent 
	 * variables. When more independent variables become necessary,
	 * set this value accordingly.
	 * @param[in] state Plastic state variables
	 * @param[in] delGamma Plastic multiplier
	 * @param[out] state_T Plastic state variables
	 * @param[out] delGamma_T Plastic multiplier
	 * @param[in] nVars number of independent variables
	 */
	template<class T>
	void InitializePlasticFAD(const TPZPlasticState<REAL> & state, const TPZVec<REAL> &delGamma,
							  TPZPlasticState<T> &state_T,         TPZVec<T> &delGamma_T,
							  const int nVars = 7 + YC_t::NYield)const;
	
	/**
	 * @brief Initializes the fValidEqs booleans indication whether to consider the plastic surfaces
	 * @param[in] res_T Valid equations are only those leading to negative phi's
	 * @param[out] validEqs Validated equations
	 */
	template <class T>
	int InitializeValidEqs(TPZVec<T> &res_T, TPZVec<int> & validEqs);
	
	/**
	 * @brief After a complete plasticLoop, plsatic surface equations related to negative plastic multipliers are discarded.
	 * @param delGamma_T plastic multiplier
	 * @param res_T Valid equations are only those leading to negative phi's
	 * @param validEqs Validated equations
	 */
    template <class T>
    int RemoveInvalidEqs(TPZVec<T> & delGamma_T, TPZVec<T> & res_T, TPZVec<int> &validEqs);
	
    /**
     * @brief Proposes an update to the plastic variables and estimates the relative error
	 * comitted in this update. Neither internal variable are used nor changed.
	 * In the Np1 variables, EpsT is imposed [in] and the Alpha and EpsP are evaluated.
	 * It returns 1 if suceeded of 0 if tolerance not achieved.
     * @param N [in] Plastic state variables at time N
	 * @param Np1 [in/out] Plastic state variables at time N+1
     * @param delGamma [in/out] plastic multipliers
     */
    void InitialGuess(
                      const TPZPlasticState<REAL> &N,
                      TPZPlasticState<REAL> &Np1,
                      TPZVec<REAL> &delGamma,
                      TPZVec<int> &validEqs
                      );
    /**
     * @brief Proposes an update to the plastic variables and estimates the relative error
	 * comitted in this update. Neither internal variable are used nor changed.
	 * In the Np1 variables, EpsT is imposed [in] and the Alpha and EpsP are evaluated.
	 * It returns 1 if suceeded of 0 if tolerance not achieved.
     * @param N [in] Plastic state variables at time N
	 * @param Np1 [in/out] Plastic state variables at time N+1
     * @param delGamma [in/out] plastic multipliers
	 * @param normEpsPErr [out] estimated norm of the relative error deltaEpsP/deltaEpsTotal
     * @param lambda [out] Line search multiplier
	 * @param validEqs [out] Valid set of plastic flow eqs
	 */
    int PlasticLoop(const TPZPlasticState<REAL> &N,
					TPZPlasticState<REAL> &Np1,
					TPZVec<REAL> &delGamma,
					REAL &normEpsPErr,
					REAL &lambda,
					TPZVec<int> & validEqs);
	
	/**
     * @brief Proposes an update to the plastic variables respecting an integration
	 * with error control.
	 * Neither internal variable are used nor changed.
	 * In the Np1 variables, EpsT is imposed \[input\] and the Alpha and EpsP are evaluated.
     * @param N [in] Plastic state variables at time N
	 * @param Np1 [in/out] Plastic state variables at time N+1
	 * @param TolEpsPErr [in] relative error deltaEpsP/deltaEpsTotal norm to respect
	 */
    int PlasticIntegrate(
						 const TPZPlasticState<REAL> &N,
						 TPZPlasticState<REAL> &Np1,
						 const REAL &TolEpsPErr);
	
	
	/**
	 * @brief Extracts the tangent matrix and residual vector from the FAD variables 
	 * and according to the preconditioning and fValidEqs booleans.
	 * @param[in] epsRes_FAD Residual of plastic eqs in FAD fashion
	 * @param[out] ResVal Residual of plastic eqs in vector fashion, preconditioned if requested
	 * @param[out] resnorm L2 norm of unpreconditioned ResVal
	 * @param[out] tangent Residual Jacobian in matrix fashion, preconditioned if requested
	 * @param[in] validEqs boolean indicating the valid plastic flow eqs
	 * @param[in] precond Boolean indicating whether or not to precondition the matrix/residual vector
	 * @param[in] resetInvalidEqs Resets the plastic surfaces assigned to be discarded
	 */ 
    template<class T1, class T_VECTOR, class T_MATRIX>//T1:input residual fad type (FAD FAD or FAD), T_MATRIX: output matrix &vector of type (FAD or REAL, respectvly)
    void ExtractTangent(
						const TPZVec<T1> & epsRes_FAD,
						T_VECTOR & ResVal, //TPZFMatrix<REAL> for the T1=fad<real> type
						REAL & resnorm, // REAL 
						T_MATRIX & tangent, // TPZFMatrix<REAL> for T1=fad<real> type
                        TPZVec<int> & validEqs,
                        TPZVec<REAL> &pivots,
						const int precond = 1,
						const int resetInvalidEqs = 1);
	
	/**
	 * @brief Stores the whole content of a plastic integration step in order to allow its 
	 * further reevaluation.
	 * @param state [in] plastic state variables
	 * @param k [in] deltaEpsTotal multiplier
	 * @param lambda [in] line search multiplier in the newton's method within plasticloop
	 * @param delGamma [in] Plastic flow multiplier
	 * @param validEqs [in] valid equations involved in one step of the integration process
	 * @param forceYield force id
	 */
	template <int N>
	void PushPlasticMem(const TPZPlasticState<REAL> & state,
						const REAL & k,
						const REAL & lambda,
						const TPZManVector<REAL,N> & delGamma,
						const TPZVec<int> & validEqs,
						const int forceYield);
	
public:
	void SetMaterialElasticOrPlastic(int choice=1/*plastic*/);
    
    void SetResidualTolerance(STATE tol)
    {
        fResTol = tol;
    }
	
    void Write(TPZStream& buf, int withclassid) const override;
    void Read(TPZStream& buf, void* context) override;
public:
	
    /** @brief Object which represents the yield criterium */
    YC_t fYC;
    /** @brief Object representing the thermodynamical force */
    TF_t fTFA;
    /** @brief Object representing the elastic response */
    ER_t fER;
	
	
protected:
	
	/** @brief Residual tolerance accepted in the plastic loop processes */
	REAL fResTol;
	
public:
	/** @brief Tolerance desired in the Plastic integration processes */
	REAL fIntegrTol;
	
	/** @brief Maximum number of Newton interations allowed in the nonlinear solvers */
	int fMaxNewton;	
	
	/**
	 * @brief Minimum multiplicaton in the Plastic Loop line search. 1.0 avoids line searching;
	 * Very Small values do not allow the code to skip out of local minima.
	 */
	REAL fMinLambda;
	
	/**
	 * @brief Minimum fraction of the full load substep proposed accepted in the
	 * plastic integration. Too low values may lead to extremely slow integration
	 * in some cases while values as high as 1.0 may prevent the plastic integrator
	 * error control from adjusting the advisable plastic substeps. Values between
	 * 1.e-3 and 1.e-2 are advisable.
	 */
	REAL fMinStepSize;
	
protected:
	
	/** @brief Plastic State Variables (EpsT, EpsP, Alpha) at the current time step */	
    TPZPlasticState<REAL> fN;
	
	/**
	 * @brief Stores the plastic evolution in the last evaluated PlasticIntegration call,
	 * It includes the N-1 data, the elastic step until yield when it exists,
	 * the plastic substeppings and the N step.
	 */
	TPZStack< TPZPlasticIntegrMem<REAL, YC_t::NYield> > fPlasticMem;
	
	/** @brief The tension sign in the convention used in the implementation of the material */
	int fMaterialTensionSign;
	
	/** @brief The tension sign in the convention defined by the external user */
	int fInterfaceTensionSign;
	
    
    
    //ofstream fOutfile(std::string &str);
	
public:
	
	//////////////////CheckConv related methods/////////////////////
	
    /** @brief Number of types of residuals */
    int NumCases() 
    {
		return 1;
    }
	
    static TPZTensor<REAL> gRefDeform;
	
    /** @brief LoadState will keep a given state as static variable of the class */
    void LoadState(TPZFMatrix<REAL> &state)
    {
		int i;
		for(i=0; i<6; i++) gRefDeform[i] = state(i,0);
    }
	
    void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase)
    {
#ifdef PZ_LOG
		TPZLogger logger("plasticity.plasticstep");
#endif
		TPZTensor<REAL> sigma;
		switch(icase)
		{
			case 0:
				TPZTensor<REAL> sig;
				tangent.Redim(6,6);
			{
				ProcessStrain(gRefDeform);
				ComputeDep(sig, tangent);
			}
				//Sigma(gRefDeform,sig,tangent);
#ifdef PZ_LOG
				std::stringstream sout;
				sout << "matriz tangent for checkconv " << tangent;
				LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
				break;
		}
    }
	
    void Residual(TPZFMatrix<REAL> &res,int icase)
    {
#ifdef PZ_LOG
		TPZLogger logger("plasticity.plasticstep");
#endif
		TPZTensor<REAL> sigma;
		switch(icase)
		{
			case 0:
				TPZTensor<REAL> sig;
				TPZFMatrix<REAL> tangent(6,6);
				res.Redim(6,1);
			{
				ProcessStrain(gRefDeform);
				ComputeDep(sig, tangent);
			}
				//Sigma(gRefDeform,sig,tangent);
#ifdef PZ_LOG
				std::stringstream sout;
				sout << "sigma for residual " << sig;
				LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
				int i;
				for(i=0; i<6; i++)
				{
					res(i,0) = sig[i];
				}
				break;
		}
    }
	
	//////////////////CheckConv related methods/////////////////////

};

template <class YC_t, class TF_t, class ER_t>
int TPZPlasticStep<YC_t, TF_t, ER_t>::ClassId() const{
    return Hash("TPZPlasticStep") ^ TPZPlasticBase::ClassId() << 1 ^ YC_t().ClassId() << 2 ^  TF_t().ClassId() << 3 ^ ER_t().ClassId() << 4;
}

#endif //TPZPLASTICSTEP_H
