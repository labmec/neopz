#ifndef TPZCOHESIVEBCH
#define TPZCOHESIVEBCH
/**
 * @file
 * @brief Contains the TPZCohesiveBC class which implements a cohesive boundary condition
 */


#include <iostream>
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "TPZMatWithMem.h"
#include "pznlfluidstructureData.h"

const int TPZCohesiveBCID = 400;

/**
 * @ingroup material
 * @brief This class implements a cohesive bc, which has the stress dependent of the displacement.
 * Remembering that the std::pair memory holds the pair (DeltaT,SigmaT) at each time step.
 * @author Nathan Shauer in 01/04/2014
 */
class  TPZCohesiveBC : public TPZMatWithMem<TPZFMatrix<REAL> >
{
private:

	/// SigmaT (Maximum Traction) of the cohesive equation
	REAL fSigmaT;
	/// DeltaC (Critical Displacement) of the cohesive equation
	REAL fDeltaC;
	/// DeltaT (diplacement for maximum traction). The point (DeltaT,SigmaT)
	REAL fDeltaT;
	
	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	EState gState;	
	
public:
	
	
	/** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(int id);
	
	/** @brief Default constructor */
	TPZCohesiveBC();
	
	/** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(const TPZCohesiveBC &mat);
  
	/** @brief Default destructor */
	virtual ~TPZCohesiveBC();
	
	/** @brief Sets the Cohesive Data */
	void SetCohesiveData(const REAL &SigmaT, const REAL &DeltaC, const REAL &DeltaT);
	
	/** @brief Calculates Sigma for determined solution */
	void CalculateSigma(REAL &w,REAL &DeltaT, REAL &SigmaT, REAL &sigma, REAL &propageted) const;

	/** @brief Calculates DerivSigma for determined solution */	
	void CalculateCohesiveDerivative(REAL &w,REAL &DeltaT, REAL &SigmaT, REAL &deriv, REAL &propageted) const;	

	/** @brief Updates the cohesive curve acording to the calculated w of the time step */
	void UpdateCohesiveCurve(TPZMaterialData &data);

	void SetLastState(){ gState = ELastState; }
	void SetCurrentState(){ gState = ECurrentState; }
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override
	{
		datavec[0].SetAllRequirements(false);
		datavec[0].fNeedsSol = true;
	}
	
	/** @brief Returns the name of the material */
	virtual std::string Name()  override { return "TPZCohesiveBC"; }
	
	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() const override {
		return 0;
	}
  
  /** @brief Returns the number of state variables associated with the material */
  virtual int NStateVariables() const  override {return 2;}
  
	/** @brief Prints out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout) override;
	
public:

	
	/**
	 * @brief It computes a contribution to the residual vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the residual vector
	 * @since April 16, 2007
	 */
	//virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
	
 	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
	 * @param data stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
  
  /**
   * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
   * @param data [in] stores all input data
   * @param weight [in] is the weight of the integration rule
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition material
   * @since October 07, 2011
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
  
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
	 * @param datavec [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
	
	/** @brief Unique identifier for serialization purposes */
	int ClassId() const override {
		return TPZCohesiveBCID;
	}
	
};

#endif