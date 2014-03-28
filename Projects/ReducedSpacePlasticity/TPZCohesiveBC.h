#ifndef TPZCOHESIVEBCH
#define TPZCOHESIVEBCH
/**
 * @file
 * @brief Contains the TPZCohesiveBC class which implements a cohesive boundary condition
 */


#include <iostream>
#include "pzmaterial.h"
#include "pzbndcond.h"

const int TPZCohesiveBCID = 400;

/**
 * @ingroup material
 * @brief This class implements a cohesive bc, which has the stress dependent of the displacement
 *        
 */
class  TPZCohesiveBC : public TPZMaterial
{
private:

	/// SigmaT (Traction) of the cohesive equation Sigma = SigmaT(1 - Delta/DeltaC) whete Delta is the displacement opening 
	REAL fSigmaT;
	/// DeltaC (Critical) of the cohesive equation Sigma = SigmaT(1 - Delta/DeltaC) whete SigmaT is the maximum traction 
	REAL fDeltaC;
	
public:
	
	/** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(int id);
	
	/** @brief Default constructor */
	TPZCohesiveBC();
	
	/** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
	TPZCohesiveBC(const TPZMaterial &mat);
	/** @brief Default destructor */
	virtual ~TPZCohesiveBC();
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
	{
		datavec[0].SetAllRequirements(false);
	}
	
	/** @brief Returns the name of the material */
	virtual std::string Name() { return "TPZCohesiveBC"; }
	
	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() {
		return 1;
	}
	

	/** @brief Returns the number of state variables associated with the material */
	virtual int NStateVariables() {
		return 1;
	}
		
	/** @brief Prints out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout);
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
	virtual int NSolutionVariables(int var);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
public:

	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
	 * to multiphysics simulation.
	 * @param datavec [in]  stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since October 18, 2011
	 */
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to the residual vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the residual vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const
	{
		return TPZCohesiveBCID;
	}
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	
};

#endif