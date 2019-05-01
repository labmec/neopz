/**
 * @file pzmaterial.h
 * @brief Header file for uncoupled TPZUncoupledMultiPhysics.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef PZUNCOUPLEDMULTIPHYSICSHPP
#define PZUNCOUPLEDMULTIPHYSICSHPP

#include "TPZMaterial.h"

#include <iostream>
#include <string>

class TPZBndCond;
class TPZMaterialData;
class TPZIntPoints;

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the TPZUncoupledMultiPhysics class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * TPZUncoupledMultiPhysics objects also need to implement the interface for post processing the results
 */
class  TPZUncoupledMultiPhysics : public TPZMaterial
{
private:
    int fId;
    
protected:

    TPZVec<TPZMaterial *> fReferredMaterials;

    
public:
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZUncoupledMultiPhysics(int id);
    
    /** @brief Default constructor */
    TPZUncoupledMultiPhysics();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZUncoupledMultiPhysics(const TPZUncoupledMultiPhysics &mat);
    /** @brief Default destructor */
    virtual ~TPZUncoupledMultiPhysics();
    
    /**
     * @brief Set the materials to which this material refers
     */
    void SetReferredMaterials(TPZVec<TPZMaterial *> &refer)
    {
        fReferredMaterials = refer;
    }
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** 
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data) override ;
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override ;
    
    
    /** @brief Returns the name of the material */
    virtual std::string Name()  override { return "TPZUncoupledMultiPhysics"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override 
    {
        return fReferredMaterials[0]->Dimension();
    }
    
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override
    {
        return 0;
    }
    
    
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout) override ;
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name) override ;
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var) override ;
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override ;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override ;
    
public:
    
    /** @brief Creates an object TPZBndCond derived of TPZUncoupledMultiPhysics*/
    virtual TPZBndCond *CreateBC(TPZMaterial *reference, int id, int typ, TPZFMatrix<STATE> &val1,
                                 TPZFMatrix<STATE> &val2) override ;
    
    /** @name Contribute methods
	 * @{
	 */
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override ;
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override ;
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) = 0;
	
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
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override ;
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override ;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override ;
	
    /** @} */
	
    
    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const override ;
	
	/** @brief Gets the order of the integration rule necessary to integrate an element multiphysic */
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const override ;
	
    /** @brief To create another material of the same type*/
    virtual TPZUncoupledMultiPhysics * NewMaterial() override ;
    
    /** @brief Reads data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data) override ;
    
    /** @brief Creates a copy of the material object and put it in the vector which is passed on */
    virtual void Clone(std::map<int, TPZMaterial * > &matvec) override ;
    
    
    /** @brief Unique identifier for serialization purposes */
    public:
virtual int ClassId() const override ;

    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const override ;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context) override ;
    
    
};


#endif

