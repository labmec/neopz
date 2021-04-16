/**
 * @file pzmaterial.h
 * @brief Header file for abstract class TPZMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef PZVECL2HPP
#define PZVECL2HPP

#include "TPZMaterial.h"

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the TPZMaterial class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * TPZMaterial objects also need to implement the interface for post processing the results
 */
class  TPZVecL2 : public TPZMaterial
{
    
    
public:
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZVecL2(int id);
    
    /** @brief Default constructor */
    TPZVecL2();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZVecL2(const TPZVecL2 &mat);
    
    TPZVecL2 &operator=(const TPZVecL2 &mat);
    
    /** @brief Default destructor */
    virtual ~TPZVecL2();
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */

    /** @brief Returns the name of the material */
    virtual std::string Name()  override { return "TPZVecL2"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const  override
    {
        return fDim;
    }
    
    virtual void SetDimension(int dim){
        fDim = dim;
    }
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override
    {
        return fNState;
    }

    virtual void SetNStateVariables(int nstate)
    {
        fNState = nstate;
    }
    
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout) override;
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name) override;
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var) override;
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * left, TPZCompEl * ritgh) override;
    
protected:
    /** @deprecated Deprecated interface for Solution method which must use material data. */
    virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
    
    /** @brief Problem dimension */
    int fDim;
    
    /// Number of state variables
    int fNState = 1;
    
public:
    
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override
    {
        TPZFMatrix<STATE> ek(ef.Rows(),ef.Rows(),0.);
        Contribute(datavec, weight, ek, ef);
    }
    
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
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
	 * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override
    {
        
    }
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
	
    /** @} */


        /** @brief To create another material of the same type*/
        virtual TPZMaterial *NewMaterial() override;
protected:
    void Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                TPZFMatrix<REAL> &axes, TPZVec<STATE> &u_exact,
                TPZFMatrix<STATE> &curlU_exact, TPZVec<REAL> &val) override;
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    public:
  virtual int ClassId() const override;

    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context) override;
    
    /** @} */
	
};

#endif

