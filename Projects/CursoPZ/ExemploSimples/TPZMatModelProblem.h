/**
 * @file TPZMatModelProblem.h
 * @brief Header file for class TPZMatModelProblem.\n
 * It implements the weak statement of the model problem from Oden's book, Chapter 1 within the PZ environment.
 */

#ifndef TPZMATMODELPROBLEM_H
#define TPZMATMODELPROBLEM_H

#include "TPZMaterial.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatModelProblem : public TPZMaterial
{
    
protected:
   
	
public:
    public:
int ClassId() const override;

    
    void Read(TPZStream& buf, void* context) override {
        TPZMaterial::Read(buf,context);
    }
    
    void Write(TPZStream &buf, int withclassid) const override{
        TPZMaterial::Write(buf,withclassid);
    }
    
    TPZMatModelProblem(int id);
    
    /** @brief Default constructor */
    TPZMatModelProblem();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatModelProblem(const TPZMatModelProblem &mat);
    /** @brief Default destructor */
    virtual ~TPZMatModelProblem();
	
    /** @brief Returns the name of the material */
    virtual std::string Name() override { return "TPZMatModelProblem"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const  override {return 1;}
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const  override { return 1;}
    
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
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
  	
};

#endif

