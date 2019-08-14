/**
 * @file TPZMatExSimples2D.h
 * @brief Header file for class TPZMatExSimples2D.\n
 * It implements the weak statement of the model problem from Oden's book, Chapter 1 within the PZ environment.
 */

#ifndef TPZMATEXSIMPLES2D_H
#define TPZMATEXSIMPLES2D_H

#include "TPZMaterial.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement of the model problem from Oden's book, Chapter 1, within the PZ environment
 */
class  TPZMatExSimples2D : public TPZMaterial
{
    
protected:
   
	
public:
    
    TPZMatExSimples2D(int id);
    
    /** @brief Default constructor */
    TPZMatExSimples2D();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatExSimples2D(const TPZMatExSimples2D &mat);
    /** @brief Default destructor */
    virtual ~TPZMatExSimples2D();
	
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMatExSimples2D"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const { return 1;}
    
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
		
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
  	
};

#endif

