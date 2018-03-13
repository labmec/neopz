//
//  TISSMaterial.h
//  PZ
//
//  Created by Thiago Dias dos Santos on 8/3/15.
//
//

#ifndef __PZ__TISSMaterial__
#define __PZ__TISSMaterial__


#include "TPZMaterial.h"
#include "pzmaterialdata.h"

/**
 * @ingroup material
 * @brief This class implements a simple material without formulation; it is used to load the solution from ISSM
 */
class  TISSMaterial : public TPZMaterial
{
    
protected:
    
    
public:
    
    TISSMaterial(int matid);
    
    /** @brief Default constructor */
    TISSMaterial();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TISSMaterial(const TISSMaterial &mat);
    
    /** @brief Default destructor */
    virtual ~TISSMaterial();
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TISSMaterial"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() { return 8;}
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
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

#endif /* defined(__PZ__TISSMaterial__) */
