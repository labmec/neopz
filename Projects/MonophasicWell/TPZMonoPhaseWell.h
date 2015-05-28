//
//  TPZMonoPhaseWell.h
//  PZ
//
//  Created by omar duran on 25/05/2015.
//
//

#ifndef __PZ__TPZMonoPhaseWell__
#define __PZ__TPZMonoPhaseWell__

#include <stdio.h>
#include "pzmaterial.h"


/**
 * @ingroup material
 * @brief Implements a Multiphycis mixed formulation for monophasic wel
 */
class TPZMonoPhaseWell : public TPZMaterial {

public:
    
    /**
     * @name Constructors and Destructors methods
     * @{
     */
    
    /** @brief Constructor for given convection */
    TPZMonoPhaseWell(int id);
    
    /** @brief Destructor */
    ~TPZMonoPhaseWell();
    
    /** @brief Copy constructor */
    TPZMonoPhaseWell(TPZMonoPhaseWell & copy);
    
    /** @} */
    
private:
    
    /** @brief Returns the name of the material */
    std::string Name() { return "TPZMonoPhaseWell";}
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return 1;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 1;}
    
    
    /**
     * @name Contribute methods (weak formulation)
     * @{
     */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
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
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /** @} */
    
    /** @name Contribute methods not used here
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
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
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }

    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }

    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }

    /** @} */
    
    
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * @since April 10, 2007
     */
    /**
     * Contribute method. Here, in base class, all requirements are considered as necessary.
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    void FillDataRequirements(TPZMaterialData &data);
    
    /**
     * @brief Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered as necessary.
     * Each derived class may optimize performance by selecting only the necessary data.
     */
    void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);

	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /** @brief Reads data of the material from a istream (file data) */
    void SetData(std::istream &data);
    
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const;
    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context);
    
    /** @} */
    
    
};


#endif /* defined(__PZ__TPZMonoPhaseWell__) */
