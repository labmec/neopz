//
//  TPZReynoldsFlow.h
//  PZ
//
//  Created by Cesar Lucci on 06/02/13.
//
//

#ifndef __PZ__TPZReynoldsFlow__
#define __PZ__TPZReynoldsFlow__

#include <iostream>


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
class  TPZReynoldsFlow : public TPZMaterial
{
public:
    
    TPZReynoldsFlow();
    TPZReynoldsFlow(int matId, REAL visc, REAL deltaT, REAL StaticPotential);
    TPZReynoldsFlow(const TPZReynoldsFlow &cp);
    virtual ~TPZReynoldsFlow();
    
    void SetNComputation()
    {
        f_nplus1Computation = false;
    }
    void SetNplus1Computation()
    {
        f_nplus1Computation = true;
    }

    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override;
    
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override;
    
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
//    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    //Sera utilizado qdo for multifisico!!!

    
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
   // virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
     //Sera utilizado qdo for multifisico!!!
    

    
	    
    /** @brief To create another material of the same type*/
    virtual TPZMaterial * NewMaterial() override;
    
    /** @brief Unique identifier for serialization purposes */
    public:
virtual int ClassId() const override;

    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context) override;
    
private:
    REAL f_visc;
    REAL f_deltaT;
    REAL f_staticPotential;
    bool f_nplus1Computation;
};

#endif



