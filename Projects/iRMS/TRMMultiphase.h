//
//  TRMMultiphase.h
//  PZ
//
//  Created by Omar on 5/15/16.
//
//

#ifndef TRMMultiphase_h
#define TRMMultiphase_h

#include <stdio.h>
#include "pzmatwithmem.h"
#include "TRMMemory.h"

#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TRMSimulationData.h"
#include "pzlog.h"


class TRMMultiphase : public TPZMatWithMem<TRMMemory,TPZDiscontinuousGalerkin> {
    
    
private:
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
public:
    
    /** @brief Default constructor */
    TRMMultiphase();
    
    /** @brief Constructor based on a material id */
    TRMMultiphase(int matid);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMMultiphase(const TRMMultiphase &mat);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMMultiphase &operator=(const TRMMultiphase &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TRMMultiphase();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TRMMultiphase";
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return 3;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 1;} // Deprecated, must to be removed

    /** @brief Returns material copied form this object */
    virtual TPZMaterial *NewMaterial()
    {
        return new TRMMultiphase(*this);
    }
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** @brief Returns the number of variables associated with varindex */
    int NSolutionVariables(int var);
    
    /** @brief Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);
    
    /** @brief Returns the solution associated with the var index */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData){
        fSimulationData = SimulationData;
    }
    
    // Contribute Methods being used
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    // one phase flow case
    
    virtual void Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    // two phase flow case
    
    virtual void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    // three phase flow case
    
    virtual void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    virtual void ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const;
    
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid);
    
    /** @brief Read object data from a stream */
    void Read(TPZStream &buf, void *context);
    
    //** @brief Copy the n+1 data to the n data */
    void UpdateMemory();
    
};


#endif /* TRMMultiphase_h */
