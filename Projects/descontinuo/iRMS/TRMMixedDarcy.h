//
//  TRMMixedDarcy.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMMixedDarcy__
#define __PZ__TRMMixedDarcy__

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TRMMemory.h"


#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TRMSimulationData.h"
#include "pzlog.h"

class TRMMixedDarcy : public TPZMatWithMem<TRMMemory,TPZMaterial> {
    
private:
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;

    /** @brief material dimension */
    int fdimension;

    /** @brief Material parameter for non symetric system  */
    REAL fnon_symetric;
    
public:
    
    
    /** @brief Default constructor */
    TRMMixedDarcy();
    
    /** @brief Constructor based on a material id */
    TRMMixedDarcy(int matid, int dimension);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMMixedDarcy(const TRMMixedDarcy &mat);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMMixedDarcy &operator=(const TRMMixedDarcy &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TRMMixedDarcy();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TRMMixedDarcy";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 1;} // for hdiv are 3
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TRMMixedDarcy(*this);
    }
    
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    // Contribute Methods being used
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    // one phase flow case
    
    void Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC_a(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void Solution_a(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    // two phase flow case
    
    void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    void ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    void Solution_ab(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    // three phase flow case
    
    void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    void ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    void Solution_abc(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const override;
    
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Read object data from a stream */
    void Read(TPZStream &buf, void *context) override;
    
    //** @brief Copy the n+1 data to the n data */
    void UpdateMemory();
};

#endif /* defined(__PZ__TRMMixedDarcy__) */
