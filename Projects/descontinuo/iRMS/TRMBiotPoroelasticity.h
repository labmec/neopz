//
//  TRMBiotPoroelasticity.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMBiotPoroelasticity__
#define __PZ__TRMBiotPoroelasticity__

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "TRMMemory.h"


#include "pzfmatrix.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TRMSimulationData.h"
#include "pzlog.h"

class TRMBiotPoroelasticity : public TPZMatWithMem<TRMMemory,TPZMaterial> {
    
private:
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief material dimension */
    int fdimension;
    
public:
    
    /** @brief Default constructor */
    TRMBiotPoroelasticity();
    
    /** @brief Constructor based on a material id */
    TRMBiotPoroelasticity(int matid, int dimension);
    
    /** @brief Default desconstructor */
    ~TRMBiotPoroelasticity();
    
    /** @brief Copy constructor $ */
    TRMBiotPoroelasticity(const TRMBiotPoroelasticity& other);
    
    /** @brief Copy assignemnt operator $ */
    TRMBiotPoroelasticity & operator = (const TRMBiotPoroelasticity& other);
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TRMBiotPoroelasticity";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {
        return fdimension;
    } // for 2d Cases Plane Strain is being considered
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TRMBiotPoroelasticity(*this);
    }
    
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
   
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

    // Contribute Methods being used
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const override;
    
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Read object data from a stream */
    void Read(TPZStream &buf, void *context) override;
    
    //** @brief Copy the n+1 data to the n data */
    void UpdateMemory();
    
    
};

#endif /* defined(__PZ__TRMBiotPoroelasticity__) */
