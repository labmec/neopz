//
//  TRMPhaseTransport.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMPhaseTransport__
#define __PZ__TRMPhaseTransport__

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TRMPhaseMemory.h"
#include "TRMPhaseTransport.h"
#include "TRMSimulationData.h"
#include "TRMBuildTransfers.h"


class TRMPhaseTransport : public TPZMatWithMem<TRMPhaseMemory,TPZMaterial> {
    
private:
    
    /** @brief Autopointer of Simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief define the transfer matrices */
    TPZAutoPointer<TRMBuildTransfers> fTransfer;
    
    /** @brief material dimension */
    int fdimension;
    
public:
    
    /** @brief Default constructor */
    TRMPhaseTransport();
    
    /** @brief Constructor based on a material id */
    TRMPhaseTransport(int matid, int dimension);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMPhaseTransport(const TRMPhaseTransport &mat);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMPhaseTransport &operator=(const TRMPhaseTransport &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TRMPhaseTransport();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TRMPhaseTransport";
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 1;} // Deprecated, must to be removed
    
    /** @brief Returns material copied form this object */
    virtual TPZMaterial *NewMaterial()
    {
        return new TRMPhaseTransport(*this);
    }
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** @brief Returns the number of variables associated with varindex */
    int NSolutionVariables(int var);
    
    
    /** @brief Set the simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    /** @brief Get the space generator */
    TPZAutoPointer<TRMSimulationData> SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set the transfer object */
    void SetTransfer(TPZAutoPointer<TRMBuildTransfers> &Transfer)
    {
        fTransfer = Transfer;
    }
    
    /** @brief Get the transfer object */
    TPZAutoPointer<TRMBuildTransfers> Transfer()
    {
        return fTransfer;
    }
    
    
    /** @brief Not used contribute methods */
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){ DebugStop();}
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){DebugStop();}
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){DebugStop();}
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){ DebugStop(); }
    
    
    
    // Contribute Methods being used
    
    /** @brief Returns the solution associated with the var index */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);

    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    
    // Two phase case
    
    void Solution_ab(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_ab(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    
    // Three phase case
    
    void Solution_abc(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute_abc(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    

    
    /**
     * Unique identifier for serialization purposes
     */
    int ClassId() const override;
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context) override;
    
    
    /// Copy the n+1 data to the n data
    void UpdateMemory();


};

#endif /* defined(__PZ__TRMPhaseTransport__) */
