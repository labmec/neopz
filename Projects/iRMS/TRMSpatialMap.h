//
//  TRMSpatialMap.hpp
//  PZ
//
//  Created by Omar on 2/19/17.
//
//

#ifndef TRMSpatialMap_hpp
#define TRMSpatialMap_hpp

#include <stdio.h>


#include <stdio.h>
#include "pzdiscgal.h"

class TRMSpatialMap : public TPZDiscontinuousGalerkin {
    
private:
    
    /** @brief material dimension */
    int fdimension;
    
public:
    
    /** @brief Default constructor */
    TRMSpatialMap();
    
    /** @brief Constructor based on a material id */
    TRMSpatialMap(int matid, int dimension);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMSpatialMap(const TRMSpatialMap &mat);
    
    /** @brief Constructor based on a TRMMultiphase object */
    TRMSpatialMap &operator=(const TRMSpatialMap &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TRMSpatialMap();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZMaterialData &data);

    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TRMSpatialMap";
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** @brief Returns the number of state variables associated with the material */
    int NStateVariables() {return 4;}
    
    /** @brief Returns material copied form this object */
    virtual TPZMaterial *NewMaterial()
    {
        return new TRMSpatialMap(*this);
    }
    
    /** @brief Print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** @brief Returns the number of variables associated with varindex */
    int NSolutionVariables(int var);    
    
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
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){ DebugStop(); }
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){ DebugStop(); }
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){ DebugStop(); }
    
    void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
    
    
    
    /**
     * Unique identifier for serialization purposes
     */
    int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);
    
    
};

#endif /* TRMSpatialMap_h */
