//
//  TPZElasticBiot.hpp
//  PZ
//
//  Created by Omar on 3/5/17.
//
//

#ifndef TPZElasticBiot_h
#define TPZElasticBiot_h

#include <stdio.h>
#include <iostream>
#include <cmath>
#include "pzmatwithmem.h"
#include "pzdiscgal.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzvec.h"
#include "TPZSimulationData.h"
#include "TPZElasticBiotMemory.h"
#include "pzaxestools.h"


class TPZElasticBiot : public TPZMatWithMem< TPZElasticBiotMemory,TPZDiscontinuousGalerkin> {
    
private:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief material dimension */
    int fdimension;
    
    /** @brief first Lame Parameter */
    REAL flambda;
    
    /** @brief undrained first Lame Parameter */
    REAL flambdau;
    
    /** @brief Bulk modulus */
    REAL fK;
    
    /** @brief undrained Bulk modulus */
    REAL fKu;
    
    /** @brief Second Lame Parameter */
    REAL fmu;
    
    /** @brief constants Biot poroelasticity */
    REAL falpha;
    
    /** @brief Storage coefficient poroelasticity */
    REAL fSe;
    
    /** @brief Intact rock porosity */
    REAL fporosity_0;
    
public:
    
    
    /** @brief Default constructor */
    TPZElasticBiot();
    
    /** @brief Constructor based on a material id */
    TPZElasticBiot(int matid, int dimension);
    
    /** @brief Constructor based on a Biot Poroelasticity  object */
    TPZElasticBiot(const TPZElasticBiot &mat);
    
    /** @brief Constructor based on a Biot Poroelasticity  object */
    TPZElasticBiot &operator=(const TPZElasticBiot &mat)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief Default destructor */
    ~TPZElasticBiot();
    
    /** @brief Set the required data at each integration point */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Set the required data at each integration point */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief Returns the name of the material */
    std::string Name() {
        return "TPZElasticBiot";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return fdimension;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 1;}
    
    virtual TPZMaterial *NewMaterial()
    {
        return new TPZElasticBiot(*this);
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
    
    // Attribute
    
    void SetPorolasticParameters(REAL l, REAL mu, REAL l_u)
    {
        flambda = l;
        fmu = mu;
        flambdau = l_u;
        fK = flambda + (2.0/3.0)*fmu;
        fKu = flambdau + (2.0/3.0)*fmu;
    }
    
    void SetBiotParameters(REAL alpha, REAL Se)
    {
        if(alpha==0){
            std::cout << "Biot constan should be at leats equal to the intact porosity, alpha = " << alpha  << std::endl;
            DebugStop();
        }
        
        falpha = alpha;
        fSe = Se;
    }

    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    void Compute_Sigma(TPZFMatrix<REAL> & S,TPZFMatrix<REAL> & Grad_v);
    
    // Contribute Methods being used
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    // GC basis methods
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    // Reduce basis methods
    void ContributeRB(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeRB_BC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const;
    
    /** @brief Save object data to a stream */
    void Write(TPZStream &buf, int withclassid);
    
    /** @brief Read object data from a stream */
    void Read(TPZStream &buf, void *context);
    
};


#endif /* TPZElasticBiot_h */
