//
//  TNFRElasticMaterial.h
//  pz
//
//  Created by Omar Durán on 7/2/18.
//

#ifndef TNFRElasticMaterial_h
#define TNFRElasticMaterial_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TNRFElasticMemory.h"
#include "TPZTensor.h"

class TNFRElasticMaterial : public TPZMatWithMem<TNRFElasticMemory,TPZMaterial> {
  
public:
    
    /// Constructor
    TNFRElasticMaterial();
    
    /// Constructor based on material identifier
    TNFRElasticMaterial(int id);
    
    /// Destructor
    ~TNFRElasticMaterial();
    
    /// material dimension
    virtual int Dimension() const ;
    
    /// number of state variables
    virtual int NStateVariables();
    
    /// setting data on domain
    void FillDataRequirements(TPZMaterialData &data);
    
    ///  setting data on boundary
    void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data);
    
    /// jacobian contribution
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /// residual contribution
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name);
    
    /// Returns the number of variables associated with provided index var
    int NSolutionVariables(int var);
    
    /// post-processing
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    /// Set drained Lamé parameters
    void SetLameParameters(REAL lambda, REAL mu);
    
private:
    
    /// First lame parameter
    REAL m_lambda;
    
    ///  Second lame parameter
    REAL m_mu;
    
};


#endif /* TNFRElasticMaterial_h */
