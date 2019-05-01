//
//  TPZMixedDarcyFlow.h
//  HDiv
//
//  Created by Omar Durán on 3/29/19.
//

#ifndef TPZMixedDarcyFlow_h
#define TPZMixedDarcyFlow_h

#include <stdio.h>
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

class TPZMixedDarcyFlow : public TPZMaterial {
    
private:
    
    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> m_kappa;
    
    /** @brief Inver of absolute permeability */
    TPZFNMatrix<9,REAL> m_kappa_inv;
    
    /** @brief Dimensional factor */
    REAL m_d;
    
    /** @brief Material dimension */
    int m_dim;
    
public:
    
    /** @brief default constructor  */
    TPZMixedDarcyFlow();
    
    /** @brief default desconstructor  */
    ~TPZMixedDarcyFlow();
    
    /** @brief constructor based on material id and dimension */
    TPZMixedDarcyFlow(int mat_id, int dim);
    
    /** @brief constructor based on object copy */
    TPZMixedDarcyFlow(const TPZMixedDarcyFlow &other);
    
    /** @brief return a material object from a this object */
    TPZMaterial * NewMaterial();
    
    /** @brief assignment operator */
    TPZMixedDarcyFlow &operator=(const TPZMixedDarcyFlow &other);
    
    /** @brief return the euclidean dimension of the weak statement */
    int Dimension() const;
    
    /** @brief Sets material dimension */
    void SetDimension(int dim) { m_dim = dim; }
    
    /** @brief return the number of state variables associated with each trial function */
    virtual int NStateVariables() const override;
    
    /** @brief print all material information */
    void Print(std::ostream & out);
    
    /** @brief print material name */
    std::string Name();
    
    /** @brief fill requirements for volumetric contribute methods multiphsycis mesh */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    /** @brief fill requirements for boundary contribute methods multiphsycis mesh */
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec);
    
    /** @brief unique class identifier */
    virtual int ClassId() const;
    
    /** @brief write class in disk */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /** @brief write class from disk */
    void Read(TPZStream &buf, void *context);
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    /** @brief Volumetric contribute with jacobian matrix */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /** @brief Volumetric contribute */
    void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    /** @brief Volumetric contribute without jacobian matrix */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /** @brief Volumetric contribute */
    void Contribute(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef);
    
    /** @brief Boundary contribute without jacobian matrix */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Boundary contribute */
    void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Variable index based on variable naming */
    int VariableIndex(const std::string &name);
    
    /** @brief size of the current variable (1 -> scalar, 3-> vector, 9 ->  Tensor ) */
    int NSolutionVariables(int var);
    
    /** @brief Postprocess required variables multiphysics */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    /// @TODO:: JV document me please
    void SetPermeability(REAL kappa){
        int n_data = m_kappa.Rows();
        m_kappa.Zero();
        m_kappa_inv.Zero();
        for (int i = 0; i < 3; i++) {
            m_kappa(i,i) = kappa;
            m_kappa_inv(i,i) = 1.0/kappa;
        }
    }
    
    /// @TODO:: JV document me please
    void SetPermeability(TPZFNMatrix<9,REAL> & kappa){
        
        if (kappa.Rows() != 3 && kappa.Cols() != 3) {
            DebugStop();
        }
        
        m_kappa = kappa;
        m_kappa.Inverse(m_kappa_inv, ELU);
    }
    
    /// @TODO:: JV document me please
    void SetDimensionalFactor(REAL d){
        m_d = d;
    }
    
};

#endif /* TPZMixedDarcyFlow_h */
