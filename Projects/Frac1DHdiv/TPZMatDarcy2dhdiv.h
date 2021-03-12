//
//  TPZMatDarcy2dhdiv.h
//  PZ
//
//  Created by nathan and omar on 9/3/14.
//
//

#ifndef __PZ__TPZMatDarcy2dhdiv__
#define __PZ__TPZMatDarcy2dhdiv__

#include <iostream>
#include "TPZMaterial.h"
#include "tpzautopointer.h"
#include "TPZFracData.h"

/**
 * @ingroup material
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Material to solve a 2d mixed formulation for darcy flow
 * @brief Here is used Hdiv for flux and L2 for pressure
 * @brief DOCUMENTATION OF WEAK FORMULATION IN LYX LOCATED AT THE SVN REPOSITORY
 */
class TPZMatDarcy2dhdiv : public TPZMaterial {
    
protected:
    
    /** @brief Problem dimension */
    int fDim;
  
    /** @brief flag to set if contribute doesnt do anything */
    bool fNotContribute;
  
public:
    
    /** @brief Default Constructor */
    TPZMatDarcy2dhdiv();
    
    /** @brief Constructor with matid */    
    TPZMatDarcy2dhdiv(int matid);
    
    /** @brief Destructor */      
    virtual ~TPZMatDarcy2dhdiv();
    
    /** @brief copy constructor */
    TPZMatDarcy2dhdiv(const TPZMatDarcy2dhdiv &copy);
    
    /** @brief operator equal */    
    TPZMatDarcy2dhdiv &operator=(const TPZMatDarcy2dhdiv &copy);
    
    /** @brief Print Method */
    virtual void Print(std::ostream & out);
    
    /** @brief Name of the material */
    virtual std::string Name() { return "TPZMatDarcy2dhdiv"; }
    
    /** @brief Returns the integrable dimension */    
    virtual int Dimension() const;
    
    /** @brief Return the number of state variables */
    virtual int NStateVariables() const { return 1; }
    
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef);
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Used contribute methods */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec,TPZVec<TPZMaterialData> &dataleftvec,TPZVec<TPZMaterialData> &datarightvec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /** @brief Bc Contribution for flux state variable */
    virtual void ApplyQnD       (TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /** @brief Bc Contribution for pressure state variable */
    virtual void ApplyPN        (TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    /** @brief Fill material data parameter with necessary requirements for the Contribute method*/    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Fill material data parameter with necessary requirements for the ContributeBC method*/    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    /** @brief Returns the variable index associated with the name */    
    virtual int VariableIndex(const std::string &name);
    
    /** @brief Returns the number of variables associated with the variable indexed by var */
    virtual int NSolutionVariables(int var);
    
    /** @brief Calculates a solution given datavec*/    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
private:
    
    /** @brief Data of the simulation */
    TPZAutoPointer<TPZFracData> fData;
    
public:
    
    /** @brief Sets data of the simulation */
    void SetSimulationData(TPZAutoPointer<TPZFracData> Data) { fData = Data;}
  
    /** @brief Sets if doesnt wanna calculate the contribution. Used to speed up the memory update of the fracture */
    void SetNotContribute(bool setNotCont = true);
  
};


#endif /* defined(__PZ__TPZMatDarcy2dhdiv__) */


