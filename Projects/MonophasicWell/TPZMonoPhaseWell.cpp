//
//  TPZMonoPhaseWell.cpp
//  PZ
//
//  Created by omar duran on 25/05/2015.
//
//

#include "TPZMonoPhaseWell.h"



/** @brief Constructor */
TPZMonoPhaseWell::TPZMonoPhaseWell(int id): TPZMaterial(id)
{
    
}

/** @brief Destructor */
TPZMonoPhaseWell::~TPZMonoPhaseWell()
{
    
}


/** @brief Copy constructor */
TPZMonoPhaseWell::TPZMonoPhaseWell(TPZMonoPhaseWell & copy){
    
}


/**
 * @name Contribute methods (weak formulation)
 * @{
 */

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
 * @param datavec [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
void TPZMonoPhaseWell::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
}

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
void TPZMonoPhaseWell::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
}

/** @} */

/**
 * @brief Fill material data parameter with necessary requirements for the
 * @since April 10, 2007
 */
/**
 * Contribute method. Here, in base class, all requirements are considered as necessary.
 * Each derived class may optimize performance by selecting only the necessary data.
 */
void TPZMonoPhaseWell::FillDataRequirements(TPZMaterialData &data){
    
}

/**
 * @brief Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered as necessary.
 * Each derived class may optimize performance by selecting only the necessary data.
 */
void TPZMonoPhaseWell::FillDataRequirements(TPZVec<TPZMaterialData > &datavec){
    
}




/** @brief Print out the data associated with the material */
void TPZMonoPhaseWell::Print(std::ostream &out){
    
}

/** @brief Returns the variable index associated with the name */
int TPZMonoPhaseWell::VariableIndex(const std::string &name){
    
    return 0;
    
}

int TPZMonoPhaseWell::NSolutionVariables(int var){
    
    return 0;
}

void TPZMonoPhaseWell::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
}

/** @brief Reads data of the material from a istream (file data) */
void TPZMonoPhaseWell::SetData(std::istream &data)
{
    
}

/** @{
 * @name Save and Load methods
 */

/** @brief Unique identifier for serialization purposes */
int TPZMonoPhaseWell::ClassId() const {

    return 9999888;
}

/** @brief Saves the element data to a stream */
void TPZMonoPhaseWell::Write(TPZStream &buf, int withclassid){
    
}

/** @brief Reads the element data from a stream */
void TPZMonoPhaseWell::Read(TPZStream &buf, void *context){
    
}

/** @} */