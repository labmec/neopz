//
//  pzelasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelasticSest2D.h"

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D()
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}
/**
 * @brief Creates an elastic material with:
 * @param id material id
 * @param E elasticity modulus
 * @param nu poisson coefficient
 * @param fx forcing function \f$ -x = fx \f$
 * @param fy forcing function \f$ -y = fy \f$
 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
 */
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id, REAL E, REAL nu, REAL fx, REAL fy, int plainstress)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** @brief Copies the data of one TPZElasticityMaterial object to another */
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(const TPZElasticityMaterial &copy)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** @brief Returns the number of state variables associated with the material */
int TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D::NStateVariables()
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** @name Contribute methods */
/** @{ */

/** @brief Calculates the element stiffness matrix */
void TPZElasticityMaterialSest2D::Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

void TPZElasticityMaterialSest2D::ContributeVecShape(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/** @brief Applies the element boundary conditions */
void TPZElasticityMaterialSest2D::ContributeBC(TPZMaterialData &data,REAL weight,
                          TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

void TPZElasticityMaterialSest2D::ContributeVecShapeBC(TPZMaterialData &data,REAL weight,
                          TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


//virtual void FillDataRequirements(TPZMaterialData &data);
void TPZElasticityMaterialSest2D::FillDataRequirements(TPZMaterialData &data)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}
void TPZElasticityMaterialSest2D::FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/** @} */

/** @brief Returns the variable index associated with the name */
int TPZElasticityMaterialSest2D::VariableIndex(const std::string &name)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * @brief Returns the number of variables associated with the variable indexed by var.
 */
int TPZElasticityMaterialSest2D::NSolutionVariables(int var)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZElasticityMaterialSest2D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/** @brief Computes the value of the flux function to be used by ZZ error estimator */
void TPZElasticityMaterialSest2D::Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * @brief Computes the error due to the difference between the interpolated flux \n
 * and the flux computed based on the derivative of the solution
 */
void TPZElasticityMaterialSest2D::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
            TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
            TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** @brief Set PresStress Tensor */
void SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy, REAL Sigzz)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

int TPZElasticityMaterialSest2D::ClassId() const
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

void TPZElasticityMaterialSest2D::Read(TPZStream &buf, void *context)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

void TPZElasticityMaterialSest2D::Write(TPZStream &buf, int withclassid)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


