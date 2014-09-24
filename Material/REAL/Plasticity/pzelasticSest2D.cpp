//
//  pzelasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelasticSest2D.h"

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D():TPZElasticityMaterial()
{
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
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id, REAL E, REAL nu, REAL fx, REAL fy, int plainstress):TPZElasticityMaterial(id,E,nu,fx,fy,plainstress)
{
 
    
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id):TPZElasticityMaterial(id)
{
}

/** @brief Copies the data of one TPZElasticityMaterial object to another */
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(const TPZElasticityMaterial &copy):TPZElasticityMaterial(copy)
{
}


