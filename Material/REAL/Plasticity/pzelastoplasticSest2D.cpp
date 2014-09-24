//
//  pzelastoplasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelastoplasticSest2D.h"
//enum SOLUTIONVARS{ENone = -1};
/**
 * Default constructor
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D():TPZMatElastoPlastic2D<T, TMEM>()
{
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh. Upon return vectorindex
 *  contains the index of the material object within the
 *  vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress):TPZMatElastoPlastic2D<T, TMEM>(id,PlaneStrainOrPlaneStress)
{
}

/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 *  Upon return vectorindex contains the index of the material
 *  object within the vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(const TPZMatElastoPlastic2D<T,TMEM> &mat):TPZMatElastoPlastic2D<T, TMEM>(mat)
{
    
}

