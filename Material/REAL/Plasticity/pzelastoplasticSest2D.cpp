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
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D()
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh. Upon return vectorindex
 *  contains the index of the material object within the
 *  vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(int id ,  int PlaneStrainOrPlaneStress)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 *  Upon return vectorindex contains the index of the material
 *  object within the vector
 */
template <class T, class TMEM>
TPZMatElastoPlasticSest2D<T,TMEM>::TPZMatElastoPlasticSest2D(const TPZMatElastoPlastic2D<T,TMEM> &mat)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** returns the name of the material*/
template <class T, class TMEM>
string TPZMatElastoPlasticSest2D<T,TMEM>::Name()
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/**returns the solution associated with the var index based on
 * the finite element approximation*/
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}
/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
 * @param DeltaStrain [out]
 * @param data [in]
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/** Calls the plasticity template aggregate applyStrainComputeDep method
 *  @param data [in]
 *  @param DeltaStrain [in]
 *  @param Stress [out]
 *  @param Dep [out]
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
                                TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/** Calls the plasticity template aggregate applyStrainComputeDep method
 *  @param data [in]
 *  @param DeltaStrain [in]
 *  @param Stress [out]
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
                      TPZFMatrix<REAL> & Stress)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/**To create another material of the same type*/
template <class T, class TMEM>
TPZMaterial * TPZMatElastoPlasticSest2D<T,TMEM>::NewMaterial()
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}


/**
 * Unique identifier for serialization purposes
 */
template <class T, class TMEM>
int TPZMatElastoPlasticSest2D<T,TMEM>::ClassId() const
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * Save the element data to a stream
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Write(TPZStream &buf, int withclassid)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}

/**
 * Read the element data from a stream
 */
template <class T, class TMEM>
void TPZMatElastoPlasticSest2D<T,TMEM>::Read(TPZStream &buf, void *context)
{
    std::cout << "\n this method is not implemented " <<std::endl;
    DebugStop();
}