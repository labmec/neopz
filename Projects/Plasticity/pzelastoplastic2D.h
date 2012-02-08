//$Id: pzelastoplastic.h,v 1.17 2009-10-04 05:44:22 diogo Exp $

#ifndef PZELASTOPLASTIC2D_H
#define PZELASTOPLASTIC2D_H


#include "pzmaterial.h"
#include "pzmatwithmem.h"
#include "pzelastoplasticmem.h"
#include "pzporoelastoplasticmem.h"
#include "pzelastoplastic.h"
#include "pzmaterial.h"







/**
 * Implements an elastoplastic material and uses the memory feature to store the damage variables
 * This material works only together with the Plasticity Library.
 */

//typedef TPZMatWithMem<TPZElastoPlasticMem> BASE_MATWITHMEM;

template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlastic2D : public TPZMatElastoPlastic<T,TMEM> //, TPZMatWithMem<TMEM>
{
public:

//enum SOLUTIONVARS{ENone = -1};
/**
 * Default constructor
 */
TPZMatElastoPlastic2D();		

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh. Upon return vectorindex
 *  contains the index of the material object within the
 *  vector
 */
TPZMatElastoPlastic2D(int id ,  int PlaneStrainOrPlaneStress);

/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 *  Upon return vectorindex contains the index of the material
 *  object within the vector
 */
TPZMatElastoPlastic2D(const TPZMatElastoPlastic2D<T,TMEM> &mat);

virtual ~TPZMatElastoPlastic2D();

/** returns the name of the material*/
virtual std::string Name();

/**returns the integrable dimension of the material*/
virtual int Dimension() { return 2; }

/** returns the number of state variables associated with the material*/
virtual int NStateVariables() { return 2; }

/** print out the data associated with the material*/
virtual void Print(std::ostream &out = std::cout, const int memory = 0);


/**returns the solution associated with the var index based on
 * the finite element approximation*/
virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);


/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 */
virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 */
virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);


/** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
 * @param Strain[out] 
 * @param DSol[in]
 */
void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix &DeltaStrain);


/** Calls the plasticity template aggregate applyStrainComputeDep method
 *  @param data[in]
 *  @param DeltaStrain[in]
 *  @param Stress[out]
 *  @param Dep[out]
 */
void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix & DeltaStrain, 
								TPZFMatrix & Stress, TPZFMatrix & Dep);


/**To create another material of the same type*/
virtual TPZAutoPointer<TPZMaterial> NewMaterial();

/**Read data of the material from a istream (file data)*/
//     virtual void SetData(std::istream &data);


/**
 * Unique identifier for serialization purposes
 */
virtual int ClassId() const;

/**
 * Save the element data to a stream
 */
virtual void Write(TPZStream &buf, int withclassid);

/**
 * Read the element data from a stream
 */
virtual void Read(TPZStream &buf, void *context);

/**
 * Defining what parameters the material needs. In particular this material needs the
 * evaluation of normal vector for the sake of boundary conditions
 */
//virtual void FillDataRequirements(TPZMaterialData &data);

protected:


//	  REAL fDeltaT;


int fPlaneStrain;




};

#endif
