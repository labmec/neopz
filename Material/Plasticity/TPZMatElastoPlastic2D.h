/**
 * @file
 */

#ifndef PZELASTOPLASTIC2D_H
#define PZELASTOPLASTIC2D_H

#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "TPZMatElastoPlastic.h"
#include "TPZMaterial.h"

/**
 * Implements an elastoplastic material and uses the memory feature to store the damage variables
 * This material works only together with the Plasticity Library.
 */

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
	virtual std::string Name() const override;
	
	/**returns the integrable dimension of the material*/
	virtual int Dimension() const override { return 2; }
	
	/** returns the number of state variables associated with the material*/
	virtual int NStateVariables() const override { return 2; }
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out) const override;

	/** print out the data associated with the material*/
	virtual void Print(std::ostream &out, const int memory) const override;
	
	
	/**returns the solution associated with the var index based on
	 * the finite element approximation*/
	virtual void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ef) override;
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 */
	virtual void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) override;
	
    /**
     * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
     */
    virtual void FillBoundaryConditionDataRequirements(int type,TPZMaterialData &data) const override;
    
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 */
	virtual void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     */
    virtual void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

	
	/** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
	 * @param DeltaStrain [out]
	 * @param data [in]
	 */
	virtual void ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain);
    
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 *  @param Dep [out]
	 */
	virtual void ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);
	
	/** Calls the plasticity template aggregate applyStrainComputeDep method
	 *  @param data [in]
	 *  @param DeltaStrain [in]
	 *  @param Stress [out]
	 */
	virtual void ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
									TPZFMatrix<REAL> & Stress);
	
	
	/**To create another material of the same type*/
	virtual TPZMaterial * NewMaterial() const override;
	
	
	/**
	 * Unique identifier for serialization purposes
	 */
	public:
virtual int ClassId() const override;

	
	/**
	 * Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/**
	 * Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context) override;
    
	
protected:
	
	
	int fPlaneStrain;
    
};

template <class T, class TMEM>
int TPZMatElastoPlastic2D<T, TMEM>::ClassId() const{
    return Hash("TPZMatElastoPlastic2D") ^ TPZMatElastoPlastic<T,TMEM>::ClassId() << 1;
}

#endif
