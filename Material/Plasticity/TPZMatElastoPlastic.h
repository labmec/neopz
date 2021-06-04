/**
 * @file
 */

#ifndef PZELASTOPLASTIC_H
#define PZELASTOPLASTIC_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatWithMem.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"
#include "TPZPorousElasticResponse.h"

/**
* Implements an elastoplastic material and uses the memory feature to store the hardening variables and strain state.
*/
template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlastic : public TPZMatBase<STATE,
                                               TPZMatSingleSpaceT<STATE>,
                                               TPZMatWithMem<TMEM>,
                                               TPZMatErrorSingleSpace<STATE>>
{
    using TBase = TPZMatBase<STATE,
                             TPZMatSingleSpaceT<STATE>,
                             TPZMatWithMem<TMEM>,
                             TPZMatErrorSingleSpace<STATE>>;
public:
    
    /**
    * Default constructor
    */
    TPZMatElastoPlastic();

    /**
    * Constructor based on material identifier
    */
    TPZMatElastoPlastic(int id);

    /**
    * Copy Constructor
    */
    TPZMatElastoPlastic(const TPZMatElastoPlastic &other);

    /**
    * Desconstructor
    */
    virtual ~TPZMatElastoPlastic();

    /** Sets the plasticity model already with proper parameters */
    virtual void SetPlasticityModel(T & plasticity);

    virtual void UpdateMaterialCoeficients(const TPZVec<REAL> &x,T & plasticity);

    /** Sets the material bulk density */
    virtual void SetBulkDensity(REAL & RhoB);

    /** returns the name of the material*/
    virtual std::string Name() const override;

    /**returns the integrable dimension of the material*/
    virtual int Dimension() const override { return 3; }

    /** returns the number of state variables associated with the material*/
    virtual int NStateVariables() const override{ return 3; }

    /** print out the data associated with the material*/
    virtual void Print(std::ostream &out, const int memory) const;

    /** print out the data associated with the material*/
    virtual void Print(std::ostream &out) const override;

    /**returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name) const override;

    /** returns the number of variables associated with the variable
    indexed by var.  var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var) const override;

    /**returns the solution associated with the var index based on
    * the finite element approximation*/
    virtual void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<REAL> &Solout) override;

    /** Evaluate error between approximate (FEM) and exact solutions.
    *  Method not implemented
    */
    void Errors(const TPZMaterialDataT<STATE>&data,
                        TPZVec<REAL> &values) override;

    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override
    {u_len=3;du_row=3;du_col=3;}
    /**
    * Returns the number of norm errors: 3 (Semi H1, L2 and H1)
    * Method not implemented
    */
    virtual int NEvalErrors() override {return 3;}

    /**
    * It computes a contribution to the stiffness matrix and load vector at one integration point.
    */
    virtual void Contribute(const TPZMaterialDataT<STATE> &data,
                            REAL weight,
                            TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) override;

    /**
    * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
    */
    virtual void ContributeBC(const TPZMaterialDataT<STATE> &data,
                              REAL weight,
                              TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
                              TPZBndCondT<STATE> &bc) override;

    /**
    * It computes a contribution to the residual vector at one integration point.
    */
    virtual void Contribute(const TPZMaterialDataT<STATE> &data,
                            REAL weight, TPZFMatrix<REAL> &ef) override;

    /**
    * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
    */
    virtual void ContributeBC(const TPZMaterialDataT<STATE> &data,
                              REAL weight,
                              TPZFMatrix<REAL> &ef, TPZBndCondT<STATE> &bc) override;

    /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
    * @param data [in]
    * @param Strain [out]
    */
    void ComputeStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &Strain);

    /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
    * @param DeltaStrain [out]
    * @param data [in]
    */
    void ComputeDeltaStrainVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &DeltaStrain);

    /** Evaluates the Stress vector based on an available DSol (solution derivatives set) vector.
    * @param data [in]
    * @param Stress [out]
    */
    void ComputeStressVector(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> &Stress);

    /** Method that checks DEP consistency
    *  @param data [in]
    *  @param DeltaStrain [in]
    */
    void CheckConvergence(const TPZMaterialDataT<STATE> & data,TPZFMatrix<REAL> & DeltaStrain);

    /** Calls the plasticity template aggregate applyStrainComputeDep method
    *  @param data [in]
    *  @param DeltaStrain [in]
    *  @param Stress [out]
    *  @param Dep [out]
    */
    void ApplyDeltaStrainComputeDep(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
                                            TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);

    /** Calls the plasticity template aggregate applyStrain method
    *  @param data [in]
    *  @param DeltaStrain [in]
    *  @param Stress [out]
    */
    void ApplyDeltaStrain(const TPZMaterialDataT<STATE> & data, TPZFMatrix<REAL> & DeltaStrain,
                                            TPZFMatrix<REAL> & Stress);

    /** Applies the tensor in vectorial form to an internally stored direction
    * @param vectorTensor [in]
    * @param Out [out]
    */
    void ApplyDirection(TPZFMatrix<REAL> &vectorTensor, TPZVec<REAL> &Out);

    /** Converts the stress vector onto a symmetric stress tensor
    * @param vectorTensor [in]
    * @param Tensor [out]
    */
    void vectorToTensor(const TPZFMatrix<REAL> & vectorTensor, TPZFMatrix<REAL> & Tensor);

    /** Evaluates the eigenvalues of the tensor vectorTensor (in compact vectorial form)
    * @param vectorTensor [in] compact vectorial form of a symmetrical tensor
    * @param ev [out] evaluated eigenvalues
    */
    void EigenValues(TPZFMatrix<REAL> & vectorTensor, TPZVec<REAL> & ev);

    /** Evaluates the eigenvectors of the tensor vectorTensor (in compact vectorial form)
    * @param vectorTensor [in] compact vectorial form of a symmetrical tensor
    * @param Solout [out] evaluated eigenvector
    * @param direction [in] selected direction (0 to 2)
    */
    void EigenVectors(TPZFMatrix<REAL> &vectorTensor, TPZVec< REAL > &Solout, int direction);

    /**To create another material of the same type*/
    virtual TPZMaterial * NewMaterial() const override;

    /**
    * Unique identifier for serialization purposes
    */
    virtual int ClassId() const override;

    /**
    * Save the element data to a stream
    */
    virtual void Write(TPZStream &buf, int withclassid) const override;

    /**
    * Read the element data from a stream
    */
    virtual void Read(TPZStream &buf, void *context) override;

    /**
    * Sets the tolerance value for post-processing purposes
    */
    virtual void SetTol(const REAL & tol);

    /**
    * Sets the SetBulkDensity of the material
    */
    virtual void SetBulkDensity(const REAL & bulk);
    
    /**
     * Sets the nonlinear elastic response (Porous Elastic Response PER) as predictor during elastoplastic process
     */
    virtual void SetPorousElasticity(TPZPorousElasticResponse & PER);
    
    /**
     * Sets the plasticity model
     */
    void SetPlasticModel(T & plasticity_model);
    
    /**
     * Gets the plasticity model
     */
    virtual T & GetPlasticModel();
    
    /**
     * Gets the nonlinear elastic response (Porous Elastic Response PER) as predictor during elastoplastic process
     */
    virtual TPZPorousElasticResponse & GetPorousElasticity(TPZPorousElasticResponse & PER);

    /**
    * Defining what parameters the material needs. In particular this material needs the
    * evaluation of normal vector for the sake of boundary conditions
    */
    virtual void FillDataRequirements(TPZMaterialData &data) const override;

    /**
     * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
     */
    void FillBoundaryConditionDataRequirements(int type,
                                              TPZMaterialData &data) const override;

    enum ESolutionVar {
        ENone = -1,
        EDisplacementDoF  = 0,
        EDisplacement     = 1,
        EStrain           = 2,
        EStress           = 3,
        EStrainElastic    = 4,
        EStrainPlastic    = 5,
        EYield            = 6,
        EVolHardening     = 7,
        EStrainPValues    = 8,
        EStressPValues    = 9,
        EStrainElasticPValues    = 10,
        EStrainPlasticPValues    = 11,
        EStrainI1           = 12,
        EStressI1           = 13,
        EStrainElasticI1    = 14,
        EStrainPlasticI1    = 15,
        EStrainJ2           = 16,
        EStressJ2           = 17,
        EStrainElasticJ2    = 18,
        EStrainPlasticJ2    = 19,
        EFailureType    = 20,
    };

protected:

    /**
    * gravity acceleration
    */
    TPZManVector<REAL, 3> m_force;

    /**
    * bulk density of rock
    */
    REAL m_rho_bulk;

    /**
    * Post Processing direction
    */
    TPZManVector<REAL,3> m_PostProcessDirection;

    /**
    * Elastoplastic material object instantiation
    * this instantiation avoids several instantiations
    * for each use of this object
    */
    T m_plasticity_model;

    /**
    * Tolerance for post-processing purposes
    */
    REAL m_tol;

    /**
     * Directive that stands for the use of nonlinear elasticity
     */
    bool m_use_non_linear_elasticity_Q;
    
    /**
     * Nonlinear elastic response (Porous Elastic Response PER)
     */
    TPZPorousElasticResponse m_PER;
	
};

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::ClassId() const{
    return Hash("TPZMatElastoPlastic") ^ TPZMatWithMem<TMEM>::ClassId() << 1 ^ T().ClassId() << 2;
}

#endif
