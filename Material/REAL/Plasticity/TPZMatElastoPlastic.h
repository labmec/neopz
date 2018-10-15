/**
 * @file
 */

#ifndef PZELASTOPLASTIC_H
#define PZELASTOPLASTIC_H


#include "TPZMaterial.h"
#include "TPZMatWithMem.h"
#include "TPZElastoPlasticMem.h"
#include "pzporoelastoplasticmem.h"

/**
* Implements an elastoplastic material and uses the memory feature to store the hardening variables and strain state.
*/
template <class T, class TMEM = TPZElastoPlasticMem>
class  TPZMatElastoPlastic : public TPZMatWithMem<TMEM>
{

public:

        enum SOLUTIONVARS{ENone = -1,
          
          EDisplacement=0,//0
          EDisplacementMem,//1
          EPrincipalStress,//2
          ENormalStress,//3
          EShearStress,//4
          ENormalStrain,//5
          EShearStrain,//6
          ENormalPlasticStrain,//7
          EPrincipalStrain,//8
          
          EVolElasticStrain,//9
          EVolPlasticStrain,//10
          EVolTotalStrain,//11
          EAlpha,//12
          EPlasticSteps,//13
          EPlasticSqJ2,//14
          EYield,//15
          EMisesStress,//12
          EI1Stress,//13
          EJ2Stress,//14
          EXStress,//15
          EYStress,//16
          EZStress,//17
          EDisplacementX,//18
          EDisplacementY,//19
          EDisplacementZ,//20
          EDisplacementMemX,//21
          EDisplacementMemY,//22
          EDisplacementMemZ//23
        };

    
        /**
        * Default constructor
        */
        TPZMatElastoPlastic();

        /** Creates a material object and inserts it in the vector of
        *  material pointers of the mesh. Upon return vectorindex
        *  contains the index of the material object within the
        *  vector
        */
        TPZMatElastoPlastic(int id);

        /** Creates a material object based on the referred object and
        *  inserts it in the vector of material pointers of the mesh.
        *  Upon return vectorindex contains the index of the material
        *  object within the vector
        */
        TPZMatElastoPlastic(const TPZMatElastoPlastic &mat);

        virtual ~TPZMatElastoPlastic();

        /** Sets the plasticity model already with proper parameters */
        void SetPlasticity(T & plasticity);

        virtual void UpdateMaterialCoeficients(TPZVec<REAL> &x,T & plasticity);

        /** Sets the material bulk density */
        void SetBulkDensity(REAL & RhoB);

        /** returns the name of the material*/
        virtual std::string Name();

        /**returns the integrable dimension of the material*/
        virtual int Dimension() const { return 3; }

        /** returns the number of state variables associated with the material*/
        virtual int NStateVariables() { return 3; }

        /** print out the data associated with the material*/
        virtual void Print(std::ostream &out, const int memory);

        /** print out the data associated with the material*/
        virtual void Print(std::ostream &out);

        /**returns the variable index associated with the name*/
        virtual int VariableIndex(const std::string &name);

        /** returns the number of variables associated with the variable
        indexed by var.  var is obtained by calling VariableIndex*/
        virtual int NSolutionVariables(int var);

        /**returns the solution associated with the var index based on
        * the finite element approximation*/
        virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

        /** Return the number of components which form the flux function
        * Method not implemented.
        */
        virtual int NFluxes()
        {
         PZError << "TPZMatElastoPlastic::NFluxes() - Method not implemented\n";
         return 0;
        }

        /** Compute the value of the flux function to be used by ZZ error estimator.
        * Method not implemented.
        */
        virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux)
        {
         PZError << "TPZMatElastoPlastic::Flux - Method not implemented\n";
        }

        /** Evaluate error between approximate (FEM) and exact solutions.
        *  Method not implemented
        */
        virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
                          TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                          TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);
        /**
        * Returns the number of norm errors: 3 (Semi H1, L2 and H1)
        * Method not implemented
        */
        virtual int NEvalErrors() {return 3;}

        /**
        * It computes a contribution to the stiffness matrix and load vector at one integration point.
        * @param data [in] stores all input data
        * @param weight [in] is the weight of the integration rule
        * @param ek [out] is the stiffness matrix
        * @param ef [out] is the load vector
        */
        virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);

        /**
        * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
        * @param data [in] stores all input data
        * @param weight [in] is the weight of the integration rule
        * @param ek [out] is the stiffness matrix
        * @param ef [out] is the load vector
        * @param bc [in] is the boundary condition material
        */
        virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

        /**
        * It computes a contribution to the residual vector at one integration point.
        * @param data [in] stores all input data
        * @param weight [in] is the weight of the integration rule
        * @param ef [out] is the residual vector
        */
        virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);

        /**
        * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
        * @param data [in] stores all input data
        * @param weight [in] is the weight of the integration rule
        * @param ef [out] is the load vector
        * @param bc [in] is the boundary condition material
        */
        virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

        /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
        * @param data [in]
        * @param Strain [out]
        */
        void ComputeStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &Strain);

        /** Evaluates the Strain vector based on an available DSol (solution derivatives set) vector.
        * @param DeltaStrain [out]
        * @param data [in]
        */
        void ComputeDeltaStrainVector(TPZMaterialData & data, TPZFMatrix<REAL> &DeltaStrain);

        /** Evaluates the Stress vector based on an available DSol (solution derivatives set) vector.
        * @param data [in]
        * @param Stress [out]
        */
        void ComputeStressVector(TPZMaterialData & data, TPZFMatrix<REAL> &Stress);

        /** Method that checks DEP consistency
        *  @param data [in]
        *  @param DeltaStrain [in]
        */
        void CheckConvergence(TPZMaterialData & data,TPZFMatrix<REAL> & DeltaStrain);

        /** Calls the plasticity template aggregate applyStrainComputeDep method
        *  @param data [in]
        *  @param DeltaStrain [in]
        *  @param Stress [out]
        *  @param Dep [out]
        */
        void ApplyDeltaStrainComputeDep(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
                                                TPZFMatrix<REAL> & Stress, TPZFMatrix<REAL> & Dep);

        /** Calls the plasticity template aggregate applyStrain method
        *  @param data [in]
        *  @param DeltaStrain [in]
        *  @param Stress [out]
        */
        void ApplyDeltaStrain(TPZMaterialData & data, TPZFMatrix<REAL> & DeltaStrain,
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
        virtual TPZMaterial * NewMaterial();

        /**
        * Unique identifier for serialization purposes
        */
        virtual int ClassId() const;

        /**
        * Save the element data to a stream
        */
        virtual void Write(TPZStream &buf, int withclassid) const;

        /**
        * Read the element data from a stream
        */
        virtual void Read(TPZStream &buf, void *context);

        /**
        * Sets the tolerance value for post-processing purposes
        */
        void SetTol(const REAL & tol);

        /**
        * Sets the SetBulkDensity of the material
        */
        void SetBulkDensity(const REAL & bulk);

        /**
        * Defining what parameters the material needs. In particular this material needs the
        * evaluation of normal vector for the sake of boundary conditions
        */
        virtual void FillDataRequirements(TPZMaterialData &data);

        /**
         * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
         */
        virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);

    protected:

        /**
        * gravity acceleration
        */
        TPZManVector<REAL, 3> fForce;

        /**
        * bulk density of rock
        */
        REAL fRhoB;

        /**
        * Post Processing direction
        */
        TPZManVector<REAL,3> fPostProcessDirection;

        /**
        * Elastoplastic material object instantiation
        * this instantiation avoids several instantiations
        * for each use of this object
        */

        T fPlasticity;

        /**
        * Tolerance for post-processing purposes
        */
        REAL fTol;
	
};

template <class T, class TMEM>
int TPZMatElastoPlastic<T,TMEM>::ClassId() const{
    return Hash("TPZMatElastoPlastic") ^ TPZMatWithMem<TMEM>::ClassId() << 1 ^ T().ClassId() << 2;
}

#endif
