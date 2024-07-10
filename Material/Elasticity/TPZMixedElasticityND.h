/**
 * @file
 * @brief Contains the TPZMixedElasticityND class which implements an N-dimensional elastic material.
 */

#ifndef MIXEDELASMATNDHPP
#define MIXEDELASMATNDHPP

#include <iostream>

#include "TPZMaterial.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @ingroup material
 * @brief This class implements an N-dimensional elastic material
 */
class TPZMixedElasticityND : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
TPZMatErrorCombinedSpaces<STATE>> {

// type alias to improve constructor readability
using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,
    TPZMatErrorCombinedSpaces<STATE> >;

public:
    
    struct TElasticityAtPoint
    {
        /** @brief Elasticity modulus */
        REAL fE;
        
        /** @brief Poison coefficient */
        REAL fnu;
        
        /** @brief first Lame Parameter */
        REAL flambda;
        
        /** @brief Second Lame Parameter */
        REAL fmu;

        TElasticityAtPoint(REAL E, REAL nu) : fE(E), fnu(nu)
        {
            flambda = (E * nu) / ((1 + nu)*(1 - 2 * nu));
            fmu = E / (2 * (1 + nu));
        }
    };



public:

    enum MVoigt {
        Exx, Eyy, Exy, Eyx, Exz, Ezx, Eyz, Ezy, Ezz
    };

    /** @brief Default constructor */
    TPZMixedElasticityND();
    /** 
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$ 
     * @param fy forcing function \f$ -y = fy \f$
     * @param planestress \f$ planestress = 1 \f$ indicates use of planestress
     * @param Dimension of the problem
     */
    TPZMixedElasticityND(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1, int Dimension = 1);

    TPZMixedElasticityND(int id, int dimension = 2);

    /** @brief Copies the data of one TPZMixedElasticityND object to another */
    TPZMixedElasticityND(const TPZMixedElasticityND &copy);


    virtual int VariableIndex(const std::string &name) const override;


    virtual int NSolutionVariables(int var) const override;

    /** index of Stress */
    int SIndex() {
        return 0;
    }

    /** index of Displacement */
    int UIndex() {
        return 1;
    }

    /** index of Rotation */
    int PIndex() {
        return 2;
    }

    virtual int NEvalErrors() const override;
    
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialDataT<STATE>> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);

    void ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elastb);

    void ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress,  TElasticityAtPoint &elast);

    void ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast);

    /** @brief Creates a new material from the current object   ??*/
    virtual TPZMaterial * NewMaterial()  const override {
        return new TPZMixedElasticityND(*this);
    }

    /** @brief Default destructor */
    virtual ~TPZMixedElasticityND();

    /**
     * @brief Set parameters of elastic material:
     * @param Eyoung Young's elasticity modulus
     * @param nu Poisson's coefficient
     */
    void SetElasticParameters(REAL Eyoung, REAL nu) {
        this->SetElasticity(Eyoung, nu);
    }

    /** @brief Set elasticity parameters */
    void SetElasticity(REAL E, REAL nu) {
        fE_const = E; // Young modulus
        fnu_const = nu; // poisson coefficient
        flambda_const = (E * nu) / ((1 + nu)*(1 - 2 * nu));
        fmu_const = E / (2 * (1 + nu));

    }
    
    void GetElasticity(REAL &E, REAL &nu)
    {
        E = fE_const;
        nu = fnu_const;
    }

    /// Set a variable elasticity and poisson coefficient
//    void SetElasticityFunction(TPZAutoPointer<TPZFunction<STATE> > func)
//    {
//        fElasticity = func;
//    }
    void SetElasticityFunction(std::function<void(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)>& func)
    {
        fElasticity = func;
    }
    
    /** @brief Set elasticity parameters in the form of Lame's parameters*/
    void SetLameParameters(REAL Lambda, REAL mu) {
        fE_const = (mu * (3.0 * Lambda + 2.0 * mu)) / (Lambda + mu);
        fnu_const = (Lambda) / (2 * (Lambda + mu));

        flambda_const = Lambda;
        fmu_const = mu;
    }

    /// set the material configuration to AxisSymmetric

    void SetAxisSymmetric() {
        if(fDimension == 3) DebugStop();
        fAxisSymmetric = 1;
    }
    /// Set the material configuration to plane strain

    void SetPlaneStrain() {
        if(fDimension != 2) DebugStop();
        fPlaneStress = 0;
        // fPlaneStrain = 1;
    }

    /// Set the material configuration to plane stress

    void SetPlaneStress() {
        if(fDimension != 2) DebugStop();
        fPlaneStress = 1;
        // fPlaneStrain = 0;
    }

    /** @brief Set forcing function */
    void SetBodyForce(REAL fx, REAL fy, REAL fz = 0.) {
        fForce[0] = fx;
        fForce[1] = fy;
        fForce[2] = fz;
    }

    /** @brief Returns the model dimension */
    int Dimension() const  override {
        return fDimension;
    }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override;

    /** @brief Print the material data*/
    virtual void Print(std::ostream & out = std::cout) const override;

    /** @brief Returns the material name*/
    virtual std::string Name() const override {
        return "TPZMixedElasticityND";
    }
    
    CSTATE GetMaxComplianceEigenvalue(TPZVec<REAL> &x) const;
    
    /** @name Contribute methods */
    /** @{ */

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;


    /** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
    virtual void Contribute_3spaces(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /** @brief Calculates the element stiffness matrix using 5 spaces - 3 from Contribute_3spaces() + Rigid body motions, and distributed forces */
    virtual void Contribute_5spaces(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

    // contribute a distributed flux and average displacement/rotation
    void ContributeRigidBodyMode(const TPZVec<TPZMaterialDataT<STATE>> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, int first_eq, int RB_space);
    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    //    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ef)
    //    {
    //        DebugStop();
    //    }


    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;


    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    STATE Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);

    /// Transform a tensor to a Voigt notation
    void ToVoigt(const TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt) const;

    /// Transform a Voigt notation to a tensor
    void FromVoigt(const TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S) const;

    /// Transforms a second-order tensor to a Mandel notation
    void MandelFromSecondOrderTensor(const TPZFMatrix<STATE> &S, TPZVec<STATE> &SMandel) const;

    /// Transforms a Mandel notation to a second-order tensor
    void SecondOrderTensorFromMandel(const TPZVec<STATE> &SMandel, TPZFMatrix<STATE> &S) const;

    /// Transforms a symmetric second-order tensor to a Mandel notation
    void MandelFromSymmetricSecondOrderTensor(const TPZFMatrix<STATE> &S, TPZVec<STATE> &SMandel) const;

    /// Transforms a Mandel notation to a second-order symmetric tensor
    void SymmetricSecondOrderTensorFromMandel(const TPZVec<STATE> &SMandel, TPZFMatrix<STATE> &S) const;

    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    STATE InnerVec(const TPZVec<STATE> &S, const TPZVec<STATE> &T);

    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU);

    /** transpose of the tensor GradU */
    STATE Transpose(TPZFMatrix<REAL> &GradU);

    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi);

    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);

public:


    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] datavec material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) override;


    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;

    virtual int ClassId() const override;

    virtual void Read(TPZStream &buf, void *context) override;

    virtual void Write(TPZStream &buf, int withclassid) const override;



protected:
    /** @brief Forcing vector */
    TPZManVector<REAL, 3> fForce = TPZManVector<REAL, 3>(3, 0.);

    /** @brief Plane stress flag */
    int fPlaneStress = 1;
    
    /** @brief Plane strain flag */
    // int fPlaneStrain = 0;
    
    /// dimension of the material
    int fDimension = 2;

    /** Elasticity function */
//    TPZAutoPointer<TPZFunction<STATE> > fElasticity;
    std::function<void(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> fElasticity{nullptr};
    
    /** @brief Elasticity modulus */
    REAL fE_const;
    
    /** @brief Poison coefficient */
    REAL fnu_const;
    
    /** @brief first Lame Parameter */
    REAL flambda_const;
    
    /** @brief Second Lame Parameter */
    REAL fmu_const;
    
    // Matrix A
    TPZFMatrix<STATE> fMatrixA;

    /// flag indicates axix-AxisSymmetric 
    bool fAxisSymmetric = false;

};

#endif
