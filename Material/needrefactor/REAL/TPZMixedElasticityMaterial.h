/**
 * @file
 * @brief Contains the TPZMixedElasticityMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef MIXEDELASMATHPP
#define MIXEDELASMATHPP

#include <iostream>

#include "TPZMaterial.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZMixedElasticityMaterial : public TPZMaterial {

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
        Exx, Eyy, Exy, Eyx
    };

    /** @brief Default constructor */
    TPZMixedElasticityMaterial();
    /** 
     * @brief Creates an elastic material with:
     * @param id material id
     * @param E elasticity modulus
     * @param nu poisson coefficient
     * @param fx forcing function \f$ -x = fx \f$ 
     * @param fy forcing function \f$ -y = fy \f$
     * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
     */
    TPZMixedElasticityMaterial(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress = 1, int fDimension = 1);

    /// dimension of the material

    TPZMixedElasticityMaterial(int id);

    /** @brief Copies the data of one TPZMixedElasticityMaterial object to another */
    TPZMixedElasticityMaterial(const TPZMixedElasticityMaterial &copy);


    virtual int VariableIndex(const std::string &name) override;


    virtual int NSolutionVariables(int var) override;

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

    virtual int NEvalErrors() override;
    
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);

    void ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elastb);

    void ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress,  TElasticityAtPoint &elast);

    void ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast);

    /** @brief Creates a new material from the current object   ??*/
    virtual TPZMaterial * NewMaterial()  override {
        return new TPZMixedElasticityMaterial(*this);
    }

    /** @brief Default destructor */
    virtual ~TPZMixedElasticityMaterial();

    /**
     * @brief Set parameters of elastic material:
     * @param First  Lame Parameter Lambda
     * @param Second Lame Parameter Mu -> G
     * @param fx forcing function \f$ -x = 0 \f$
     * @param fy forcing function \f$ -y = 0 \f$
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

    /// Set a variable elasticity and poisson coefficient
    void SetElasticityFunction(TPZAutoPointer<TPZFunction<STATE> > func)
    {
        fElasticity = func;
    }
    
    /** @brief Set elasticity parameters */
    void SetLameParameters(REAL Lambda, REAL mu) {
        fE_const = (mu * (3.0 * Lambda + 2.0 * mu)) / (Lambda + mu);
        fnu_const = (Lambda) / (2 * (Lambda + mu));

        flambda_const = Lambda;
        fmu_const = mu;
    }

    /// set the material configuration to AxisSymmetric

    void SetAxisSymmetric() {
        fAxisSymmetric = 1;
    }
    /// Set the material configuration to plane strain

    void SetPlaneStrain() {
        fPlaneStress = 0;
    }

    /// Set the material configuration to plane stress

    void SetPlaneStress() {
        fPlaneStress = 1;
    }

    /** @brief Set forcing function */
    void SetBodyForce(REAL fx, REAL fy) {
        fForce[0] = fx;
        fForce[1] = fy;
        fForce[2] = 0.;
    }

    /** @brief Returns the model dimension */
    int Dimension() const  override {
        return 2;
    }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override;

    /** @brief Print the material data*/
    virtual void Print(std::ostream & out = std::cout) override;

    /** @brief Returns the material name*/
    virtual std::string Name()  override {
        return "TPZMixedElasticityMaterial";
    }

    /** @brief Returns the number of components which form the flux function */
    virtual short NumberOfFluxes() {
        return 3;
    }

    /** @name Contribute methods */
    /** @{ */

    /** @brief Calculates the element stiffness matrix */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    /** @brief Calculates the element stiffness matrix - simulate compaction as aditional variable */
    //    virtual void Contribute(TPZVec<TPZMaterialData> &data, REAL weight, TPZFMatrix<STATE> &ef)
    //    {
    //        DebugStop();
    //    }

    /** @brief Calculates the element stiffness matrix */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override {
        TPZMaterial::Contribute(data, weight, ef);
    }


    /** @brief Applies the element boundary conditions Mixed */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;


    /** @brief Applies the element boundary conditions */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight,
            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;

    /** @brief Applies the element boundary conditions */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight,
            TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
        TPZMaterial::ContributeBC(data, weight, ef, bc);
    }

    //virtual void FillDataRequirements(TPZMaterialData &data);
    virtual void FillDataRequirements(TPZMaterialData &data) override;
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;

    virtual void FillBoundaryConditionDataRequirement(int type, TPZMaterialData &data) override;

    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)  override {
        PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
    }

    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    STATE Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);

    /// Transform a tensor to a Voigt notation
    static void ToVoigt(TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt);

    /// Transform a Voigt notation to a tensor
    static void FromVoigt(TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S);

    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template<class TVar>
    TVar InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T);

    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU);

    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<REAL> &GradU);

    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi);

    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);

public:

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &data, int var, TPZVec<STATE> &Solout) override;

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)  {
        TPZMaterial::SolutionDisc(data, dataleft, dataright, var, Solout) ;
    }
    /** 
     * @brief Computes the error due to the difference between the interpolated flux \n
     * and the flux computed based on the derivative of the solution
     */
protected:
    void Errors(TPZVec<REAL> &x, TPZVec<STATE> &u,
            TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
                TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values) override;
public:

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors) override;

    virtual int ClassId() const override;

    virtual void Read(TPZStream &buf, void *context) override;

    virtual void Write(TPZStream &buf, int withclassid) const override;



protected:
    /** @brief Forcing vector */
    TPZManVector<REAL, 3> fForce = TPZManVector<REAL, 3>(2, 0.);

    /** @brief Uses plain stress */
    int fPlaneStress = 1;

    /// dimension of the material
    int fDimension = 2;

    /** Elasticity function */
    TPZAutoPointer<TPZFunction<STATE> > fElasticity;
    
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
