//
// Created by Giovane Avancini and Nathan Shauer on 18/12/23.
//

#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZMatInterfaceCombinedSpaces.h>
#include <TPZMatErrorCombinedSpaces.h>
#include <math.h>

#ifndef TPZELASTICITYTH
#define TPZELASTICITYTH

enum AnalysisType {EGeneral, EPlaneStrain, EPlaneStress};

class TPZElasticityTH : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>> {
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>>;
    
protected:
    /// @brief  Material dimension
    int fdimension;
    
    /// Approximation Space for Velocity
    int fSpace = 1;
    
    /// Young modulus
    REAL fyoung;

    /// Poisson coeficient
    REAL fpoisson;

    /// Lame first constant
    REAL flambda;

    /// Lame second constant
    REAL fmu;

    /// Bulk modulus
    REAL fbulk;

    /// Thickness in case of a plane analysis
    REAL fthickness;

    AnalysisType fAnalysisType;
    
    enum SpaceIndex {EUindex, EPindex, EVMindex, EPMindex};
    
    /// Big number for penalization method
    REAL fBigNumber = pow(10,std::numeric_limits<STATE>::max_digits10*2/3);
    
public:
    /// Empty Constructor
    TPZElasticityTH();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZElasticityTH(int matID, int dimension, REAL young_modulus, REAL poisson, AnalysisType analysisType = AnalysisType::EGeneral, REAL thickness = 1.0);
    
    /// Destructor
    ~TPZElasticityTH();
    
    /// returns the solution associated with the var index based on the finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>>&datavec, int var, TPZVec<STATE>& Solout) override;
    
    /// returns the number of variables associated with the variable indexed by var. Var is obtained by calling VariableIndex
    int NSolutionVariables(int var) const override;
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    /** returns the integrable dimension of the material */
    int Dimension() const override {return fdimension;}
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const override {return 1;}
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one domain integration point.
     * @param datavec[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     */
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors) override;
    
    int NEvalErrors() const override {return 8;}

    virtual void DeviatoricElasticityTensor(TPZFNMatrix<36,REAL>& D);

    virtual void ElasticityTensor(TPZFNMatrix<36,REAL>& D);

    virtual void StrainTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& epsilon);

    virtual void DeviatoricStressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma);

    virtual void StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma); //This function computes sigma directly from epsilon

    virtual void StressTensor(const TPZFNMatrix<10, STATE>& gradU, TPZFNMatrix<6,REAL>& sigma, REAL pressure); //This function computes sigma through its deviatoric and hidrostatic counterparts

    enum SolutionVars {ENone = -1, EPressure = 0, EDisplacement = 1, EForce = 2, EStress = 3, EStrain = 4, EVonMises = 5};
};

#endif
