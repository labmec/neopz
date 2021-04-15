/**
 * \file
 * @brief Contains the TPZElasticity3D class which implements a 3D isotropic elasticity material.
 */

#ifndef PZELAST3D
#define PZELAST3D

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <vector>


/**
 * @ingroup material
 * @brief This class implements a 3D isotropic elasticity material.
 * @since Aug 31, 2005.
 */
class TPZElasticity3D : public TPZMaterial {
    
    private:
    
    enum PLASTICPOSTPROC{ENonePlasticProc = -1, EVonMises = 0, EMohrCoulomb = 1};
    
    STATE VonMisesPlasticFunction( TPZFMatrix<STATE> & StressTensor) const;
    
    STATE MohrCoulombPlasticFunction( TPZFMatrix<STATE> & StressTensor) const;
    
    void Invariants( TPZFMatrix<STATE> & A, STATE & I1, STATE & I2, STATE & I3) const;
    
    void StressDecomposition( TPZFMatrix<STATE> & StressTensor, TPZFMatrix<STATE> & Deviator, STATE & p) const;
	
	public :
	
	enum SOLUTIONVARS {ENone = -1, EDisplacement = 0, EDisplacementX, EDisplacementY, EDisplacementZ,
		EPrincipalStress, EPrincipalStrain, EPrincipalDirection1, EPrincipalDirection2, EPrincipalDirection3,
		EVonMisesStress, EStress, EStrain, EStrain1, EStress1, ENormalStress, ENormalStrain, EStressX, EStressY, EStressZ,
        EI1, EI2, EI3, EPlasticFunction};
    

	
	/** 
	 * @brief Class constructor
	 * @param nummat - material ID.
	 * @param E - Young's modulus.
	 * @param poisson - poisson's ratio
	 * @param force - external forces
	 */
	TPZElasticity3D(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                    STATE preStressXX = 0., STATE preStressYY = 0., STATE preStressZZ = 0.);
    
	/** 
	 * @brief Constructor
	 * @param nummat - material ID.
	 */
	TPZElasticity3D(int nummat);
    
	/** @brief Default constructor */
	TPZElasticity3D();
    
	/** @brief Copy constructor */
	TPZElasticity3D(const TPZElasticity3D &cp);
	
	/** @brief Default destructor */
	virtual ~TPZElasticity3D();
	
	/** @brief Returns model dimension */
	int Dimension() const  override { return 3;}
	
	/** @brief Number of state variables */
	int NStateVariables() const override { return 3;}
	
    /** 
     * @brief Gets the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const override
    {
        return 2*(elPMaxOrder);
    }

	/** @brief Print material report */
	virtual void Print(std::ostream & out) override;
	
	/** 
	 * @brief Direction to post process stress and strain. \n
	 * Result of post processing is (Stress.Direction) or (Strain.Direction)
	 */
	void SetPostProcessingDirection(TPZVec<REAL> &Direction){
		if (Direction.NElements() != 3){
			PZError << __PRETTY_FUNCTION__ << " - ERROR!\n";
		}
		this->fPostProcessDirection.Resize(3);
		for(int i = 0; i < 3; i++) this->fPostProcessDirection[i] = Direction[i];
	}
	
	void SetVonMises(REAL fy){
      this->fFy = fy;
        this->fFrictionAngle = 0.;
        this->fCohesion = 0.;
        this->fPlasticPostProc = EVonMises;
    }
	
    void SetMohrCoulomb(REAL fc, REAL ft){
        this->fFy = 0.;
        this->fFrictionAngle = asin( (fc-ft)/(fc+ft) );
        this->fCohesion = fc*ft/(fc-ft) * tan(fFrictionAngle);
        this->fPlasticPostProc = EMohrCoulomb;
    }
    
	/** @brief Material name */
	virtual std::string Name() override { return "TPZElasticity3D"; }
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override;
    
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
    
    virtual void ContributeVecShape(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	/** @brief Implements Dirichlet and Neumann boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
    
    virtual void ContributeVecShapeBC(TPZMaterialData & data, REAL weight,
                                      TPZFMatrix<STATE> & ek, TPZFMatrix<STATE> & ef,TPZBndCond &bc);
    
	/** @brief Implements Dirichlet and Neumann boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	/** @brief Returns index of post-processing variable */
	virtual int VariableIndex(const std::string &name) override;
	
	/** @brief Number of data of variable var */
	virtual int NSolutionVariables(int var) override;
protected:
	/** @brief Post-processing method. Based on solution Sol and its derivatives DSol, it computes the post-processed variable var */
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** @brief Evaluate error between approximate (FEM) and exact solutions */
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
						TPZFMatrix<REAL> &axes,
						TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	/** @brief Returns the number of norm errors: 3 (Semi H1, L2 and H1) */
	virtual int NEvalErrors() override {return 3;}
	
	/** @brief Fill material data parameter with necessary requirements for the Contribute method. */
	void FillDataRequirements(TPZMaterialData &data) override;
	
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override
    {
        // default is no specific data requirements
        if(type == 50)
        {
            data.fNeedsSol = true;
        }
        data.fNeedsNormal = true;
    }

	void SetMaterialDataHook(REAL Ela, REAL poisson)
	{
		fE = Ela;
		fPoisson = poisson;
        SetC();
	}
	
	void SetMaterialDataLame(REAL lambda, REAL mu)
	{
		fE = mu*((REAL(3))*lambda+((REAL(2))*mu))/(lambda + mu);
		fPoisson = lambda/((REAL(2))*(lambda + mu));
        SetC();
	}
	
	void SetC()
	{
		C1 = fE / (2.+ 2.*fPoisson);
		C2 = fE * fPoisson / (-1. + fPoisson + 2.*fPoisson*fPoisson);
		C3 = fE * (fPoisson - 1.) / (-1. + fPoisson +2. * fPoisson * fPoisson);
	}
	
	void SetForce(TPZVec <STATE> force)
	{
		for(int i = 0; i < 3; i++) this->fForce[i] = force[i];
	}
    
    STATE GetE()
    {
        return fE;
    }
    
    STATE GetPoisson()
    {
        return fPoisson;
    }
    
    STATE GetLambda()
    {
        STATE lambda = (fPoisson*fE)/((1.+fPoisson)*(1.-2.*fPoisson));
        return lambda;
    }
    
    STATE GetMU()
    {
        STATE mu = fE/(2.*(1.+fPoisson));
        return mu;
    }
    
    STATE GetPrestress(int index)
    {
        return this->fPreStress[index];
    }
	
	protected :
	
	/** @brief Young's modulus */
	STATE fE;
	
	/** @brief Poisson's ratio */
	STATE fPoisson;
	
	REAL C1; /**< \f$ C1 = E / (2.+ 2.*nu) \f$ */
	REAL C2; /**< \f$ C2 = E * nu / (-1. + nu + 2.*nu*nu) \f$ */
	REAL C3; /**< \f$ C3 = E * (nu - 1.) / (-1. + nu +2. * nu * nu) \f$ */
	/** @brief External forces */
	TPZManVector<STATE,3> fForce;
	
	/** @brief Direction to compute stress and strain */
	TPZManVector<REAL,3> fPostProcessDirection;
	
	/** @brief Yeilding stress for von mises post processing */
	REAL fFy;
    
    /** @brief Mohr-Coulomb parameters */
    REAL fFrictionAngle, fCohesion;
    
    /** @brief Plastic model for post-processing */
    PLASTICPOSTPROC fPlasticPostProc;
    
    TPZManVector<REAL> fPreStress;//just prestress XX, YY and ZZ
	
public:
	virtual void ComputeStressVector(TPZFMatrix<STATE> &Stress, TPZFMatrix<STATE> &DSol) const;
	void ComputeStrainVector(TPZFMatrix<STATE> &Strain, TPZFMatrix<STATE> &DSol) const;
	virtual void ComputeStressTensor(TPZFMatrix<STATE> &Stress, TPZMaterialData &data) const;
    void ComputeStressTensor(TPZFMatrix<STATE> &Stress, TPZFMatrix<STATE> &DSol) const;
	void ComputeStrainTensor(TPZFMatrix<STATE> &Strain, TPZFMatrix<STATE> &DSol) const;
	void ApplyDirection(TPZFMatrix<STATE> &StrVec, TPZVec<STATE> &Out) const;
	void PrincipalDirection(TPZFMatrix<STATE> &DSol, TPZVec< STATE > &Solout, int direction) const;
	
public:
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context) override;
	public:
virtual int ClassId() const override;

	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial()  override { return new TPZElasticity3D(*this);}
	
	static STATE gTolerance;
	
};

#endif
