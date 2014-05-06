/**
 * \file
 * @brief Contains the TPZElasticity3D class which implements a 3D isotropic elasticity material.
 */

#ifndef PZELAST3D
#define PZELAST3D

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <vector>

/**
 * @ingroup material
 * @brief This class implements a 3D isotropic elasticity material.
 * @since Aug 31, 2005.
 */
class TPZElasticity3D : public TPZMaterial {
	
	public :
	
	enum SOLUTIONVARS {ENone = -1, EDisplacement = 0, EDisplacementX, EDisplacementY, EDisplacementZ,
		EPrincipalStress, EPrincipalStrain, EPrincipalDirection1, EPrincipalDirection2, EPrincipalDirection3,
		EVonMisesStress, EStress, EStrain, EStrain1, EStress1, ENormalStress, ENormalStrain, EStressX, EStressY, EStressZ};
	
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
	int Dimension() const { return 3;}
	
	/** @brief Number of state variables */
	int NStateVariables(){ return 3;}
	
    /** 
     * @brief Gets the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const
    {
        return 2*(elPMaxOrder);
    }

	/** @brief Print material report */
	virtual void Print(std::ostream & out);
	
	/** 
	 * @brief Direction to post process stress and strain. \n
	 * Result of post processing is (Stress.Direction) or (Strain.Direction)
	 */
	void SetPostProcessingDirection(TPZVec<STATE> &Direction){
		if (Direction.NElements() != 3){
			PZError << __PRETTY_FUNCTION__ << " - ERROR!\n";
		}
		this->fPostProcessDirection.Resize(3);
		for(int i = 0; i < 3; i++) this->fPostProcessDirection[i] = Direction[i];
	}
	
	void SetYieldingStress(STATE fy){ this->fFy = fy; }
	
	/** @brief Material name */
	virtual std::string Name() { return "TPZElasticity3D"; }
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);
    
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
    
    virtual void ContributeVecShape(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	/** @brief Implements Dirichlet and Neumann boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc);
    
    virtual void ContributeVecShapeBC(TPZMaterialData & data, REAL weight,
                                      TPZFMatrix<STATE> & ek, TPZFMatrix<STATE> & ef,TPZBndCond &bc);
    
	/** @brief Implements Dirichlet and Neumann boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	/** @brief Returns index of post-processing variable */
	virtual int VariableIndex(const std::string &name);
	
	/** @brief Number of data of variable var */
	virtual int NSolutionVariables(int var);
protected:
	/** @brief Post-processing method. Based on solution Sol and its derivatives DSol, it computes the post-processed variable var */
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** 
	 * @brief Return the number of components which form the flux function
	 * @note Method not implemented.
	 */
	virtual int NFluxes() {
		PZError << "\nTPZElasticity3D::NFluxes() - Method not implemented\n";
		return 0;
	}
	
	/**
	 * @brief Compute the value of the flux function to be used by ZZ error estimator.
	 * @note Method not implemented.
	 */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux){
		PZError << "\nTPZElasticity3D::Flux - Method not implemented\n";
	}
	
	/** @brief Evaluate error between approximate (FEM) and exact solutions */
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u, TPZFMatrix<STATE> &dudx,
						TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
						TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	/** @brief Returns the number of norm errors: 3 (Semi H1, L2 and H1) */
	virtual int NEvalErrors() {return 3;}
	
	/** @brief Fill material data parameter with necessary requirements for the Contribute method. */
	void FillDataRequirements(TPZMaterialData &data);
	
	void SetMaterialDataHook(STATE Ela, STATE poisson)
	{
		fE = Ela;
		fPoisson = poisson;
        SetC();
	}
	
	void SetMaterialDataLame(STATE lambda, STATE mu)
	{
		fE = mu*((REAL(3))*lambda+((REAL(2))*mu))/(lambda + mu);
		fPoisson = lambda/((REAL(2))*(lambda + mu));
        SetC();
	}
	
	void SetC()
	{
#ifndef CODE1
		C1 = fE / (2.+ 2.*fPoisson);
		C2 = fE * fPoisson / (-1. + fPoisson + 2.*fPoisson*fPoisson);
		C3 = fE * (fPoisson - 1.) / (-1. + fPoisson +2. * fPoisson * fPoisson);
#endif
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
	
#ifndef CODE1
	STATE C1; /**< \f$ C1 = E / (2.+ 2.*nu) \f$ */
	STATE C2; /**< \f$ C2 = E * nu / (-1. + nu + 2.*nu*nu) \f$ */
	STATE C3; /**< \f$ C3 = E * (nu - 1.) / (-1. + nu +2. * nu * nu) \f$ */
#endif
	/** @brief External forces */
	TPZManVector<STATE,3> fForce;
	
	/** @brief Direction to compute stress and strain */
	TPZManVector<STATE,3> fPostProcessDirection;
	
	/** @brief Yeilding stress */
	STATE fFy;
    
    TPZManVector<STATE> fPreStress;//just prestress XX, YY and ZZ
	
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
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	virtual int ClassId() const;
	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial * NewMaterial() { return new TPZElasticity3D(*this);}
	
	static STATE gTolerance;
	
};

#endif
