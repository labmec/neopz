/**
 * \file
 * @brief Contains the TPZElasticity3DGD class which implements a 3D isotropic elasticity material.
 */

#ifndef PZELAST3D
#define PZELAST3D

#include <iostream>
#include "pzdiscgal.h"
#ifdef _AUTODIFF
#include "fadType.h"
#endif
#include "pzfmatrix.h"
#include "pzvec.h"
#include "tpzintpoints.h"
#include "tpzautopointer.h"
#include <vector>
//#include "TPZPlasticModel.h"

/**
 * @ingroup material
 * @brief This class implements a 3D isotropic elasticity material.
 * @since Aug 31, 2005.
 */
class TPZElasticity3DGD : public TPZDiscontinuousGalerkin {
	
	public :
	
	enum SOLUTIONVARS {ENone = -1, EDisplacement = 0, EDisplacementX, EDisplacementY, EDisplacementZ,
		EPrincipalStress, EPrincipalStrain, EPrincipalDirection1, EPrincipalDirection2, EPrincipalDirection3,
		EVonMisesStress, EStress, EStrain, EStrain1, EStress1, ENormalStress, ENormalStrain, EStressX, EStressY, EStressZ,
    EI1, EI2, EI3, EStressVector,EStressTensor,EDeviatoricStressTensor,EDeviatoricStrainTensor,EI1Strain,EStrainTensor,EStressFlac,
    EMatId, EYielding, EDamage, EStrainZ, EPlasticStrainZ, EPlasticFunctionVal, EElementSize,
    EGeoElIndex, EPrincipalPlasticStrain};
	
	/** 
	 * @brief Class constructor
	 * @param nummat - material ID.
	 * @param E - Young's modulus.
	 * @param poisson - poisson's ratio
	 * @param force - external forces
	 */
	TPZElasticity3DGD(int nummat, REAL E, REAL poisson, TPZVec<REAL> &force,
                    REAL preStressXX = 0., REAL preStressYY = 0., REAL preStressZZ = 0.);

	/** 
	 * @brief Constructor
	 * @param nummat - material ID.
	 */
	TPZElasticity3DGD(int nummat);
    
	/** @brief Default constructor */
	TPZElasticity3DGD();
    
	/** @brief Copy constructor */
	TPZElasticity3DGD(const TPZElasticity3DGD &cp);
	
	/** @brief Default destructor */
	virtual ~TPZElasticity3DGD();
	
	/** @brief Returns model dimension */
	virtual int Dimension() const { return 3;}

  virtual bool UsesBigNumberForBC(int bcType) const;
	
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
	void SetPostProcessingDirection(TPZVec<REAL> &Direction){
		if (Direction.NElements() != 3){
			PZError << __PRETTY_FUNCTION__ << " - ERROR!\n";
		}
		this->fPostProcessDirection.Resize(3);
		for(int i = 0; i < 3; i++) this->fPostProcessDirection[i] = Direction[i];
	}
	
	void SetYieldingStress(REAL fy){ this->fFy = fy; }
	
	/** @brief Material name */
	virtual std::string Name() { return "TPZElasticity3DGD"; }
	
	/** Contribute to stiff matrix and load vector.
	 *  See base class to more informations.
	 */
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);

  private:

	void ContributeOriginal(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);

#ifdef _AUTODIFF
	void ContributeDifFinita(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef);
#endif

  public:

	/** Contribute to stiff matrix and load vector.
	 *  See base class to more informations.
	 */
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ef);

	/** @brief Implements Dirichlet and Neumann boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc);

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

  ///Modelo elastoplastico
//  TPZAutoPointer< TPZPlasticModel > fPlasticModel;

  ///Modelo de dano continuo.
//  TPZAutoPointer< TPZDamageModel > fDamageModel;

  ///Para fraturas coesivas. Usado em ContributeInterface
//  TPZAutoPointer< TPZCohesiveDamageModel > fCohesiveDamageModel;

  void Vector2Tensor(const TPZVec<REAL> &Vec, TPZFMatrix<STATE> &Tensor) const;

  void Tensor2Vector(const TPZFMatrix<STATE> &Tensor, TPZVec<REAL> &Vec) const;

//  REAL GetDamage(/*const TPZVec<REAL> &Strain,*/ int elindex, int elIntPoint) const;

public:

	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

	/** 
	 * @brief Return the number of components which form the flux function
	 * @note Method not implemented.
	 */
	virtual int NFluxes() {
		PZError << "\nTPZElasticity3DGD::NFluxes() - Method not implemented\n";
		return 0;
	}
	
	/**
	 * @brief Compute the value of the flux function to be used by ZZ error estimator.
	 * @note Method not implemented.
	 */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<STATE> &axes, TPZVec<REAL> &flux){
		PZError << "\nTPZElasticity3DGD::Flux - Method not implemented\n";
	}

	/** @brief Evaluate error between approximate (FEM) and exact solutions */
	virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix<STATE> &dudx,
		TPZFMatrix<STATE> &axes, TPZVec<REAL> &flux,
		TPZVec<REAL> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &values);
	
	/** @brief Returns the number of norm errors: 3 (Semi H1, L2 and H1) */
	virtual int NEvalErrors() {return 3;}
	
	/** @brief Fill material data parameter with necessary requirements for the Contribute method. */
	void FillDataRequirements(TPZMaterialData &data);
	
	void SetMaterialDataHook(REAL Ela, REAL poisson)
	{
		fE = Ela;
		fPoisson = poisson;
	}
	
	void SetMaterialDataLame(REAL lambda, REAL mu)
	{
		fE = mu*((REAL(3))*lambda+((REAL(2))*mu))/(lambda + mu);
		fPoisson = lambda/((REAL(2))*(lambda + mu));
	}

	
	void SetForce(TPZVec <REAL> force)
	{
		for(int i = 0; i < 3; i++) this->fForce[i] = force[i];
	}
    
    REAL GetE()
    {
        return fE;
    }
    
    REAL GetPoisson()
    {
        return fPoisson;
    }
    
    REAL GetLambda()
    {
        REAL lambda = (fPoisson*fE)/((1.+fPoisson)*(1.-2.*fPoisson));
        return lambda;
    }
    
    REAL GetMU()
    {
        REAL mu = fE/(2.*(1.+fPoisson));
        return mu;
    }
    
    REAL GetPrestress(int index)
    {
        return this->fPreStress[index];
    }
	
	protected :
	
	/** @brief Young's modulus */
	REAL fE;
	
	/** @brief Poisson's ratio */
	REAL fPoisson;

	/** @brief External forces */
	TPZManVector<REAL,3> fForce;
	
	/** @brief Direction to compute stress and strain */
	TPZManVector<REAL,3> fPostProcessDirection;
	
	/** @brief Yeilding stress */
	REAL fFy;
    
    TPZManVector<REAL> fPreStress;//just prestress XX, YY and ZZ

	virtual void ComputeStressVector(TPZVec<REAL> &Stress, bool &yielding, REAL &damage, const TPZFMatrix<STATE> &DSol, int elindex, int elIntPoint, const REAL elementCharacteristicSize) const;
  void ComputeStressVector(TPZVec<REAL> &Stress, bool &yielding, REAL &damage, const TPZVec<REAL> &Strain, int elindex, int elIntPoint, const REAL elementCharacteristicSize) const;
  void ComputeStrainVector(TPZVec<REAL> &Strain, const TPZFMatrix<STATE> &DSol) const;
  void ComputeStressTensor(TPZFMatrix<STATE> &Stress, bool &yielding, REAL &damage, const TPZFMatrix<STATE> &DSol, int elindex, int elIntPoint, const REAL elementCharacteristicSize) const;
  void ComputeStrainTensor(TPZFMatrix<STATE> &Strain, const TPZFMatrix<STATE> &DSol) const;
	void ApplyDirection(const TPZVec<REAL> &StrVec, TPZVec<REAL> &Out) const;
	void PrincipalDirection(const TPZFMatrix<STATE> &DSol, TPZVec< REAL > &Solout, int direction) const;


	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
		public:
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

		protected :
	public:
virtual int ClassId() const;
public:

	/** @brief Creates a new material from the current object   ??*/
	virtual TPZMaterial* NewMaterial() { return new TPZElasticity3DGD(*this);}

  /** Set the integration rule order based on the element
   *   p order of interpolation, its dimension and the characteristics
   *   of the material
   */
  virtual void SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
                                  int elPMaxOrder,
                                  int elDimension);

  /** Set the integration rule order based on the element
   *   p order of interpolation, its dimension and the characteristics
   *   of the material for a boundary condition element
   */
  virtual void SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
                                  int elPMaxOrder,
                                  int elDimension,
                                  TPZBndCond &bc);
	
	static REAL gTolerance;

  /** calculando jacobiana por diferenca finita
    * jac = D[Stress,Strain]
    */
	void EstimateJacobian(TPZFMatrix<STATE> &dsol, int elindex, int elIntPoint, const REAL elementCharacteristicSize,
                        TPZVec< TPZVec<REAL> > &jac);

  // *****    TPZDiscontinuousGalerkin  **********//

    /** Fill material data parameter with necessary requirements for the
    * ContributeInterface method. Here, in base class, all requirements are considered
    * as necessary. Each derived class may optimize performance by selecting
    * only the necessary data.
    * @since April 10, 2007
    */
  virtual void FillDataRequirementsInterface(TPZMaterialData &facedata);

  /**
   * It computes a contribution to stiffness matrix and load vector at one integration point
   * @param data [in]
   * @param weight [in]
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @since April 16, 2007
   */
  virtual void ContributeInterface(TPZMaterialData &facedata,
                                   TPZMaterialData &leftdata, TPZMaterialData &rightdata,
								   REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

  /**
   * It computes a contribution to stiffness matrix and load vector at one BC integration point
   * @param data [in]
   * @param weight [in]
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition object
   * @since April 16, 2007
   */
  virtual void ContributeBCInterface(TPZMaterialData &facedata,
                                     TPZMaterialData &leftdata,
									 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                     TPZBndCond &bc);

  protected:

  //Returns gamma coeficient for ContributeInterface method, which implements a solution penalty flux
  void InterfaceGamma(const int elindex, const int elIntPoint,
                      REAL h, REAL normalOpening,
                      REAL &gamma, REAL &maxTensionStress, bool &limitedTensionStress) const;

};

#endif
