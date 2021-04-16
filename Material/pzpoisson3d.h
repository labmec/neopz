/**
 * @file
 * @brief Contains the TPZMatPoisson3d class.
 */

#ifndef MATPOISSON3DH
#define MATPOISSON3DH

#include <iostream>

#include "pzfmatrix.h"
#include "TPZMaterial.h"

#include "fadType.h"


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 * \f$ -fK Laplac(u) + fC * div(fConvDir*u) = - fXf  \f$
 */
class TPZMatPoisson3d : public TPZMaterial {
	
	protected :
	
	/** @brief Forcing function value */
	STATE fXf;
	
	/** @brief Problem dimension */
	int fDim;

    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
    TPZFNMatrix<9,STATE> fTensorPerm;

    /** @brief inverse of the permeability tensor.*/
    TPZFNMatrix<9,STATE> fInvPerm;
		
	/** @brief Variable which multiplies the convection term of the equation */
	REAL fC;
	
	/** @brief Direction of the convection operator */
	REAL fConvDir[3];
	
	/** @brief Symmetry coefficient of elliptic term */
	/** 
	 * Symmetrical formulation - Global element method - has coefficient = -1. \n
	 * Non-symmetrical formulation - Baumann's formulation - has coefficient = +1.
	 */
	REAL fSymmetry;
	
	/** @brief Multiplication value for the streamline diffusion term */
	STATE fSD;
	
	/** @brief Enumerate for penalty term definitions */
	enum EPenaltyType {ENoPenalty = 0, EFluxPenalty = 1, ESolutionPenalty, EBoth};
	
	/** @brief Penalty term definition */
	EPenaltyType fPenaltyType;
   
	
public:
	
	/** @brief Constant multiplyer of penalty term, when required is set. */
	REAL fPenaltyConstant;
	
	/** @brief Defines no penalty terms in ContributeInterface */
	void SetNoPenalty(){ this->fPenaltyType = ENoPenalty;}
	
	/** @brief Defines flux penalty terms in ContributeInterface */
	void SetFluxPenalty(){ this->fPenaltyType = EFluxPenalty; }
	
	/** @brief Defines solution penalty terms in ContributeInterface */
	void SetSolutionPenalty(){ this->fPenaltyType = ESolutionPenalty; }
	
	/** @brief Defines solution and flux penalty terms in ContributeInterface */
	void SetBothPenalty(){ this->fPenaltyType = EBoth; }
	
	/** @brief Using in InterfaceErrors */
	static STATE gAlfa;
    /** to solve neumann problem for error estimation*/
     bool fNeumann;
	
	TPZMatPoisson3d(int nummat, int dim);
    
    TPZMatPoisson3d(int matid) : TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
    TPZMaterial(matid), fXf(0.), fDim(-1), fC(0.), fSymmetry(0.), fSD(0.)
        ,fPenaltyType(ENoPenalty)
    {
    
    }
	
	TPZMatPoisson3d();
	
	TPZMatPoisson3d(const TPZMatPoisson3d &copy);
	
	virtual ~TPZMatPoisson3d();
	
	TPZMatPoisson3d &operator=(const TPZMatPoisson3d &copy);
	
	/** @brief Set material elliptic term as the global element method, i.e. the symmetrical formulation */
	void SetSymmetric(){
		this->fSymmetry = -1.0;
	}
    
    void SetNeumannProblem(){
        fNeumann=true;
    }
    
	
	/** @brief Set material elliptic term as the Baumann's formulation, i.e. the non-symmetrical formulation */
	void SetNonSymmetric() {
		this->fSymmetry = +1.0;
	}
	
	bool IsSymetric(){
		if (fSymmetry == -1.0) return true;
		if (fSymmetry == +1.0) return false;
		PZError << __PRETTY_FUNCTION__ << "\n Comparacao de numeros reais da errado\n";
		return false;
	}
	
	virtual TPZMaterial * NewMaterial() override {
		return new TPZMatPoisson3d(*this);
	}
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** 
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
    }
	    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override
    {
        data.SetAllRequirements(false);
        if (type == 50) {
            data.fNeedsSol = true;
        }
        if (type == 3) {
            data.fNeedsNormal = true;
        }
    }
    
	
	int Dimension() const override{ return fDim;}
	
    /** @brief Returns the number of state variables associated with the material */
	virtual int NStateVariables() const override
    {
        return 1;
    }
	
	virtual void SetParameters(STATE diff, REAL conv, TPZVec<REAL> &convdir);

    virtual void SetParameters(TPZFNMatrix<9,STATE> tensorPerm,TPZFNMatrix<9,STATE> invPerm, REAL conv, TPZVec<REAL> &convdir);
	
	void GetParameters(TPZFNMatrix<9,STATE> &tensorPerm, TPZFNMatrix<9,STATE> &invPerm, REAL &conv, TPZVec<REAL> &convdir);

	void GetParameters(STATE &diff, REAL &conv, TPZVec<REAL> &convdir);

    void SetPermeability(STATE perm);

    void SetPermeabilityTensor(const TPZFNMatrix<9,STATE> K, const TPZFNMatrix<9,STATE> invK);

    void GetPermeability(TPZFNMatrix<9,STATE> &K);

    void GetInvPermeability(TPZFNMatrix<9,STATE> &invK);

    void GetPermeabilities(TPZVec<REAL> &x, TPZFNMatrix<9,STATE> &PermTensor, TPZFNMatrix<9,STATE> &InvPermTensor);
    
    void SetDimension(int dim)
    {
        if(dim<0 || dim >3)
        {
            DebugStop();
        }
        fDim = dim;
    }
	
	void SetInternalFlux(STATE flux)
	{
		fXf = flux;
	}
	
	void SetSD(STATE sd)
	{
		fSD = sd;
	}
	
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override { return "TPZMatPoisson3d"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) override {
		TPZMaterial::Contribute(data,weight,ef);
	}
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
        TPZMaterial::Contribute(datavec,weight,ek,ef);
    }
	
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        TPZMaterial::ContributeBC(datavec,weight,ek,ef,bc);
    }

   virtual void LocalNeumanContribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	/** @brief Computes contribution to the energy at an integration point */
	void ContributeEnergy(TPZVec<REAL> &x,
						  TPZVec<FADFADREAL> &sol,
						  TPZVec<FADFADREAL> &dsol,
						  FADFADREAL &U,
						  REAL weight);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	

	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBCInterface(TPZMaterialData &data,TPZMaterialData &dataleft,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data,TPZMaterialData &dataleft,TPZMaterialData &dataright,REAL weight,
									 TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::ContributeInterface(data,dataleft,dataright,weight,ef);
	}

	virtual void ContributeBCEnergy(TPZVec<REAL> &x,TPZVec<FADFADREAL> &sol, FADFADREAL &U,
									REAL weight, TPZBndCond &bc);

	/** @} */
    
	STATE AVGK();

	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override {
        TPZMaterial::Solution(datavec,var,Solout);
    }
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec)  override {
        TPZMaterial::FillDataRequirements(datavec);
    }
    
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;

	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
                        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
    
	
	
	virtual int NEvalErrors()  override {return 6;}
	
	/**
	 * @brief Compute square of residual of the differential equation at one integration point.
	 * @param X is the point coordinate (x,y,z)
	 * @param sol is the solution vector
	 * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
	 */   
	virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) override;
	
	
	void InterfaceErrors(TPZVec<REAL> &/*x*/,
						 TPZVec<STATE> &leftu, TPZFMatrix<STATE> &leftdudx, /* TPZFMatrix<REAL> &leftaxes,*/
						 TPZVec<STATE> &rightu, TPZFMatrix<STATE> &rightdudx, /* TPZFMatrix<REAL> &rightaxes,*/
						 TPZVec<STATE> &/*flux*/,
						 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values,
						 TPZVec<STATE> normal, STATE elsize);
	
	/** 
	 * @brief Computes interface jump from element to Dirichlet boundary condition
	 * @return Returns sol-u_dirichlet
	 * @since Mar 08, 2006
	 */
	virtual void BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump) override;
	
	virtual int IsInterfaceConservative() override { return 1;}
	
        public:
virtual int ClassId() const override;

	
	virtual void Write(TPZStream &buf, int withclassid) const override;
	
	virtual void Read(TPZStream &buf, void *context) override;

};

#endif
