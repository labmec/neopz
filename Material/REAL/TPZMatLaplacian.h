/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#ifndef MATLAPLACDH
#define MATLAPLACDH

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief \f$ -fK Laplac(u) = fXf  \f$
 */
/**
 * \f$ -fK Laplac(u) = fXf  \f$
 */
class TPZMatLaplacian : public TPZDiscontinuousGalerkin {
	
	protected :
	
	/** @brief Forcing function value */
	STATE fXf;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Coeficient which multiplies the Laplacian operator. */
	STATE fK;
		
	/** @brief Symmetry coefficient of elliptic term */
	/** 
	 * Symmetrical formulation - Global element method - has coefficient = -1. \n
	 * Non-symmetrical formulation - Baumann's formulation - has coefficient = +1.
	 */
	REAL fSymmetry;
	
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

	TPZMatLaplacian(int nummat, int dim);

  TPZMatLaplacian(int matid) : TPZDiscontinuousGalerkin(matid), fXf(0.), fDim(1), fK(1.),
     fSymmetry(0.), fPenaltyType(ENoPenalty)
  {

  }

	TPZMatLaplacian();

	TPZMatLaplacian(const TPZMatLaplacian &copy);

	virtual ~TPZMatLaplacian();

	TPZMatLaplacian &operator=(const TPZMatLaplacian &copy);

	/** @brief Set material elliptic term as the global element method, i.e. the symmetrical formulation */
	void SetSymmetric(){
		this->fSymmetry = -1.0;
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

	virtual TPZMaterial * NewMaterial(){
		return new TPZMatLaplacian(*this);
	}

    /**
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/**
	 * Contribute method. Here, in base class, all requirements are considered as necessary.
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
    }

    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        if (type == 50) {
            data.fNeedsSol = true;
        }
        if (type == 3) {
            data.fNeedsNormal = true;
        }
    }


	virtual int Dimension() const { return fDim;}

	int NStateVariables();

	void SetParameters(STATE diff, STATE f);

	void GetParameters(STATE &diff, STATE &f);

  void SetDimension(int dim)
  {
      if(dim<0 || dim >3)
      {
          DebugStop();
      }
      fDim = dim;
  }


	virtual void Print(std::ostream & out);

	virtual std::string Name() { return "TPZMatLaplacian"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */

    virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
        TPZDiscontinuousGalerkin::Contribute(datavec,weight,ek,ef);
    }


	virtual void ContributeBCHDiv(TPZMaterialData &data,REAL weight,
								  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual void ContributeHDiv(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	/** @} */

	virtual int VariableIndex(const std::string &name);

	virtual int NSolutionVariables(int var);

	virtual int NFluxes(){ return 3;}

protected:

	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);

  public:

	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
        DebugStop();
    }

	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);

	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	void ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);


	virtual int NEvalErrors() {return 3;}

	/**
	 * @brief Compute square of residual of the differential equation at one integration point.
	 * @param X is the point coordinate (x,y,z)
	 * @param sol is the solution vector
	 * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
	 */
	virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);


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
	virtual void BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump);

	virtual int IsInterfaceConservative(){ return 1;}

    virtual int ClassId() const {
        return TPZMatLaplacianID;
    }

	virtual void Write(TPZStream &buf, int withclassid);

	virtual void Read(TPZStream &buf, void *context);

};

#endif

