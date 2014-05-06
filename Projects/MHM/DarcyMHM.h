/**
 * @file
 * @brief Contains the TPZMatDarcyMHM class.
 */

#ifndef MATDARCYMHM
#define MATDARCYMHM

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"


#ifdef _AUTODIFF
#include "fadType.h"
#endif


/**
 * @ingroup material
 * @brief Implements the MHM formulation to Darcy problem
 */
/**
 */
class TPZMatDarcyMHM : public TPZDiscontinuousGalerkin {
	
	protected :
	
	/** @brief Forcing function value */
	STATE fXf;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Coeficient which multiplies the Laplacian operator. */
	STATE fK;
    
    /** @brief Flag indicating whether the interface contributes positively or negatively */
    STATE fMultiplier;
    
			
public:
	
	
	TPZMatDarcyMHM(int nummat, int dim);
    
    TPZMatDarcyMHM(int matid) : TPZDiscontinuousGalerkin(matid), fXf(0.), fDim(-1), fK(1.), fMultiplier(1.)
    {
    
    }
	
	TPZMatDarcyMHM();
	
	TPZMatDarcyMHM(const TPZMatDarcyMHM &copy);
	

	virtual ~TPZMatDarcyMHM();
        
	TPZMatDarcyMHM &operator=(const TPZMatDarcyMHM &copy);
		
	virtual TPZMaterial * NewMaterial(){
		return new TPZMatDarcyMHM(*this);
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
    
    /**
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data)
    {
        data.SetAllRequirements(false);
        data.fNeedsNeighborCenter = true;
        data.fNeedsNormal = true;
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
    
	
	int Dimension() const { return fDim;}
	
	int NStateVariables();
	
    void SetDimension(int dim)
    {
        if(dim<0 || dim >3)
        {
            DebugStop();
        }
        fDim = dim;
    }
	
    void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }
    
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMatDarcyMHM"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
//	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) {
//		TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
//	}
    
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
//    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
//        TPZDiscontinuousGalerkin::Contribute(datavec,weight,ek,ef);
//    }
	
//    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
//        TPZDiscontinuousGalerkin::ContributeBC(datavec,weight,ek,ef,bc);
//    }

    
//	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
//		TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
//	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	

	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
//	virtual void ContributeBCInterface(TPZMaterialData &data,TPZMaterialData &dataleft,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
//	{
//		TPZDiscontinuousGalerkin::ContributeBCInterface(data,dataleft,weight,ef,bc);
//	}
	
//	virtual void ContributeInterface(TPZMaterialData &data,TPZMaterialData &dataleft,TPZMaterialData &dataright,REAL weight,
//									 TPZFMatrix<STATE> &ef)
//	{
//		TPZDiscontinuousGalerkin::ContributeInterface(data,dataleft,dataright,weight,ef);
//	}
    
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
	
protected:
//    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
//        TPZDiscontinuousGalerkin::Solution(datavec,var,Solout);
//    }
    
//    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) {
//        TPZDiscontinuousGalerkin::FillDataRequirements(datavec);
//    }
    
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
	
//	virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
	
	
	virtual int NEvalErrors() {return 6;}
	
	virtual int ClassId() const {
        return TPZMATDARCYMHM;
    }
	
	virtual void Write(TPZStream &buf, int withclassid);
	
	virtual void Read(TPZStream &buf, void *context);

};


#endif
