/**
 * @file
 * @brief Contains the TPZDiscontinuousGalerkin class which implements the interface for discontinuous Galerkin formulation.
 */

#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZMaterialData;

/**
 * @ingroup material
 * @brief Defines the interface which material objects need to implement for discontinuous Galerkin formulations
 */
class TPZDiscontinuousGalerkin : public TPZMaterial {
	
	public :
	/** @brief Simple constructor */
	TPZDiscontinuousGalerkin();
	/** @brief Constructor with the index of the material object within the vector */
	TPZDiscontinuousGalerkin(int nummat);
	
	/** @brief Copy constructor */
	TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy);
	/** @brief Destructor */
	virtual ~TPZDiscontinuousGalerkin();
	
	virtual std::string Name();
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the ContributeInterface method.
     * @since April 10, 2007
	 */
	/** 
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data);
	
    /**
     * @{
     * @name Contribute methods
     * @}
     */
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) = 0;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
        TPZMaterial::Contribute(datavec,weight,ek,ef);
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
    {
        TPZMaterial::Contribute(datavec,weight,ef);
    }


    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) = 0;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
        TPZMaterial::ContributeBC(datavec,weight,ek,ef,bc);
    }

	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) {
		TPZMaterial::Contribute(data,weight,ef);
	}
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        TPZMaterial::ContributeBC(datavec,weight,ef,bc);
    }

	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
                                     REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

	
	/**
	 * @brief It computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef);
	
    
	/**
	 * @brief Computes a contribution to residual vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @since June 5, 2012
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef);
	
    
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point 
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) = 0;
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point to multiphysics simulation
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since February 21, 2013
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
	
	/**
	 * @brief It computes a contribution to residual vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<STATE> &Solout)
    {
        std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
    }
	
    /** @} */
	
	/**
	 * @brief Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
	 * @since Feb 05, 2004
	 */
	/** 
	 * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
	 * when using multigrid pre-conditioner. \n Conservative interfaces into agglomerate elements do not
	 * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.\n
	 * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
	 * does not ruin the solution.
	 */
	virtual int IsInterfaceConservative();
	
	/** 
	 * @brief Computes interface jump = leftu - rightu
	 * @since Feb 14, 2006
	 */
	virtual void InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump);
	
	
	/** 
	 * @brief Computes interface jump from element to Dirichlet boundary condition.
	 * It has to reimplemented
	 * @since Mar 08, 2006
	 */
	virtual void BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump);
	
	
	virtual int NStateVariables() = 0;
	
	
	virtual void ContributeInterfaceErrors(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
										   REAL weight,
										   TPZVec<STATE> &nkL,
										   TPZVec<STATE> &nkR,
										   int &errorid) {
		PZError << "Method not implemented\n";
	}
	
	virtual void ContributeInterfaceBCErrors(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZVec<STATE> &nk,
											 TPZBndCond &bc,
											 int &errorid) {
		PZError << "Method not implemented\n";
	}
	
    virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                        TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
        TPZMaterial::Errors(x,u, dudx, axes, flux, u_exact,du_exact,values);

    }

    virtual void Errors(TPZMaterialData &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);

    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors);

    /** @{
     * @name Save and Load methods
     */

	/** @brief Unique identifier for serialization purposes */
	public:
virtual int ClassId() const;

	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
    /**
     * @}
     */
};


#endif
