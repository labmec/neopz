/**
 * \file
 * @brief Contains the TPZDiscontinuousGalerkin class which implements the interface for discontinuous Galerkin formulation.
 */
// $Id: pzdiscgal.h,v 1.20 2009-11-04 14:04:43 fortiago Exp $
#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZMaterialData;

/**
 * @ingroup material
 * @brief Defines the interface which material objects need to implement for discontinuous Galerkin formulations
 */
class TPZDiscontinuousGalerkin  : public TPZMaterial {
	
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
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) = 0;
	
	
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
    void SolutionDisc(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, int var, TPZVec<REAL> &Solout)
    {
        std::cout << __PRETTY_FUNCTION__ << " should never be called\n";
    }
	
	
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
										   TPZVec<REAL> &nkL,
										   TPZVec<REAL> &nkR,
										   int &errorid){
		PZError << "Method not implemented\n";
	}
	
	virtual void ContributeInterfaceBCErrors(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZVec<REAL> &nk,
											 TPZBndCond &bc,
											 int &errorid){
		PZError << "Method not implemented\n";
	}
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const;
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
};

#endif
