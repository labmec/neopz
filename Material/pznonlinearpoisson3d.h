/**
 * \file
 * @brief Contains the TPZNonLinearPoisson3d class.
 */

//$Id: pznonlinearpoisson3d.h,v 1.7 2009-09-01 19:44:48 phil Exp $

#ifndef MATNLPOISSON3DH
#define MATNLPOISSON3DH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZNonLinearPoisson3d : public TPZMatPoisson3dReferred {
	
protected:
	
	/** @brief Definitions of stabilization terms */
	enum EStabilizationType {ENoStabilization = 0, ESUPG = 1, EGradient = 2};
	
	/** @brief Stabilization term definition */
	int fStabilizationType;
	
public:
	
	TPZNonLinearPoisson3d(int nummat, int dim);
	
	TPZNonLinearPoisson3d(const TPZNonLinearPoisson3d &cp);
	
	virtual ~TPZNonLinearPoisson3d();
	
	bool IsReferred(){ return this->fIsReferred;}
	
	void SetReferred(bool Is){ this->fIsReferred = Is; }
	
	/** @brief Define SUPG stabilization term. */
	void SetSUPGStab(REAL sd = 1.0);
	
	/** @brief Define gradient stabilization term. */
	void SetGradientStab(REAL sd = 1.0);
	
	/** @brief Define no stabilization term. */
	void SetNoStabilizationTerm();
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
	
	/**
	 * @brief It computes a contribution to the residual vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the residual vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMatPoisson3dReferred::ContributeBC(data,weight,ef,bc);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix &ek,
									 TPZFMatrix &ef);
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix &ef)
	{
		TPZMatPoisson3dReferred::ContributeInterface(data,dataleft,dataright,weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc);
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
	{
		TPZMatPoisson3dReferred::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
	
protected:
    bool fIsReferred;
	
};

#endif
