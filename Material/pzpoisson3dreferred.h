/**
 * \file
 * @brief Contains the TPZMatPoisson3dReferred class which implements a version of TPZMatPoisson3d \n
 * (convection term is given at each integration point)
 */

//$Id: pzpoisson3dreferred.h,v 1.6 2009-09-01 19:44:48 phil Exp $

#ifndef MATPOISSON3DREFERREDH
#define MATPOISSON3DREFERREDH

#include <iostream>
#include "pzpoisson3d.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief This class implements a version of TPZMatPoisson3d where the convection term is given at each integration point
 * from a previous calculation.
 */
/**
 * The convection term fC * fConvDir = - fAlpha * grad(sol) where grad(sol) is the gradient of the previous solution.
 */
class TPZMatPoisson3dReferred : public TPZMatPoisson3d {
	
protected:
	
	REAL falpha;
	
	/** SetConvectionTerm
	 */
	void SetConvectionTerm(TPZFMatrix &dsol, TPZFMatrix &axes);
	
	/** SeConvectionTerm for ContributeInterface methods
	 * It expect dsolL and dsolR to be dSol/dX, i.e. the derivatives 
	 * with respect to the global coordinates.
	 */
	void SetConvectionTermInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR);
	
public:
	
	TPZMatPoisson3dReferred(int nummat, int dim);
	
	virtual ~TPZMatPoisson3dReferred();
	
	//  virtual int HasForcingFunction() {return true;}
	
	TPZMatPoisson3dReferred(const TPZMatPoisson3dReferred &copy) : TPZMatPoisson3d(copy){
		this->falpha = copy.falpha;
	}
	
	virtual TPZAutoPointer<TPZMaterial> NewMaterial(){
		return new TPZMatPoisson3dReferred(*this);
	}
	
	void SetAlpha(REAL alpha){
		this->falpha = alpha;
	}
	
	REAL GetAlpha(){
		return this->falpha;
	}
    
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek, 
                            TPZFMatrix &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
							  TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ek,
									 TPZFMatrix &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc);
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef)
	{
		TPZMatPoisson3d::Contribute(data,weight,ef);
	}
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMatPoisson3d::ContributeBC(data,weight,ef,bc);
	}
	virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ef)
	{
		TPZMatPoisson3d::ContributeInterface(data,weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
	{
		TPZMatPoisson3d::ContributeBCInterface(data,weight,ef,bc);
	}
	
};

#endif
