/**
 * @file
 * @brief Contains the TPZMatPoisson3dReferred class which implements a version of TPZMatPoisson3d \n
 * (convection term is given at each integration point)
 */

#ifndef MATPOISSON3DREFERREDH
#define MATPOISSON3DREFERREDH

#include <iostream>
#include "pzpoisson3d.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief This class implements a version of TPZMatPoisson3d where the convection term is given \n
 * at each integration point from a previous calculation.
 */
/**
 * The convection term \f$ fC * fConvDir = - fAlpha * grad(sol) \f$ where \f$ grad(sol) \f$ is the gradient of the previous solution.
 */
class TPZMatPoisson3dReferred : public TPZMatPoisson3d {
	
protected:
	
	REAL falpha;
	
	/** @brief Sets convection term */
	void SetConvectionTerm(TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &axes);
	
	/** @brief Sets convection term for ContributeInterface methods */
	/** It expect dsolL and dsolR to be dSol/dX, i.e. the derivatives with respect to the global coordinates. */
	void SetConvectionTermInterface(TPZFMatrix<STATE> &dsolL, TPZFMatrix<STATE> &dsolR);
	
public:
	
	TPZMatPoisson3dReferred(int nummat, int dim);
	
	virtual ~TPZMatPoisson3dReferred();
	
	TPZMatPoisson3dReferred(const TPZMatPoisson3dReferred &copy) : 
    TPZRegisterClassId(&TPZMatPoisson3dReferred::ClassId), TPZMatPoisson3d(copy){
		this->falpha = copy.falpha;
	}
	
	virtual TPZMaterial * NewMaterial() override{
		return new TPZMatPoisson3dReferred(*this);
	}
    
    virtual int NStateVariables() const override
    {
        return 1;
    }
	
	void SetAlpha(REAL alpha){
		this->falpha = alpha;
	}
	
	REAL GetAlpha(){
		return this->falpha;
	}
    
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<STATE> &ek, 
                            TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
							  TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCond &bc) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ek,
									 TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override;
    
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ef) override
	{
		TPZMatPoisson3d::Contribute(data,weight,ef);
	}
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZMatPoisson3d::ContributeBC(data,weight,ef,bc);
	}
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ef) override
	{
		TPZMatPoisson3d::ContributeInterface(data,dataleft, dataright, weight,ef);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override
	{
		TPZMatPoisson3d::ContributeBCInterface(data,dataleft,weight,ef,bc);
	}
    virtual int ClassId() const override;
};

#endif
