//
//  mixedpoisson.h
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_mixedpoisson_h
#define PZ_mixedpoisson_h

#include "TPZMaterial.h"

#include "pzpoisson3d.h"
#include "TPZMaterial.h"

/**
 * @ingroup material
 * @author Agnaldo Farias
 * @since 5/28/2012
 * @brief Material to solve a mixed poisson problem 2d by multiphysics simulation
 * @brief Pressure(p): uses H1 space.  Velocity (Q): uses Hdiv space
 */


class TPZHybridPoisson : public TPZMatPoisson3d{
    
protected:
    
public:
    TPZHybridPoisson();
    
    TPZHybridPoisson(int matid, int dim);
    
    virtual ~TPZHybridPoisson();
    
    virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZHybridPoisson"; }
    
	/** 
	 * @brief Fill material data parameter with necessary requirements for the ContributeInterface method.
     * @since April 10, 2007
	 */
	/** 
	 * Here, in base class, all requirements are considered as necessary. \n
	 * Each derived class may optimize performance by selecting only the necessary data.
	 */
	virtual void FillDataRequirementsInterface(TPZMaterialData &data)
    {
        data.fNeedsNormal = true;
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
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
    {
        
    }
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ef)
    {
        
    }

    	
};

#endif
