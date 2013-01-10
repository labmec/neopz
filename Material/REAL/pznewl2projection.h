
/**
 * @file
 * @brief Contains the TPZL2Projection class which implements an L2 projection to use in reconstruction gradient.
 */

//
//  pznewl2projection.h
//  PZ
//
//  Created by Agnaldo Farias on 1/9/13.
//
//



#ifndef PZ_pznewl2projection_h
#define PZ_pznewl2projection_h

#include "pzmaterial.h"
#include "pzdiscgal.h"

class TPZL2ProjectionForGradient : public TPZDiscontinuousGalerkin{

protected:
    
    /** @brief Problem dimension */
	int fDim;
	
	/** @brief Number of state variables */
	int fNStateVars;
	
	/**
     *@brief data of gradient vector:
     * first column: index of the element
     * second column: midpoint of the element (x0, y0, z0)
     * other columns: components of the gradient reconstructed (dudx, dudy, dudz)
     */
	TPZFMatrix<REAL> fgradients;
    
    public:
    
    /**
	 * @brief Class constructor
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param gradients data of gradient vector
	 */
	TPZL2ProjectionForGradient(int id, int dim,int nvar);
	
	/** @brief Default destructor */
	~TPZL2ProjectionForGradient();
	    
	
	void SetGradients(TPZFMatrix<REAL> &grad) {
		fgradients = grad;
	}
	
    /** @brief Returns problem dimension */
	virtual int Dimension(){ return this->fDim; }
	
	/** @brief Returns number of state variables */
	virtual int NStateVariables(){ return this->fNStateVars; }
    
};
#endif
