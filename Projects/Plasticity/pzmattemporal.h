//$Id: pzmattemporal.h,v 1.4 2009-10-05 03:49:58 erick Exp $

#ifndef PZMATTEMPORAL_H
#define PZMATTEMPORAL_H


#include "pzmaterial.h"
#include "pzconslaw.h" 

/**
 * Interface setup for all classes involving temporal behaviour
 */

class TPZMatTemporal : public TPZMaterialData
{
	public:
	TPZMatTemporal():fDeltaT(1.), fTime(Advanced_CT)
	 {; };
	
	virtual ~TPZMatTemporal() 
	 {};
	
	virtual void SetDeltaT(const REAL deltaT)
	 { fDeltaT = deltaT; };
	
	virtual void SetContributionTime(TPZContributeTime time)
	 { fTime = time; };
		
	protected:
	/**
	 * Time lapse in the temporal integrator
	 */
	REAL fDeltaT;	
	
	/**
	 * Some materials require double temporal contributions in order
	 * to contribute converged (Last_CT) and current (Advanced_CT)
	 * in the temporal scheme.
	 */
	TPZContributeTime fTime;
	/**
	 * Whether to compute the initial solution vector based on a constant initial value
	 */
};

#endif
