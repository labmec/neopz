//$Id: pzporoelastoplasticmem.h,v 1.3 2009-10-06 01:00:06 erick Exp $

#ifndef PZPOROELASTOPLASTICMEM_H
#define PZPOROELASTOPLASTICMEM_H

#include "TPZMaterial.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"
#include "TPZElastoPlasticMem.h"

  /**
   * This class defines the material memory that should be stored at each integration point
   * for the purposes of an elastoplastic material.
   */

class TPZPoroElastoPlasticMem : public TPZElastoPlasticMem
{
public:
	TPZPoroElastoPlasticMem();
	
	TPZPoroElastoPlasticMem(const TPZPoroElastoPlasticMem & source);
	
	const TPZPoroElastoPlasticMem & operator=(const TPZPoroElastoPlasticMem & source);
	
	virtual ~TPZPoroElastoPlasticMem();
	
	const std::string Name()const;
	
	public:
int ClassId() const override;

	
    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;

	virtual void Print(std::ostream &out = std::cout)const override;
	
	/**
	 * Operator<<
	 */
    friend std::ostream& operator<<( std::ostream& Out, const TPZPoroElastoPlasticMem & s )
	{
		s.Print(Out);
		return Out;
	}
	
	/**
	 * Total Pore Pressure
	 */
	REAL fPorePressure;
    
	/** 
	 * Spatial Gradient of Pore Pressure
	 */
	TPZVec<REAL> fdPorePressure;
    
    /**
     * Darcy's velocity
     */
    TPZVec<REAL> fv;
	
};



#endif
