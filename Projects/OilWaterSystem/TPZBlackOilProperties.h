#ifndef TPZProblemDATAH
#define TPZProblemDATAH
/**
 * @file   TPZBlackOilProperties.h
 * @Author Omar Duran Triana. omaryesiduran@gmail.com
 * @date   March, 2015
 * @brief  Brief this file implements a blackoil correlations in Field units.
 *
 * Define correlations and procedures to compute BlackOil model properties.
 * Fanchi, J. R., 2006. Petroleum Engineering Handbook, Volume I: General Engineering. s.l.:SPE books. 
 */

#include "tpzautopointer.h"
#include <math.h>

class TPZBlackOilProperties {
	
protected:

	/** @brief offset value for central finite difference derivative approximation */
	REAL fepsilon;	
	
private:

    
public:
	
	TPZBlackOilProperties();
	
	~TPZBlackOilProperties();

    /**
    * @defgroup Water Properties of water from few field parameters
    * @brief    Implements several correlations to water properties from the following parameters:
    *           Water Salinity
    *
    * @{
    */

    /**
    * Computes the water density @ sc
    * @param watersalinity total dissolved solids (percent).
    * @see Densitysc_w(watersalinity)
    * @return Water density (lbm/feet3)
    */
    REAL Density_w(REAL &watersalinity);    
    
    // @}
    
    /**
    * @defgroup Oil Properties of oil from few field parameters
    * @brief    Implements several correlations to oil properties from the following parameters:
    * 
    * 
    * @{
    */

    // @}
    
    /**
    * @defgroup Gas Properties of gas from few field parameters
    * @brief    Implements several correlations to gas properties from the following parameters:
    * 
    * 
    * @{
    */

    // @}
	
	
};


#endif