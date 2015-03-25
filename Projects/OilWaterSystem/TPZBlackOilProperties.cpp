/**
 * @file   TPZBlackOilProperties.h
 * @Author Omar Duran Triana. omaryesiduran@gmail.com
 * @date   March, 2015
 * @brief  Brief this file implements a blackoil correlations in Field units.
 *
 * Define correlations and procedures to compute BlackOil model properties.
 * Fanchi, J. R., 2006. Petroleum Engineering Handbook, Volume I: General Engineering. s.l.:SPE books. 
 */

#include "TPZBlackOilProperties.h"


TPZBlackOilProperties::TPZBlackOilProperties()
{
    /** @brief offset value for central finite difference derivative approximation */
    fepsilon = 1.0e-8;  
}

TPZBlackOilProperties::~TPZBlackOilProperties()
{
	
}