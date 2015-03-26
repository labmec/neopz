#ifndef TPZBlackOilPropertiesH
#define TPZBlackOilPropertiesH
/**
 * @file   TPZBlackOilProperties.h
 * @Author Omar Duran Triana. omaryesiduran@gmail.com
 * @date   March, 2015
 * @brief  Brief this file implements a blackoil correlations in Field units.
 *
 * Define correlations and procedures to compute BlackOil model properties.
 * Fanchi, J. R., 2006. Petroleum Engineering Handbook, Volume I: General Engineering. s.l.:SPE books.
 * McCain .W. 1990. The properties of petroleum fluids. Second edition. PennWell Books.
 * 
 */

#include "tpzautopointer.h"
#include "pzmanvector.h"
#include <math.h>

class TPZBlackOilProperties {
	
protected:

	/** @brief offset value for central finite difference derivative approximation */
 	REAL fepsilon;	
    
    /** @brief Field surface pressure condition */
    REAL fPressure_sc;
    
    /** @brief Field surface temperature condition */
    REAL fTemperature_sc;    
    
    /** @brief Range of operational Pressure */
    TPZManVector<REAL> fOperationalPressure;
    
    /** @brief Range of operational Pressure */
    TPZManVector<REAL> fOperationalTemperature;   
    
    /** @brief Brine density at standard conditions */
    REAL fDensity_w_sc;
    
    /** @brief Set of correlations models for water properties */
    TPZManVector<REAL> fCorrelationSet_w;
    
    /** @brief Set of correlations models for oil properties */
    TPZManVector<REAL> fCorrelationSet_o;
    
    /** @brief Set of correlations models for gas properties */
    TPZManVector<REAL> fCorrelationSet_g;    
       
	
private:
    
    /**
    * Print a error message if some condition is break
    * @return Error message.
    */    
    void TPZBlackOilPropertiesError();
    
    /**
    * Print a error message if some condition is break but acceptable
    * @return Warning message.
    */    
    void TPZBlackOilPropertiesWarning();
    
    
public:
	
	TPZBlackOilProperties();
    
	~TPZBlackOilProperties();
    
    
    /**
    * @defgroup Set_Get Set and get Methods
    * @{
    */
    
    void Setepsilon(REAL & epsilon){ fepsilon = epsilon; }  
    REAL Getepsilon(){ return fepsilon; }
    
    void SetPressure_sc(REAL & Pressure_sc){ fPressure_sc = Pressure_sc; }  
    REAL GetPressure_sc(){ return fPressure_sc; }
    
    void SetTemperature_sc(REAL & Temperature_sc){ fTemperature_sc = Temperature_sc; }  
    REAL GetTemperature_sc(){ return fTemperature_sc; }    
    
    void SetDensity_w_sc(REAL & Density_w_sc){ fDensity_w_sc = Density_w_sc; }  
    REAL GetDensity_w_sc(){ return fDensity_w_sc; }       
    
    // @}
    
    
    

    /**
    * @defgroup Water Properties of water from few field parameters
    * @brief    Implements several correlations to water properties from the following parameters:
    *           Water Pressure (Psi or Pa)
    *           Water Temperature (F  or C)
    *           Water Salinity (percent)
    *
    * @{
    */

    /**
    * Computes the formation volume factor of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @see B_w(REAL &P, REAL &T)
    * @return Water formation volume factor (feet3@PT/feet3@sc or m3@PT/m3@sc)
    */
    REAL B_w(REAL &P, REAL &T);
    
    /**
    * Computes the specific gravity of water @ sc
    * @param watersalinity total dissolved solids (percent).
    * @see gamma_w(REAL &watersalinity)
    * @return specific gravity of water (-)
    */
    REAL gamma_w(REAL &watersalinity);      
    
    /**
    * Computes the density of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see Density_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Density of water (lbm/feet3 or kg/m3)
    */
    REAL Density_w(REAL &P, REAL &T, REAL &watersalinity);     
    
    /**
    * Computes the compressibility of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see c_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Compressibility of water (1/psi or 1/Pa)
    */
    REAL c_w(REAL &P, REAL &T, REAL &watersalinity);
    
    /**
    * Computes the viscosity of water @ P standard and T conditions
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see muPsc_w(REAL &T, REAL &watersalinity)
    * @return Viscosity of water at P standard (cp or Pa plus second)
    */
    REAL muPsc_w(REAL &T, REAL &watersalinity);  
    
    /**
    * Computes the viscosity of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see mu_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Viscosity of water (cp or Pa plus second)
    */
    REAL mu_w(REAL &P, REAL &T, REAL &watersalinity);     
    
    
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