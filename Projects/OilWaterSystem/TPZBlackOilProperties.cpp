/**
 * @file   TPZBlackOilProperties.h
 * @Author Omar Duran Triana. omaryesiduran@gmail.com
 * @date   March, 2015
 * @brief  Brief this file implements a blackoil correlations in Field units.
 *
 * Define correlations and procedures to compute BlackOil model properties.
 * Fanchi, J. R., 2006. Petroleum Engineering Handbook, Volume I: General Engineering. s.l.:SPE books.
 * McCain .W. 1990. The properties of petroleum fluids. Second edition. PennWell Books.
 */

#include "TPZBlackOilProperties.h"


TPZBlackOilProperties::TPZBlackOilProperties()
{
    /** @brief offset value for central finite difference derivative approximation */
    fepsilon = 1.0e-8;
    
    /** @brief Field surface pressure condition 14.7 Psi*/
    fPressure_sc = 14.7;
    
    /** @brief Field surface temperature condition 60 F*/
    fTemperature_sc = 60.0;
    
    /** @brief Range of operational Pressure Psi [1.0 , 55000.0] (Russia -> Kola Superdeep Borehole 12.345 km P approximates to 55000 Psi )*/
    fOperationalPressure.Resize(2);
    fOperationalPressure[0] = 1.0;
    fOperationalPressure[1] = 55000.0;    
    
    /** @brief Range of operational Temperature F [1.0 , 350.0] (Russia -> Kola Superdeep Borehole 12.345 km T approximates to 350 F )*/
    fOperationalTemperature.Resize(2);
    fOperationalTemperature[0] = 1.0;
    fOperationalTemperature[1] = 350.0;    
    
    /** @brief Brine density at standard conditions Defined form McCain as 62.368 lbm/feet3 */
    fDensity_w_sc = 62.368;   
    
    /** @brief Set of correlations models for water properties */
    fCorrelationSet_w.Resize(8);
    fCorrelationSet_w.Fill(0);
    
    /** @brief Set of correlations models for oil properties */
    fCorrelationSet_o.Resize(8);
    fCorrelationSet_o.Fill(0);    
    
    /** @brief Set of correlations models for gas properties */
    fCorrelationSet_g.Resize(8);
    fCorrelationSet_g.Fill(0);    
}

TPZBlackOilProperties::~TPZBlackOilProperties()
{
	
}

    /**
    * Print a error message if some condition is break but acceptable
    * @return Warning message.
    */    
    void TPZBlackOilProperties::TPZBlackOilPropertiesError()
    {
        std::cout << "Error in blackoil properties." << std::endl; 
    }
    
    /**
    * Print a error message if some condition is break but acceptable
    * @return Warning message.
    */     
    void TPZBlackOilProperties::TPZBlackOilPropertiesWarning()
    {
        std::cout << "Warning in blackoil properties." << std::endl; 
    }
    
    
    
    
    
    
    
    
    
    /**
    * @defgroup Water Properties of water from few field parameters
    * @brief    Implements several correlations to water properties from the following parameters:
    *           Water Salinity
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
    REAL TPZBlackOilProperties::B_w(REAL &P, REAL &T)
    {
        return -1;
    }
    
    /**
    * Computes the specific gravity of water @ sc
    * @param watersalinity total dissolved solids (percent).
    * @see gamma_w(REAL &watersalinity)
    * @return specific gravity of water (-)
    */
    REAL TPZBlackOilProperties::gamma_w(REAL &watersalinity)
    {
        if((watersalinity < 0.0 ) || (watersalinity > 1.0))
        {
            this->TPZBlackOilPropertiesError();
            std::cout << "Total dissolved solid is out the range [0,1]. watersalinity = " << watersalinity << std::endl; 
            DebugStop();
        }
        
        REAL a0, a1, a2, grammabrine;
        a0 = fDensity_w_sc;
        a1 = 0.438603;
        a2 = 0.00160074;  
        grammabrine = 1 + (a1 * watersalinity + a2 * watersalinity * watersalinity) / a0;
        return grammabrine;
    }      
    
    /**
    * Computes the density of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see Density_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Density of water (lbm/feet3 or kg/m3)
    */
    REAL TPZBlackOilProperties::Density_w(REAL &P, REAL &T, REAL &watersalinity)
    {
        return -1;
    }
    
    /**
    * Computes the compressibility of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see c_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Compressibility of water (1/psi or 1/Pa)
    */
    REAL TPZBlackOilProperties::c_w(REAL &P, REAL &T, REAL &watersalinity)
    {
        return -1;
    }
    
    /**
    * Computes the viscosity of water @ P standard and T conditions
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see muPsc_w(REAL &T, REAL &watersalinity)
    * @return Viscosity of water at P standard (cp or Pa plus second)
    */
    REAL TPZBlackOilProperties::muPsc_w(REAL &T, REAL &watersalinity)
    {
        if((watersalinity < 0.0 ) || (watersalinity > 1.0))
        {
            this->TPZBlackOilPropertiesError();
            std::cout << "Total dissolved solid is out the range [0,1]. watersalinity = " << watersalinity << std::endl; 
            DebugStop();
        }
        
        if((T < fOperationalTemperature[0] ) || (T > fOperationalTemperature[1]))
        {
            this->TPZBlackOilPropertiesError();
            std::cout << "Temperature is out the range operational range. " << fOperationalTemperature << std::endl; 
            DebugStop();
        }
              
        
        REAL a0, a1, a2, a3, b0, b1, b2, b3, b4, fA, fB, fmuPsc_w;
        REAL S2,S3,S4;
        S2 = watersalinity * watersalinity;
        S3 = watersalinity * S2;
        S4 = watersalinity * S3;
        a0 = 109.527;
        a1 = -8.40564;
        a2 = 0.313314;
        a3 = 0.00872213;
        b0 = -1.12166;
        b1 = 0.0263951;
        b2 = -0.0000679461;
        b3 = -0.0000547119;
        b4 = 0.00000155586;        
        fA = a0 + a1 * watersalinity + a2 * S2 + a3 * S3;
        fB = b0 + b1 * watersalinity + b2 * S2 + b3 * S3 + b4 * S4;
        fmuPsc_w = fA * pow(T,fB);
        
        return fmuPsc_w;
    }
    
    /**
    * Computes the viscosity of water @ P and T conditions
    * @param P Water Pressure (Psi or Pa).
    * @param T Water Temperature (F or C).
    * @param watersalinity total dissolved solids (percent).
    * @see mu_w(REAL &P, REAL &T, REAL &watersalinity)
    * @return Viscosity of water (cp or Pa plus second)
    */
    REAL TPZBlackOilProperties::mu_w(REAL &P, REAL &T, REAL &watersalinity)
    {
        return -1;
    }
    
    
    
    
    
    
    
    
    
    
    
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
