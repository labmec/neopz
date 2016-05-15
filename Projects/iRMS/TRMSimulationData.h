//
//  TRMSimulationData.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMSimulationData__
#define __PZ__TRMSimulationData__

#include <stdio.h>
#include "pzreal.h"
#include "TRMBlackOilProperties.h"
#include "TRMRockProperties.h"
#include "TRMPetrophysicsProperties.h"
#include "TRMSpatialPropertiesMap.h"
#include "TRMRawData.h"


class TRMSimulationData {
    
protected:
    
    
    /** @brief current time state */
    bool fIsCurrentStateQ;
    
    /** @brief one-phase flow */
    bool fIsOnePhaseQ;
    
    /** @brief two-phase flow */
    bool fIsTwoPhaseQ;
    
    /** @brief three-phase flow */
    bool fIsThreePhaseQ;
    
    /** @brief Definition of the flow system one - two or three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief Store time values to be reported */
    TPZStack< STATE , 500 > fReportingTimes;
    
    /** @brief ntime steps */
    int fn_steps;
    
    /** @brief Initial time */
    STATE ftime_0;
    
    /** @brief Final time */
    STATE ftime_n;
    
    /** @brief Time step */
    STATE fdt;
    
    /** @brief Min time step */
    STATE fdt_min;
    
    /** @brief Max time step */
    STATE fdt_max;
    
    /** @brief Increment dt factor */
    STATE fdt_up;
    
    /** @brief Decrement dt factor */
    STATE fdt_down;
    
    /** @brief number of corrections steps */
    int fn_corrections;
    
    /** @brief residue overal tolerance */
    STATE fepsilon_res;
    
    /** @brief correction overal tolerance */
    STATE fepsilon_cor;
    
    /** @brief Autopointer of the RawData */
    TPZAutoPointer<TRMRawData> fRawData;
    
    /** @brief Stores all the petrophysics data */
    TRMPetrophysicsProperties fPetroProp;
    
    /** @brief Stores all the blackoil fluid properties */
    TRMBlackOilProperties fBlackOilProp;
    
    /** @brief Stores all the rock and geomechanic properties */
    TRMRockProperties fRockProp;
    
    /** @brief Stores the spatial information given in maps */
    TRMSpatialPropertiesMap fSpatialProp;
    
    
public:
    
    
    /** @brief default constructor */
    TRMSimulationData();
    
    /** @brief destructor */
    ~TRMSimulationData();
    
    /**
     * @defgroup Access methods
     * @brief    Implements several set/get attributes for the simulation data:
     *
     * @{
     */
    
    /** @brief current time state */
    void SetCurrentStateQ(bool state) { fIsCurrentStateQ = state; }
    
    /** @brief current time state */
    bool IsCurrentStateQ() {return fIsCurrentStateQ;}
    
    /** @brief Mono-phasic system */
    bool IsOnePhaseQ() {return fIsOnePhaseQ;}
    
    /** @brief Two-phasic system */
    bool IsTwoPhaseQ() {return fIsTwoPhaseQ;}
    
    /** @brief Three-phasic system */
    bool IsThreePhaseQ() {return fIsThreePhaseQ;}
    
    /** @brief Definition of the flow system one - two and  ... three phase */
    void SetSystemType(){
        
        switch (fSystemType.size()) {
            case 1:
            {
                fIsOnePhaseQ = true;
            }
                break;
            case 2:
            {
                fIsTwoPhaseQ = true;
                
            }
                break;
            case 3:
            {
                fIsThreePhaseQ = true;
            }
                break;
                
            default:
            {
                std::cout << "This code run just three-phasic systems" << std::endl;
                DebugStop();
            }
                break;
        }
    }
    
    /** @brief Setup reporting times and time step size */
    void SetTimeControls(int n_times, STATE dt, STATE dt_in, STATE dt_de);
    
    /** @brief Setup reporting times and time step size */
    void SetNumericControls(int n_corrections, STATE epsilon_res, STATE epsilon_cor);
    
    /** @brief Store time values to be reported */
    TPZStack< STATE , 500 > ReportingTimes(){
        return fReportingTimes;
    }
    
    /** @brief Initial time */
    STATE time_0() { return ftime_0; }
    
    /** @brief Final time */
    STATE time_n() { return ftime_n; }
    
    /** @brief Time step */
    STATE dt() { return fdt; }
    
    /** @brief Min time step */
    STATE dt_min() { return fdt_min; }
    
    /** @brief Max time step */
    STATE dt_max() { return fdt_max; }
    
    /** @brief Increment dt factor */
    STATE dt_up() { return fdt_up; }
    
    /** @brief Decrement dt factor */
    STATE dt_down() { return fdt_down; }
    
    /** @brief number of corrections steps */
    int n_corrections() { return fn_corrections; }
    
    /** @brief residue overal tolerance */
    STATE epsilon_res() { return fepsilon_res; }
    
    /** @brief correction overal tolerance */
    STATE epsilon_cor() { return fepsilon_cor; }
    
    /** @brief Set autopointer of the RawData */
    void SetRawData(TPZAutoPointer<TRMRawData> &RawData);
    
    /** @brief Get autopointer of the RawData */
    TPZAutoPointer<TRMRawData>  RawData(){
        return fRawData;
    }
    
    /** @brief Stores all the petrophysics data */
    void SetPetroProp(TRMPetrophysicsProperties PetroProp)
    {
        fPetroProp = PetroProp;
    }
    
    TRMPetrophysicsProperties & PetroProp()
    {
        return fPetroProp;
    }

    /** @brief Stores all the blackoil fluid properties */
    void SetBlackOilProp(TRMBlackOilProperties BlackOilProp)
    {
        fBlackOilProp = BlackOilProp;
    }
    
    TRMBlackOilProperties & BlackOilProp()
    {
        return fBlackOilProp;
    }
  
    /** @brief Stores all the rock and geomechanic properties */
    void SetRockProp(TRMRockProperties RockProp)
    {
        fRockProp = RockProp;
    }
    
    TRMRockProperties & RockProp()
    {
        return fRockProp;
    }
    
    /** @brief Stores the spatial information given in maps */
    void SetSpatialProp(TRMSpatialPropertiesMap SpatialProp)
    {
        fSpatialProp = SpatialProp;
    }
    
    TRMSpatialPropertiesMap & GetSpatialProp()
    {
        return fSpatialProp;
    }
    
    static STATE Gravity();
    
    // @}
    

};

#endif /* defined(__PZ__TRMSimulationData__) */
