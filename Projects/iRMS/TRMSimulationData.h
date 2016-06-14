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
#include "TRMPhaseProperties.h"
#include "TRMRockProperties.h"
#include "TRMPetrophysicsProperties.h"
#include "TRMSpatialPropertiesMap.h"
#include "TRMRawData.h"


class TRMSimulationData {
    
protected:

    /** @brief initial state */
    bool fIsInitialStateQ;
    
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
    
    /** @brief Definition gravity field */
    TPZVec<STATE> fg;
    
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
    
    /** @brief Time step */
    STATE ftime;
    
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
    
    /** @brief Autopointer of all the petrophysics data */
    TPZAutoPointer<TRMPetrophysicsProperties> fPetroPhysics;
    
    /** @brief phase alpha */
    TPZAutoPointer<TRMPhaseProperties> fPhase_alpha;
    
    /** @brief phase beta */
    TPZAutoPointer<TRMPhaseProperties> fPhase_beta;
    
    /** @brief phase gamma */
    TPZAutoPointer<TRMPhaseProperties> fPhase_gamma;
    
    /** @brief Stores all the rock and geomechanic properties */
    TRMRockProperties fRockProp;
    
    /** @brief Stores the spatial information given in maps */
    TPZAutoPointer<TRMSpatialPropertiesMap> fMap;
    
    
public:
    
    
    /** @brief default constructor */
    TRMSimulationData();
    
    /** @brief default constructor */
    TRMSimulationData(const TRMSimulationData &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMSimulationData &operator=(const TRMSimulationData &copy)
    {
        DebugStop();
        return *this;
    }
    
    /** @brief destructor */
    ~TRMSimulationData();
    
    /**
     * @defgroup Access methods
     * @brief    Implements several set/get attributes for the simulation data:
     *
     * @{
     */

    /** @brief Set initial state */
    void SetInitialStateQ(bool state) { fIsInitialStateQ = state; }
    
    /** @brief Get initial state */
    bool IsInitialStateQ() {return fIsInitialStateQ;}
    
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
    
    /** @brief Set the type of flow system one, two or three phase */
    void SetSystemType(TPZStack<std::string> &SystemType, TPZStack< TPZAutoPointer<TRMPhaseProperties> > &Phases){

        if(Phases.size() != SystemType.size()){
            DebugStop();
        }
        
        TPZAutoPointer<TRMPetrophysicsProperties> petro_physics = new TRMPetrophysicsProperties;
        petro_physics->SetSystemType(SystemType);
        switch (SystemType.size()) {
            case 1:
            {
                fIsOnePhaseQ    = true;
                SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseAlpha(Phases[0]);
            }
                break;
            case 2:
            {
                fIsTwoPhaseQ = true;
                SetPhaseAlpha(Phases[0]);
                SetPhaseBeta(Phases[1]);
        
                petro_physics->SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseBeta(Phases[1]);
            }
                break;
            case 3:
            {
                fIsThreePhaseQ = true;
                SetPhaseAlpha(Phases[0]);
                SetPhaseBeta(Phases[1]);
                SetPhaseGamma(Phases[2]);
                
                petro_physics->SetPhaseAlpha(Phases[0]);
                petro_physics->SetPhaseBeta(Phases[1]);
                petro_physics->SetPhaseGamma(Phases[2]);
            }
                break;
                
            default:
            {
                std::cout << "This code run just three-phasic systems" << std::endl;
                DebugStop();
            }
                break;
        }
        fSystemType = SystemType;
        fPetroPhysics  = petro_physics;
        
    }
    
    /** @brief Get the type of flow system one, two or three phase */
    TPZStack<std::string>  SystemType(){
        return fSystemType;
    }
    
    /** @brief Set phase alpha */
    void SetPhaseAlpha(TPZAutoPointer<TRMPhaseProperties> &alpha);
    
    /** @brief Get phase alpha */
    TPZAutoPointer<TRMPhaseProperties> & AlphaProp();
    
    /** @brief Set phase beta */
    void SetPhaseBeta(TPZAutoPointer<TRMPhaseProperties> &beta);
    
    /** @brief Get phase beta */
    TPZAutoPointer<TRMPhaseProperties> & BetaProp();
    
    /** @brief Set phase gamma */
    void SetPhaseGamma(TPZAutoPointer<TRMPhaseProperties> &gamma);
    
    /** @brief Get phase gamma */
    TPZAutoPointer<TRMPhaseProperties> & GammaProp();
    
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
    
    /** @brief Set Time step */
    void Setdt(STATE dt) { fdt = dt; }
    
    /** @brief Time step */
    STATE dt() { return fdt; }

    /** @brief Time */
    void SetTime(STATE time) { ftime = time; }
    
    /** @brief Time */
    STATE t() { return ftime; }
    
    /** @brief Min time step */
    STATE dt_min() { return fdt_min; }
    
    /** @brief Max time step */
    STATE dt_max() { return fdt_max; }
    
    /** @brief Increment dt factor */
    STATE dt_up() { return fdt_up; }
    
    /** @brief Decrement dt factor */
    STATE dt_down() { return fdt_down; }
    
    /** @brief number of corrections steps */
    int n_steps() { return fn_steps; }
    
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
    void SetPetroPhysics(TPZAutoPointer<TRMPetrophysicsProperties> &PetroPhysics)
    {
        fPetroPhysics = PetroPhysics;
    }
    
    TPZAutoPointer<TRMPetrophysicsProperties> & PetroPhysics()
    {
        return fPetroPhysics;
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
    void SetMap(TPZAutoPointer<TRMSpatialPropertiesMap> &SpatialProp)
    {
        fMap = SpatialProp;
    }
    
    TPZAutoPointer<TRMSpatialPropertiesMap> & Map()
    {
        return fMap;
    }
    
    void SetGravity(TPZVec<STATE> &g)
    {
        fg = g;
    }
    
    TPZVec<STATE> & Gravity()
    {
        return fg;
    }
    
    // @}
    

};

#endif /* defined(__PZ__TRMSimulationData__) */
