//
//  TRMPetrophysicsProperties.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMPetrophysicsProperties__
#define __PZ__TRMPetrophysicsProperties__

#include <stdio.h>
#include "pzreal.h"
#include "pzstack.h"
#include "pzfunction.h"

#include "TRMPhaseProperties.h"
#include "TRMWaterPhase.h"
#include "TRMOilPhase.h"
#include "TRMGasPhase.h"

class TRMPetrophysicsProperties{
    
    /**
     * @ingroup Petrophysics Properties
     * @brief Capillary Pressure, relative permeabilities, and one, two and three phase mobility models, Ternary phase diagram respresented by a master triangle
     * @since June 09, 2015
     */
    
private:
    
    /** @brief one-phase flow */
    bool fIsOnePhaseQ;
    
    /** @brief two-phase flow */
    bool fIsTwoPhaseQ;
    
    /** @brief three-phase flow */
    bool fIsThreePhaseQ;
    
    
    /** @brief Definition of the flow system one - two or three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief phase alpha */
    TPZAutoPointer<TRMPhaseProperties> fPhase_alpha;
    
    /** @brief phase beta */
    TPZAutoPointer<TRMPhaseProperties> fPhase_beta;
    
    /** @brief phase gamma */
    TPZAutoPointer<TRMPhaseProperties> fPhase_gamma;
    
public:
    
    /** @brief Default constructor */
    TRMPetrophysicsProperties();

    /** @brief Default constructor */
    ~TRMPetrophysicsProperties();
    
    /** @brief Copy constructor */
    TRMPetrophysicsProperties(const TRMPetrophysicsProperties &copy)
    {
        DebugStop();
    }
    
    /** @brief Copy assignemnt operator $ */
    TRMPetrophysicsProperties &operator=(const TRMPetrophysicsProperties &copy)
    {
        DebugStop();
        return *this;
    }
    
    
    /**
     * @defgroup Access methods
     * @brief    Implements several set/get attributes for the simulation data:
     *
     * @{
     */
    
    /** @brief Set the type of flow system one, two or three phase */
    void SetSystemType(TPZStack<std::string> &SystemType){

#ifdef PZDEBUG
        if(SystemType.size() == 0){
            DebugStop();
        }
#endif
        fSystemType = SystemType;
        switch (fSystemType.size()) {
            case 1:
            {
                fIsOnePhaseQ    = true;
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
    
    // @}
    
    
    /**
     * @defgroup one phase models
     * @brief    Implements Capillary pressure, relative model and mobility for one phase, (Corner of the multiphase diagram of properties)
     *
     * @{
     */
    
    // Capillary Pressure models
    
    /** @brief Alpha-Alpha Capillary Pressure - Pa $P_{caa}$ */
    void Pc(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x);
    
    // Relative permeabilities models
    
    /** @brief Alpha-Alpha phase realtive permeability $k_{ra}$ */
    void Kr(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);
    
    // Mobilities models
    
    /** @brief Alpha-Alpha phase moblity $l_{a}$ */
    void l_single(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    
    // @}
    
    
    /**
     * @defgroup two phase models
     * @brief    Implements Capillary pressure, relative model and mobility for one phase system, (Corner of the multiphase diagram of properties)
     *
     * @{
     */
    
    // Capillary Pressure models
    
    /** @brief Alpha-Alpha Capillary Pressure - Pa $P_{caa}$ */
    void Pcab(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x);
    
    // Relative permeabilities models
    
    /** @brief Alpha-Alpha phase realtive permeability $k_{ra}$ */
    void Kra(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);

    /** @brief Beta-Beta phase realtive permeability $k_{rb}$ */
    void Krb(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);
    
    // Mobilities models
    
    /** @brief Alpha-Alpha phase moblity $l_{a}$ */
    void la(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);

    /** @brief Beta-Beta phase moblity $l_{b}$ */
    void lb(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    /** @brief Alpha-Alpha phase fractional flow $f_{a}$ */
    void fa(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x);

    /** @brief Beta-Beta phase fractional flow $f_{b}$ */
    void fb(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x);
    
    // @}
    
    /**
     * @defgroup three phase models
     * @brief    Implements Capillary pressure, relative model and mobility for two phase system, (Edges of the multiphase diagram of properties)
     *
     * @{
     */
    
    // Capillary Pressure models
    
    /** @brief Alpha-Alpha Capillary Pressure - Pa $P_{caa}$ */
    void Pcab_3p(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x);
    
    void Pcbc_3p(TPZManVector<STATE,10> &pc, TPZManVector<STATE,10> &x);
    
    // Relative permeabilities models
    
    /** @brief Alpha-Alpha phase realtive permeability $k_{ra}$ */
    void Kra_3p(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);
    
    /** @brief Beta-Beta phase realtive permeability $k_{rb}$ */
    void Krb_3p(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);
    
    /** @brief Gamma-Gamma phase realtive permeability $k_{rc}$ */
    void Krc_3p(TPZManVector<STATE,10> &kr, TPZManVector<STATE,10> &x);
    
    // Mobilities models
    
    /** @brief Alpha-Alpha phase moblity $l_{a}$ */
    void la_3p(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    /** @brief Beta-Beta phase moblity $l_{b}$ */
    void lb_3p(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    /** @brief Gamma-Gamma phase moblity $l_{c}$ */
    void lc_3p(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    /** @brief Alpha-Alpha phase fractional flow $f_{a}$ */
    void fa_3p(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x);
    
    /** @brief Beta-Beta phase fractional flow $f_{b}$ */
    void fb_3p(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x);
    
    /** @brief Beta-Beta phase fractional flow $f_{c}$ */
    void fc_3p(TPZManVector<STATE,10> &f, TPZManVector<STATE,10> &x);
    
    
    // @}
    
    
    /**
     * @defgroup total mobiliy
     *
     * @{
     */
    
    /** @brief Total moblity $l_{a}$ */
    void l(TPZManVector<STATE,10> &l, TPZManVector<STATE,10> &x);
    
    // @}
    
    
};

#endif /* defined(__PZ__TRMPetrophysicsProperties__) */
