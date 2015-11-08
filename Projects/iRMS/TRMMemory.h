//
//  TRMMemory.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
//

#ifndef __PZ__TRMMemory__
#define __PZ__TRMMemory__

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"
#include "pzfmatrix.h"


/*! @brief Store the information required on a integration point.
 *         Brief description continued.
 *
 *  Store the saturation at n step.
 *  Also it can store the nonlinear part of the flux at n step.
 *  Store the xyz of the spatial properties.
 *  Gravitational Segregation terms
 */

class TRMMemory {

    /**
     * @defgroup Compute and get matrices
     * @{
     */
    
    /** @brief Total flux */
    TPZManVector<STATE> fu;

    /** @brief Total flux at the previous timestep */
    TPZManVector<STATE> fu_n;
    
    /** @brief Weighted Pressure */
    STATE fPressure;
    
    /** @brief Weighted Pressure at the previous timestep */
    STATE fPressure_n;
    
    /** @brief Water Saturation */
    STATE fSw;
    
    /** @brief Water saturation at the previous timestep */
    STATE fSw_n;
    
    /** @brief Rock Porosity */
    STATE fporosity;
    
    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> K;
    
    //@}
    
public:
    
    /** @brief Default constructor */
    TRMMemory();
    
    /** @brief Default destructor */
    ~TRMMemory();
    

    void UpdateSolutionMemory()
    {
        //update saturation and pressure and total flux (un = unp1)
        fPressure_n = fPressure;
    }
    
    /**
     * @defgroup Set Get methods
     * @{
     */
    
    /** @brief Set the total flux */
    void SetTotal_Flux(TPZManVector<STATE> u){
        fu = u;
    }
    
    /** @brief Get the total flux */
    TPZManVector<STATE> GetTotal_Flux(){
        return fu;
    }
    
    /** @brief Set the total flux at the previous timestep */
    void SetTotal_Flux_n(TPZManVector<STATE> u_n){
        fu_n = u_n;
    }
    
    /** @brief Get the total flux at the previous timestep */
    TPZManVector<STATE> GetTotal_Flux_n(){
        return fu_n;
    }
    
    /** @brief Set the weighted pressure */
    void SetPressure(STATE p){
        fPressure = p;
    }
    
    /** @brief Get the weighted pressure */
    STATE GetPressure(){
        return fPressure;
    }
    
    /** @brief Set the weighted pressure at the previous timestep */
    void SetPressure_n(STATE p_n){
        fPressure_n = p_n;
    }
    
    /** @brief Get the weighted pressure at the previous timestep */
    STATE GetPressure_n(){
        return fPressure_n;
    }
    
    /** @brief Set the water saturation */
    void SetSaturation(STATE Sw){
        fSw = Sw;
    }
    
    /** @brief Get the water saturation */
    STATE GetSaturation(){
        return fSw;
    }
    
    /** @brief Set the water saturation at the previous timestep */
    void SetSw_n(STATE Sw_n){
        fSw_n = Sw_n;
    }
    
    /** @brief Get the water saturation at the previous timestep */
    STATE GetSw_n(){
        return fSw_n;
    }
    
    /** @brief Set the rock porosity */
    void SetPorosity(STATE phi){
        fporosity = phi;
    }
    
    /** @brief Get the rock porosity */
    STATE GetPorosity(){
        return fporosity;
    }
    
    //@}
    
  
    /**
     * @defgroup Write, Read and Print methods
     * @{
     */
    
    void Write(TPZStream &buf, int withclassid)
    {
        buf.Write(&fPressure_n);
        buf.Write(&fPressure);
    }

    void Read(TPZStream &buf, void *context)
    {
        buf.Read(&fPressure_n);
        buf.Read(&fPressure);
    }

    void Print(std::ostream &out) const
    {
        out << fPressure_n;
        out << fPressure;
    }
    
    //@}
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TRMMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* defined(__PZ__TRMMemory__) */
