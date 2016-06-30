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
    TPZManVector<STATE,3> fu;

    /** @brief Total flux at the previous timestep */
    TPZManVector<STATE,3> fu_n;

    /** @brief Total flux divergence */
    STATE fdivu;
    
    /** @brief Total flux divergence at the previous timestep */
    STATE fdivu_n;
    
    /** @brief Weighted Pressure */
    STATE fPressure;
    
    /** @brief Weighted Pressure at the previous timestep */
    STATE fPressure_n;
    
    /** @brief Average Weighted Pressure */
    STATE fA_Pressure;
    
    /** @brief Average Weighted Pressure at the previous timestep */
    STATE fA_Pressure_n;
    
    /** @brief Water Saturation */
    STATE fSw;
    
    /** @brief Water saturation at the previous timestep */
    STATE fSw_n;
    
    /** @brief Rock Porosity */
    STATE fporosity;

    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> fK;
    
    /** @brief Integration weight */
    STATE fw;
    
    /** @brief Jacobian det */
    STATE fdet;
    
    /** @brief Right hand side */
    STATE frhs;
    
    /** @brief Spatial coordinate */
    TPZManVector<STATE,3> fx;

    
    //@}
    
public:
    
    /** @brief Default constructor */
    TRMMemory();
    
    /** @brief Default destructor */
    ~TRMMemory();
    
    TRMMemory(const TRMMemory &copy)
    {

        fu = copy.fu;
        fu_n = copy.fu_n;
        fdivu = copy.fdivu;
        fdivu_n = copy.fdivu_n;
        fPressure = copy.fPressure;
        fPressure_n = copy.fA_Pressure_n;
        fA_Pressure = copy.fA_Pressure;
        fA_Pressure_n = copy.fA_Pressure_n;
        fSw = copy.fSw;
        fSw_n = copy.fSw_n;
        fporosity = copy.fporosity;
        fK = copy.fK;
        fw = copy.fw;
        fdet = copy.fdet;
        frhs = copy.frhs;
        fx = copy.fx;
    }
    
    TRMMemory &operator=(const TRMMemory &cp)
    {
        /** @brief Total flux */
        fu = cp.fu;
        
        /** @brief Total flux at the previous timestep */
        fu_n = cp.fu_n;
        
        /** @brief Total flux divergence */
        fdivu = cp.fdivu;
        
        /** @brief Total flux divergence at the previous timestep */
        fdivu_n = cp.fdivu_n;
        
        /** @brief Weighted Pressure */
        fPressure = cp.fPressure;
        
        /** @brief Weighted Pressure at the previous timestep */
        fPressure_n = cp.fPressure_n;
        
        /** @brief Water Saturation */
        fSw = cp.fSw;
        
        /** @brief Water saturation at the previous timestep */
        fSw_n = cp.fSw_n;
        
        /** @brief Rock Porosity */
        fporosity = cp.fporosity;
        
        /** @brief Absolute permeability */
        fK = cp.fK;
        
        /** @brief Integration weight */
        fw = cp.fw;
        
        /** @brief Jacobian det */
        fdet = cp.fdet;
        
        /** @brief Right hand side */
        frhs = cp.frhs;
        
        /** @brief Spatial coordinate */
        fx = cp.fx;
        
        return *this;
    }

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
    void SetTotal_Flux(TPZManVector<STATE> &u){
        fu = u;
    }
    
    /** @brief Get the total flux */
    TPZManVector<STATE> GetTotal_Flux(){
        return fu;
    }
    
    /** @brief Set the total flux at the previous timestep */
    void SetTotal_Flux_n(TPZManVector<STATE> &u_n){
        fu_n = u_n;
    }
    
    /** @brief Get the total flux at the previous timestep */
    TPZManVector<STATE> GetTotal_Flux_n(){
        return fu_n;
    }
    
    /** @brief Set the total flux at the previous timestep */
    void SetDiv_Flux(STATE divu){
        fdivu = divu;
    }
    
    /** @brief Get the total flux at the previous timestep */
    STATE GetDiv_Flux(){
        return fdivu;
    }
    
    /** @brief Set the total flux at the previous timestep */
    void SetDiv_Flux_n(STATE divu_n){
        fdivu_n = divu_n;
    }
    
    /** @brief Get the total flux at the previous timestep */
    STATE GetDiv_Flux_n(){
        return fdivu_n;
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
    
    /** @brief Set integration weight */
    void SetWeight(STATE w){
        fw = w;
    }
    
    /** @brief Get integration weight */
    STATE GetWeight(){
        return fw;
    }
    
    /** @brief Set Jacobian det */
    void SetDetJac(STATE det){
        fdet = det;
    }
    
    /** @brief Get Jacobian det  */
    STATE GetDetJac(){
        return fdet;
    }
    
    /** @brief Set Right hand side */
    void SetRhs(STATE rhs){
        frhs = rhs;
    }
    
    /** @brief Get Right hand side */
    STATE GetRhs(){
        return frhs;
    }
    
    /** @brief Set Spatial coordinate */
    void SetX(TPZManVector<STATE>  x){
        fx = x;
    }
    
    /** @brief Get Spatial coordinate */
    TPZManVector<STATE>  GetX(){
        return fx;
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
