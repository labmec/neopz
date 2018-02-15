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
     * @defgroup Required Memory quantities
     * @{
     */
    
    /** @brief Total flux */
    TPZManVector<REAL,3> fu;
    
    /** @brief Total flux at the previous timestep */
    TPZManVector<REAL,3> fu_n;
    
    /** @brief Weighted pressure */
    REAL fp;
    
    /** @brief Weighted pressure at the previous timestep */
    REAL fp_n;
    
    /** @brief Average weighted pressure */
    REAL fp_avg;
    
    /** @brief Average weighted pressure at the previous timestep */
    REAL fp_avg_n;
    
    /** @brief Alpha Saturation */
    REAL fsa;
    
    /** @brief Alpha saturation at the previous timestep */
    REAL fsa_n;
    
    /** @brief Beta Saturation */
    REAL fsb;
    
    /** @brief Beta saturation at the previous timestep */
    REAL fsb_n;
    
    //@}
    
    /**
     * @defgroup Memory quantities being transfered, tentative use
     * @{
     */
    
    /** @brief Total flux divergence */
    REAL fdivu;
    
    /** @brief Total flux divergence at the previous timestep */
    REAL fdivu_n;
    
    /** @brief Rock Porosity */
    REAL fporosity;

    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> fK;
    
    /** @brief Integration weight */
    REAL fw;
    
    /** @brief Jacobian det */
    REAL fdet;
    
    /** @brief Right hand side */
    REAL frhs;
    
    /** @brief Spatial coordinate */
    TPZManVector<REAL,3> fx;

    
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
        fp = copy.fp;
        fp_n = copy.fp_n;
        fp_avg = copy.fp_avg;
        fp_avg_n = copy.fp_avg_n;
        fsa = copy.fsa;
        fsa_n = copy.fsa_n;
        fsb = copy.fsb;
        fsb_n = copy.fsb_n;
        fporosity = copy.fporosity;
        fK = copy.fK;
        fw = copy.fw;
        fdet = copy.fdet;
        frhs = copy.frhs;
        fx = copy.fx;
    }
    
    TRMMemory &operator=(const TRMMemory &cp)
    {

        fp = cp.fp;
        fp_n = cp.fp_n;
        fp_avg = cp.fp_avg;
        fp_avg_n = cp.fp_avg_n;
        fsa = cp.fsa;
        fsa_n = cp.fsa_n;
        fsb = cp.fsb;
        fsb_n = cp.fsb_n;

        fu = cp.fu;
        fu_n = cp.fu_n;
        fdivu = cp.fdivu;
        fdivu_n = cp.fdivu_n;
        fporosity = cp.fporosity;
        fK = cp.fK;
        fw = cp.fw;
        fdet = cp.fdet;
        frhs = cp.frhs;
        fx = cp.fx;
        
        return *this;
    }

    void UpdateSolutionMemory()///// @omar:: I think this method is completely useless!!
    {
        //update saturation and pressure and total flux (un = unp1)
        DebugStop();
    }
    
    /**
     * @defgroup Set and Get methods
     * @{
     */
    
    /** @brief Set total flux */
    void Set_u(TPZManVector<REAL,3> &u){
        fu = u;
    }
    
    /** @brief Get total flux */
    TPZManVector<REAL,3> u(){
        return fu;
    }

    /** @brief Set total flux at last step */
    void Set_u_n(TPZManVector<REAL,3> &u_n){
        fu_n = u_n;
    }
    
    /** @brief Get total flux at last step */
    TPZManVector<REAL,3> u_n(){
        return fu_n;
    }
    
    /** @brief Set the weighted pressure */
    void Set_p(REAL p){
        fp = p;
    }
    
    /** @brief Get the weighted pressure */
    REAL p(){
        return fp;
    }
    
    /** @brief Set the weighted pressure at the previous timestep */
    void Set_p_n(REAL p_n){
        fp_n = p_n;
    }
    
    /** @brief Get the weighted pressure at the previous timestep */
    REAL p_n(){
        return fp_n;
    }
    
    /** @brief Set the average weighted pressure */
    void Set_p_avg(REAL p_avg){
        fp_avg = p_avg;
    }
    
    /** @brief Get the average weighted pressure */
    REAL p_avg(){
        return fp_avg;
    }
    
    /** @brief Set the average weighted pressure at the previous timestep */
    void Set_p_avg_n(REAL p_avg_n){
        fp_avg_n = p_avg_n;
    }
    
    /** @brief Get the average weighted pressure at the previous timestep */
    REAL p_avg_n(){
        return fp_avg_n;
    }
    
    /** @brief Set alpha saturation */
    void Set_sa(REAL sa){
        fsa = sa;
    }
    
    /** @brief Get alpha saturation */
    REAL sa(){
        return fsa;
    }
    
    /** @brief Set alpha saturation at last step */
    void Set_sa_n(REAL sa_n){
        fsa_n = sa_n;
    }
    
    /** @brief Get alpha saturation at last step */
    REAL sa_n(){
        return fsa_n;
    }
    
    /** @brief Set beta saturation */
    void Set_sb(REAL sb){
        fsb = sb;
    }
    
    /** @brief Get beta saturation */
    REAL sb(){
        return fsb;
    }
    
    /** @brief Set beta saturation at last step */
    void Set_sb_n(REAL sb_n){
        fsb_n = sb_n;
    }
    
    /** @brief Get beta saturation at last step */
    REAL sb_n(){
        return fsb_n;
    }
    
    //@}
    
    
    
    
  
    /**
     * @defgroup Write, Read and Print methods
     * @{
     */
    
    void Write(TPZStream &buf, int withclassid) const{
//        buf.Write(&fPressure_n);
//        buf.Write(&fPressure);
        DebugStop();
    }

    void Read(TPZStream &buf, void *context)
    {
//        buf.Read(&fPressure_n);
//        buf.Read(&fPressure);
        DebugStop();
    }

    void Print(std::ostream &out) const
    {
//        out << fPressure_n;
//        out << fPressure;
        DebugStop();
    }
    
    //@}
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TRMMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* defined(__PZ__TRMMemory__) */
