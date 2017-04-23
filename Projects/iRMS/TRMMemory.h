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
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Elliptic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief elliptic functions functions */
    TPZFMatrix<STATE> f_e_phi_u;
    
    /** @brief elliptic functions functions */
    TPZFMatrix<STATE> f_e_grad_phi_u;
    
    // initial state items
    
    /** @brief gradient of u_n at intial state*/
    TPZFNMatrix<9,REAL> f_e_grad_u_0;
    
    /** @brief sigma at intial state*/
    TPZFNMatrix<9,REAL> f_e_sigma_0;
    
    // last time state items
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> f_e_u;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> f_e_grad_u;
    
    // current time state items
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> f_e_u_n;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> f_e_grad_u_n;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Parabolic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief parabolic q base functions */
    TPZFMatrix<STATE> f_p_phi_q;
    
    /** @brief parabolic div_q base functions */
    TPZFMatrix<STATE> f_p_div_phi_q;
    
    /** @brief parabolic p base functions */
    TPZFMatrix<STATE> f_p_phi_p;
    
    // initial state items
    
    /** @brief weighted pressure at intial state */
    REAL f_p_p_0;
    
    /** @brief Rock Porosity */
    REAL f_p_phi_0;
    
    /** @brief absolute permeability */
    TPZFNMatrix<9,REAL> f_p_K_0;
    
    /** @brief absolute permeability inverse */
    TPZFNMatrix<9,REAL> f_p_Kinv_0;
    
    // last time state items
    
    /** @brief total velocity at last state */
    TPZManVector<REAL,3> f_p_q;
    
    /** @brief divergence of velocity at last state */
    REAL f_p_div_q;
    
    /** @brief weighted pressure at last state */
    REAL f_p_p;
    
    // current time state items
    
    /** @brief total velocity */
    TPZManVector<REAL,3> f_p_q_n;
    
    /** @brief divergence of velocity at the current timestep */
    REAL f_p_div_q_n;
    
    /** @brief weighted pressure at the current timestep */
    REAL f_p_p_n;
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
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

    /** @brief Total flux divergence */
    REAL fdivu;
    
    /** @brief Total flux divergence at the previous timestep */
    REAL fdivu_n;
    
    /** @brief Integration weight */
    REAL fw;
    
    /** @brief Jacobian det */
    REAL fdet;
    
    /** @brief Right hand side */
    REAL frhs;
    
    /** @brief Spatial coordinate */
    TPZManVector<REAL,3> fx;

    
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
        f_p_phi_0 = copy.f_p_phi_0;
        f_p_K_0 = copy.f_p_K_0;
        f_p_Kinv_0 = copy.f_p_Kinv_0;
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
        f_p_phi_0 = cp.f_p_phi_0;
        f_p_K_0 = cp.f_p_K_0;
        f_p_Kinv_0 = cp.f_p_Kinv_0;
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
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Elliptic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief set elliptic u base functions */
    void Set_phi_u(TPZFMatrix<STATE> & phi_u){
        f_e_phi_u = phi_u;
    }

    /** @brief get elliptic u base functions */
    TPZFMatrix<STATE> & phi_u(){
        return f_e_phi_u;
    }
    
    /** @brief set elliptic grad_u base functions */
    void Set_grad_phi_u(TPZFMatrix<STATE> & grad_phi_u){
        f_e_grad_phi_u = grad_phi_u;
    }
    
    /** @brief get elliptic functions functions */
    TPZFMatrix<STATE> & grad_phi_u(){
        return f_e_grad_phi_u;
    }
    
    // initial state items
    
    /** @brief set gradient of u_n at intial state */
    void Set_grad_u_0(TPZFMatrix<STATE> & grad_u_0){
        f_e_grad_u_0 = grad_u_0;
    }
    
    /** @brief get gradient of u_n at intial state */
    TPZFMatrix<STATE> & grad_u_0(){
        return f_e_grad_u_0;
    }
    
    /** @brief set sigma at intial state */
    void Set_sigma_0(TPZFMatrix<STATE> & sigma_0){
        f_e_sigma_0 = sigma_0;
    }
    
    /** @brief get sigma at intial state */
    TPZFMatrix<STATE> & sigma_0(){
        return f_e_sigma_0;
    }

    
    // last time state items
    
    /** @brief set displacements at last time */
    void Set_u(TPZFMatrix<STATE> & u){
        f_e_u = u;
    }
    
    /** @brief get displacements at last time */
    TPZFMatrix<STATE> & u(){
        return f_e_u;
    }
    
    /** @brief set grad_u at last time */
    void Set_grad_u(TPZFMatrix<STATE> & grad_u){
        f_e_grad_u = grad_u;
    }
    
    /** @brief get grad_u at last time */
    TPZFMatrix<STATE> & grad_u(){
        return f_e_grad_u;
    }
    
    // current time state items
    
    /** @brief set displacements at current time */
    void Set_u_n(TPZFMatrix<STATE> & u_n){
        f_e_u_n = u_n;
    }
    
    /** @brief get displacements at current time */
    TPZFMatrix<STATE> & u_n(){
        return f_e_u_n;
    }
    
    /** @brief set grad_u at current time */
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n){
        f_e_grad_u_n = grad_u_n;
    }
    
    /** @brief get grad_u at current time */
    TPZFMatrix<STATE> & grad_u_n(){
        return f_e_grad_u_n;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Elliptic memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    // Basis functions
    
    /** @brief set parabolic q base functions */
    void Set_phi_q(TPZFMatrix<STATE> & phi_q){
        f_p_phi_q = phi_q;
    }
    
    /** @brief get parabolic q base functions */
    TPZFMatrix<STATE> & phi_q(){
        return f_p_phi_q;
    }
    
    /** @brief set parabolic div_q base functions */
    void Set_div_phi_q(TPZFMatrix<STATE> & div_phi_q){
        f_p_div_phi_q = div_phi_q;
    }
    
    /** @brief get parabolic div_q base functions */
    TPZFMatrix<STATE> & div_phi_q(){
        return f_p_div_phi_q;
    }
    
    /** @brief set parabolic p base functions */
    void Set_phi_p(TPZFMatrix<STATE> & phi_p){
        f_p_phi_p = phi_p;
    }
    
    /** @brief get parabolic p base functions */
    TPZFMatrix<STATE> & phi_p(){
        return f_p_phi_p;
    }
    
    
    // initial state items
    
    /** @brief Set intact rock Porosity */
    void Set_phi_0(REAL phi_0){
        f_p_phi_0 = phi_0;
    }
    
    /** @brief Get intact rock Porosity */
    REAL phi_0(){
        return f_p_phi_0;
    }
    
    /** @brief Set intact absolute permeability */
    void Set_K_0(TPZFNMatrix<9,REAL> K_0){
        f_p_K_0 = K_0;
    }
    
    /** @brief Get intact absolute permeability */
    TPZFNMatrix<9,REAL> & K_0(){
        return f_p_K_0;
    }
    
    /** @brief Set intact absolute permeability inv */
    void Set_Kinv_0(TPZFNMatrix<9,REAL> Kinv_0){
        f_p_Kinv_0 = Kinv_0;
    }
    
    /** @brief Get intact absolute permeability inv */
    TPZFNMatrix<9,REAL> & Kinv_0(){
        return f_p_Kinv_0;
    }
    
    /** @brief set weighted pressure at intial state */
    void Set_p_0(STATE & p_p_0){
        f_p_p_0 = p_p_0;
    }
    
    /** @brief get weighted pressure at intial state */
    STATE & p_0(){
        return f_p_p_0;
    }
    
    // last time state items

    /** @brief get total velocity at last state */
    void Set_q(TPZManVector<REAL,3> & q){
        f_p_q = q;
    }
    
    /** @brief get total velocity at last state */
    TPZManVector<REAL,3> & q(){
        return f_p_q;
    }
    
    /** @brief set divergence of velocity at last state */
    void Set_div_q(STATE & div_q){
        f_p_div_q = div_q;
    }
    
    /** @brief get divergence of velocity at last state */
    STATE & div_q(){
        return f_p_div_q;
    }
    
    /** @brief set weighted pressure at last state */
    void Set_p(STATE & p_p){
        f_p_p = p_p;
    }
    
    /** @brief get weighted pressure at last state */
    STATE & p(){
        return f_p_p;
    }
    
    // current time state items
    
    /** @brief get total velocity at current state */
    void Set_q_n(TPZManVector<REAL,3> & q_n){
        f_p_q_n = q_n;
    }
    
    /** @brief get total velocity at current state */
    TPZManVector<REAL,3> & q_n(){
        return f_p_q_n;
    }
    
    /** @brief set divergence of velocity at current state */
    void Set_div_q_n(STATE & div_q_n){
        f_p_div_q_n = div_q_n;
    }
    
    /** @brief get divergence of velocity at current state */
    STATE & div_q_n(){
        return f_p_div_q_n;
    }
    
    /** @brief set weighted pressure at current state */
    void Set_p_n(STATE & p_p_n){
        f_p_p_n = p_p_n;
    }
    
    /** @brief get weighted pressure at current state */
    STATE & p_n(){
        return f_p_p_n;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /** @brief Set total flux */
    void Set_u(TPZManVector<REAL,3> &u){
        fu = u;
    }
    
//    /** @brief Get total flux */
//    TPZManVector<REAL,3> u(){
//        return fu;
//    }

    /** @brief Set total flux at last step */
    void Set_u_n(TPZManVector<REAL,3> &u_n){
        fu_n = u_n;
    }
    
//    /** @brief Get total flux at last step */
//    TPZManVector<REAL,3> u_n(){
//        return fu_n;
//    }
    
//    /** @brief Set the weighted pressure */
//    void Set_p(REAL p){
//        fp = p;
//    }
    
//    /** @brief Get the weighted pressure */
//    REAL p(){
//        return fp;
//    }
    
//    /** @brief Set the weighted pressure at the previous timestep */
//    void Set_p_n(REAL p_n){
//        fp_n = p_n;
//    }
    
//    /** @brief Get the weighted pressure at the previous timestep */
//    REAL p_n(){
//        return fp_n;
//    }
    
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
    
    /** @brief Set x coordinate */
    void Set_x(TPZManVector<REAL,3> x){
        fx = x;
    }
    
    /** @brief Get x coordinate */
    TPZManVector<REAL,3> & x(){
        return fx;
    }
    
    
    void Write(TPZStream &buf, int withclassid)
    {
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
        out << fp_avg; // Dummy

    }
    
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TRMMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* defined(__PZ__TRMMemory__) */
