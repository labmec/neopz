//
//  TRMPhaseInterfaceMemory.h
//  PZ
//
//  Created by Philippe Devloo on 7/5/15.
//
//

#ifndef __PZ__TRMPhaseMemory__
#define __PZ__TRMPhaseMemory__

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"
#include "pzfmatrix.h"


class TRMPhaseMemory
{

private:
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Spatial memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief first lame undrained parameter */
    REAL f_lambda;
    
    /** @brief second lame undrained parameter */
    REAL f_mu;
    
    /** @brief Storage constrained modulus */
    REAL f_S_e;
    
    /** @brief Maurice Biot coefficient */
    REAL f_alpha;
    
    /** @brief contains the normal flux per surface area */
    REAL fun;
    
    /** @brief contains the volumetric average pressure */
    REAL fp_avg;

    /** @brief contains the volumetric average pressure at last time step */
    REAL fp_avg_n;
    
    /** @brief contains the volumetric saturation of alpha phase at initial state */
    REAL fsa_0;
    
    /** @brief contains the volumetric saturation of beta phase at initial state */
    REAL fsb_0;
    
    /** @brief contains the volumetric saturation of alpha phase */
    REAL fsa;
    
    /** @brief contains the volumetric saturation of alhpa at last time step */
    REAL fsa_n;
    
    /** @brief contains the volumetric saturation of beta phase */
    REAL fsb;

    /** @brief contains the volumetric saturation of alhpa at last time step */
    REAL fsb_n;
    
    /** @brief Rock Porosity */
    REAL fporosity;
    
    /** @brief Absolute permeability */
    TPZFNMatrix<9,REAL> fK;
    
    /** @brief Absolute permeability inverse */
    TPZFNMatrix<9,REAL> fKinv;
    
    /** @brief Spatial coordinate */
    TPZManVector<REAL,3> fx;
    
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

    /** @brief weighted pressure at intial state */
    REAL f_p_p_0;
    
public:

    
    /** @brief Default constructor */
    TRMPhaseMemory();
    
    /** @brief Default destructor */
    ~TRMPhaseMemory();
    
    TRMPhaseMemory(const TRMPhaseMemory &copy)
    {
        fp_avg      = copy.fp_avg;
        fp_avg_n    = copy.fp_avg_n;
        fsa_0         = copy.fsa_0;
        fsb_0         = copy.fsb_0;
        fsa         = copy.fsa;
        fsa_n       = copy.fsa_n;
        fsb         = copy.fsb;
        fsb_n       = copy.fsb_n;
        
        fporosity = copy.fporosity;
        fK = copy.fK;
        fKinv = copy.fKinv;
        fx = copy.fx;
        
    }
    
    TRMPhaseMemory &operator=(const TRMPhaseMemory &other)
    {
        
        if (this != & other) // prevent self-assignment
        {
            fp_avg      = other.fp_avg;
            fp_avg_n    = other.fp_avg_n;
            fsa_0       = other.fsa_0;
            fsb_0       = other.fsb_0;
            fsa         = other.fsa;
            fsa_n       = other.fsa_n;
            fsb         = other.fsb;
            fsb_n       = other.fsb_n;
            
            fporosity   = other.fporosity;
            fK          = other.fK;
            fKinv       = other.fKinv;
            fx          = other.fx;
            
        }
        return *this;
        
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Segregated Memory (\partial Gamma and Omega) :: Spatial memory items
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Set first lame undrained parameter */
    void Set_lambda(REAL lambda){
        f_lambda = lambda;
    }
    
    /** @brief Get first lame undrained parameter */
    REAL lambda(){
        return f_lambda;
    }
    
    /** @brief Set second lame undrained parameter */
    void Set_mu(REAL mu){
        f_mu = mu;
    }
    
    /** @brief Get second lame undrained parameter */
    REAL mu(){
        return f_mu;
    }
    
    /** @brief Set storage constrained modulus */
    void Set_S_e(REAL S_e){
        f_S_e = S_e;
    }
    
    /** @brief Get storage constrained modulus */
    REAL S_e(){
        return f_S_e;
    }
    
    /** @brief Set Maurice Biot coefficient */
    void Set_alpha(REAL alpha){
        f_alpha = alpha;
    }
    
    /** @brief Get Maurice Biot coefficient */
    REAL alpha(){
        return f_alpha;
    }
    
   
    /** @brief Set the average normal flux */
    void Set_un(REAL un){
        fun = un;
    }
    
    /** @brief Set the average normal flux */
    REAL un(){
        return fun;
    }
    
    /** @brief Set the average weighted pressure */
    void Set_p_avg(REAL p_avg){
        fp_avg = p_avg;
    }
    
    /** @brief Get the average weighted pressure */
    REAL & p_avg(){
        return fp_avg;
    }
    
    /** @brief Set the average weighted pressure at the previous timestep */
    void Set_p_avg_n(REAL p_avg_n){
        fp_avg_n = p_avg_n;
    }
    
    /** @brief Get the average weighted pressure at the previous timestep */
    REAL & p_avg_n(){
        return fp_avg_n;
    }
    
    /** @brief Set alpha saturation at initial state */
    void Set_sa_0(REAL sa_0){
        fsa_0 = sa_0;
    }
    
    /** @brief Get alpha saturation at initial state */
    REAL & sa_0(){
        return fsa_0;
    }
    
    /** @brief Set beta saturation at initial state */
    void Set_sb_0(REAL sb_0){
        fsb_0 = sb_0;
    }
    
    /** @brief Get beta saturation at initial state */
    REAL & sb_0(){
        return fsb_0;
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
    
    /** @brief Set intact rock Porosity */
    void Set_phi_0(REAL porosity){
        fporosity = porosity;
    }
    
    /** @brief Get intact rock Porosity */
    REAL phi_0(){
        return fporosity;
    }
    
    /** @brief Set intact absolute permeability */
    void Set_K_0(TPZFNMatrix<9,REAL> K){
        fK = K;
    }
    
    /** @brief Get intact absolute permeability */
    TPZFNMatrix<9,REAL> & K_0(){
        return fK;
    }
    
    /** @brief Set intact absolute permeability inv */
    void Set_Kinv_0(TPZFNMatrix<9,REAL> Kinv){
        fKinv = Kinv;
    }
    
    /** @brief Get intact absolute permeability inv */
    TPZFNMatrix<9,REAL> & Kinv_0(){
        return fKinv;
    }
    
    /** @brief Set x coordinate */
    void Set_x(TPZManVector<REAL,3> x){
        fx = x;
    }
    
    
    ///////// elliptic data
    
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
    
    /** @brief set grad_u at last time */
    void Set_grad_u(TPZFMatrix<STATE> & grad_u){
        f_e_grad_u = grad_u;
    }
    
    /** @brief get grad_u at last time */
    TPZFMatrix<STATE> & grad_u(){
        return f_e_grad_u;
    }
    
    // current time state items
    
    /** @brief set grad_u at current time */
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n){
        f_e_grad_u_n = grad_u_n;
    }
    
    /** @brief get grad_u at current time */
    TPZFMatrix<STATE> & grad_u_n(){
        return f_e_grad_u_n;
    }
    
    
    /** @brief set weighted pressure at intial state */
    void Set_p_0(STATE & p_p_0){
        f_p_p_0 = p_p_0;
    }
    
    /** @brief get weighted pressure at intial state */
    STATE & p_0(){
        return f_p_p_0;
    }
    
    
    /** @brief Get x coordinate */
    TPZManVector<REAL,3> & x(){
        return fx;
    }

    void Write(TPZStream &buf, int withclassid)
    {
        buf.Write(&fp_avg);
        buf.Write(&fp_avg_n);
        buf.Write(&fsa);
        buf.Write(&fsa_n);
        buf.Write(&fsb);
        buf.Write(&fsb_n);
    }

    void Read(TPZStream &buf, void *context)
    {
        buf.Read(&fp_avg);
        buf.Read(&fp_avg_n);
        buf.Read(&fsa);
        buf.Read(&fsa_n);
        buf.Read(&fsb);
        buf.Read(&fsb_n);
    }

    void Print(std::ostream &out) const
    {
        
        out << "TRMPhaseMemory item, with values ";
        out << fp_avg;
        out << fp_avg_n;
        out << fsa;
        out << fsa_n;
        out << fsb;
        out << fsb_n;
    }


};

inline std::ostream &operator<<(std::ostream &out,const TRMPhaseMemory &mem)
{
    mem.Print(out);
    return out;
}



#endif /* defined(__PZ__TRMPhaseInterfaceMemory__) */