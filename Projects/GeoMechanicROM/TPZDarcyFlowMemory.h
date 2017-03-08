//
//  TPZDarcyFlowMemory.h
//  PZ
//
//  Created by omar on 3/7/17.
//
//

#ifndef TPZDarcyFlowMemory_h
#define TPZDarcyFlowMemory_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"
#include "pzfmatrix.h"

class TPZDarcyFlowMemory {
    
    /** @brief flux space functions */
    TPZFMatrix<STATE> fphi_q;
    
    /** @brief flux space div functions */
    TPZFMatrix<STATE> fdiv_phi_q;
    
    /** @brief pressure space functions */
    TPZFMatrix<STATE> fphi_p;
    
    /** @brief pressure space functions */
    TPZFMatrix<STATE> fgrad_phi_p;
    
    /** @brief flux */
    TPZFNMatrix<3,REAL> fq_n;
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> fu_n;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> fgrad_u_n;
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> fu;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> fgrad_u;
    
    /** @brief gradient of u_n */
    REAL fporosity_0;
    
    /** @brief pressure at initial state */
    REAL fp_0;
    
    /** @brief gradient of u_n at intial state*/
    TPZFNMatrix<9,REAL> fgrad_u_0;
    
    /** @brief pressure at current state */
    REAL fp_n;
    
    /** @brief pressure at current state */
    REAL fp;
    
    /** @brief gradient of u_n at intial state*/
    TPZFNMatrix<9,REAL> fgrad_p_n;
    
    /** @brief gradient of u_n at intial state*/
    TPZFNMatrix<9,REAL> fgrad_p;
    
public:
    
    /** @brief Default constructor */
    TPZDarcyFlowMemory();
    
    /** @brief Default destructor */
    ~TPZDarcyFlowMemory();
    
    
    TPZDarcyFlowMemory(const TPZDarcyFlowMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;
        
    }
    
    TPZDarcyFlowMemory &operator=(const TPZDarcyFlowMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;
        
        return * this;
    }
    
    /** @brief Set flux space */
    void Set_phi_q(TPZFMatrix<STATE> & phi_q){
        fphi_q = phi_q;
    }
    
    /** @brief Set flux space  */
    TPZFMatrix<STATE> & phi_q(){
        return fphi_q;
    }
    
    /** @brief Set pressure space */
    void Set_phi_p(TPZFMatrix<STATE> & phi_p){
        fphi_p = phi_p;
    }
    
    /** @brief Get pressure space  */
    TPZFMatrix<STATE> & phi_p(){
        return fphi_p;
    }
    
    /** @brief Set pressure space */
    void Set_grad_phi_p(TPZFMatrix<STATE> & grad_phi_p){
        fgrad_phi_p = grad_phi_p;
    }
    
    /** @brief Get pressure space  */
    TPZFMatrix<STATE> & grad_phi_p(){
        return fgrad_phi_p;
    }
    
    /** @brief Set flux */
    void Set_q_n(TPZFMatrix<STATE> & q_n){
        fq_n = q_n;
    }
    
    /** @brief Get flux */
    TPZFMatrix<STATE> & q_n(){
        return fq_n;
    }
    
    /** @brief Set displacement at last state */
    void Set_u_n(TPZFNMatrix<3,REAL> &u_n){
        fu_n = u_n;
    }
    
    /** @brief Get displacement at last state */
    TPZFNMatrix<3,REAL> & u_n(){
        return fu_n;
    }
    
    /** @brief Set gradient of u_n */
    void Set_grad_u_n(TPZFNMatrix<9,REAL> & grad_u_n){
        fgrad_u_n = grad_u_n;
    }
    
    /** @brief Get gradient of u_n */
    TPZFNMatrix<9,REAL> & grad_u_n(){
        return fgrad_u_n;
    }
    
    /** @brief Set displacement at last state */
    void Set_u(TPZFNMatrix<3,REAL> &u){
        fu = u;
    }
    
    /** @brief Get displacement at last state */
    TPZFNMatrix<3,REAL> & u(){
        return fu;
    }
    
    /** @brief Set gradient of u */
    void Set_grad_u(TPZFNMatrix<9,REAL> & grad_u){
        fgrad_u = grad_u;
    }
    
    /** @brief Get gradient of u */
    TPZFNMatrix<9,REAL> & grad_u(){
        return fgrad_u;
    }
    
    /** @brief Set gradient of u */
    void Set_grad_u_0(TPZFNMatrix<9,REAL> & grad_u_0){
        fgrad_u_0 = grad_u_0;
    }
    
    /** @brief Get gradient of u */
    TPZFNMatrix<9,REAL> & grad_u_0(){
        return fgrad_u_0;
    }
    
    /** @brief Set porosity */
    void Set_porosity_0(REAL & porosity_0){
        fporosity_0 = porosity_0;
    }
    
    /** @brief Get porosity */
    REAL & porosity_0(){
        return fporosity_0;
    }
    
    /** @brief Set pressure */
    void Set_pressure_0(REAL & p_0){
        fp_0 = p_0;
    }
    
    /** @brief Get pressure */
    REAL & p_0(){
        return fp_0;
    }
    
    /** @brief Set pressure */
    void Set_p(REAL & p){
        fp = p;
    }
    
    /** @brief Get pressure */
    REAL & p(){
        return fp;
    }
    
    /** @brief Set pressure */
    void Set_p_n(REAL & p_n){
        fp_n = p_n;
    }
    
    /** @brief Get pressure */
    REAL & p_n(){
        return fp_n;
    }
    
    /** @brief Set gradient of p */
    void Set_grad_p_n(TPZFNMatrix<9,REAL> & grad_p_n){
        fgrad_p_n = grad_p_n;
    }
    
    /** @brief Get gradient of p_n */
    TPZFNMatrix<9,REAL> & grad_p_n(){
        return fgrad_p_n;
    }
    
    /** @brief Set gradient of p */
    void Set_grad_p(TPZFNMatrix<9,REAL> & grad_p){
        fgrad_p = grad_p;
    }
    
    /** @brief Get gradient of p */
    TPZFNMatrix<9,REAL> & grad_p(){
        return fgrad_p;
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
        out << "Fancy print!" << std::endl;
    }
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZDarcyFlowMemory &mem)
{
    mem.Print(out);
    return out;
}



#endif /* TPZDarcyFlowMemory_h */
