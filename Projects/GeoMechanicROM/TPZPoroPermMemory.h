//
//  TPZPoroPermMemory.h
//  PZ
//
//  Created by Omar on 9/6/16.
//
//

#ifndef TPZPoroPermMemory_h
#define TPZPoroPermMemory_h

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"
#include "pzfmatrix.h"

class TPZPoroPermMemory {
   
    /** @brief RB functions */
    TPZFMatrix<STATE> fphi_u;
    
    /** @brief RB functions */
    TPZFMatrix<STATE> fgrad_phi_u;
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> fu_n;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> fgrad_u_n;
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> fu;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> fgrad_u;
    
    
public:
    
    /** @brief Default constructor */
    TPZPoroPermMemory();
    
    /** @brief Default destructor */
    ~TPZPoroPermMemory();
    
    TPZPoroPermMemory(const TPZPoroPermMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;

    }
    
    TPZPoroPermMemory &operator=(const TPZPoroPermMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;
        
        return *this;
    }
    
    /** @brief Set RB i function */
    void Set_phi_u(TPZFMatrix<STATE> & phi_u){
        fphi_u = phi_u;
    }
    
    /** @brief Get RB functions  */
    TPZFMatrix<STATE> & phi_u(){
        return fphi_u;
    }
    
    /** @brief Set RB i function */
    void Set_grad_phi_u(TPZFMatrix<STATE>& grad_phi_u){
        fgrad_phi_u = grad_phi_u;
    }
    
    /** @brief Get RB functions  */
    TPZFMatrix<STATE> & grad_phi_u(){
        return fgrad_phi_u;
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

inline std::ostream &operator<<(std::ostream &out,const TPZPoroPermMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPoroPermMemory_h */
