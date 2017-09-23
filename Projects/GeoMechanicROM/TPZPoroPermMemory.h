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
    
    
    /** @brief displacements */
    TPZManVector<REAL,3> fu_n;
    
    /** @brief gradient of u_n */
    TPZFMatrix<REAL> fgrad_u_n;
    
   /** @brief elastic strain at n */
    TPZFMatrix<REAL> fepsilon_e_n;

   /** @brief plastic strain at n */
    TPZFMatrix<REAL> fepsilon_p_n;
    
public:
    
    /** @brief Default constructor */
    TPZPoroPermMemory();
    
    /** @brief Default destructor */
    ~TPZPoroPermMemory();
    
    TPZPoroPermMemory(const TPZPoroPermMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;
        fepsilon_e_n = copy.fepsilon_e_n;
        fepsilon_p_n = copy.fepsilon_p_n;

    }
    
    TPZPoroPermMemory &operator=(const TPZPoroPermMemory &copy)
    {
        
        fu_n = copy.fu_n;
        fgrad_u_n = copy.fgrad_u_n;
        fepsilon_e_n = copy.fepsilon_e_n;
        fepsilon_p_n = copy.fepsilon_p_n;
        
        return *this;
    }
    
    void UpdateSolutionMemory()///// @omar:: I think this method is completely useless!!
    {
        //update saturation and pressure and total flux (un = unp1)
        DebugStop();
    }
    
    /** @brief Set displacement at last state */
    void Set_u_n(TPZManVector<REAL,3> &u_n){
        fu_n = u_n;
    }
    
    /** @brief Get displacement at last state */
    TPZManVector<REAL,3> u_n(){
        return fu_n;
    }
    
    /** @brief Set gradient of u_n */
    void Set_grad_u_n(TPZFMatrix<REAL> & grad_u_n){
        fgrad_u_n = grad_u_n;
    }
    
    /** @brief Get gradient of u_n */
    TPZFMatrix<REAL> grad_u_n(){
        return fgrad_u_n;
    }
    
    /** @brief Set elastic strain at n */
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n){
        fepsilon_e_n = epsilon_e_n;
    }
    
    /** @brief Get elastic strain at n */
    TPZFMatrix<REAL> epsilon_e_n(){
        return fepsilon_e_n;
    }
    
    /** @brief Set plastic strain at n */
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n){
        fepsilon_p_n = epsilon_p_n;
    }
    
    /** @brief Get plastic strain at n */
    TPZFMatrix<REAL> epsilon_p_n(){
        return fepsilon_p_n;
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
        //        out << fPressure_n;
        //        out << fPressure;
        DebugStop();
    }
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZPoroPermMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPoroPermMemory_h */
