//
//  TPZPoroPermMemory.h
//  PZ
//
//  Created by Omar and Manouchehr on 9/6/16.
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
    TPZManVector<REAL,3> m_u_n;
    
    /** @brief gradient of u_n */
    TPZFMatrix<REAL> m_grad_u_n;
    
   /** @brief elastic strain at n */
    TPZFMatrix<REAL> m_epsilon_e_n;

   /** @brief plastic strain at n */
    TPZFMatrix<REAL> m_epsilon_p_n;
    
public:
    
    /** @brief Default constructor */
    TPZPoroPermMemory();
    
    /** @brief Default destructor */
    ~TPZPoroPermMemory();
    
    TPZPoroPermMemory(const TPZPoroPermMemory &copy)
    {
        
        m_u_n = copy.m_u_n;
        m_grad_u_n = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;

    }
    
    TPZPoroPermMemory &operator=(const TPZPoroPermMemory &copy)
    {
        
        m_u_n = copy.m_u_n;
        m_grad_u_n = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;
        
        return *this;
    }
    
    void UpdateSolutionMemory()///// @omar:: I think this method is completely useless!!
    {
        //update saturation and pressure and total flux (un = unp1)
        DebugStop();
    }
    
    /** @brief Set displacement at last state */
    void Set_u_n(TPZManVector<REAL,3> &u_n){
        m_u_n = u_n;
    }
    
    /** @brief Get displacement at last state */
    TPZManVector<REAL,3> u_n(){
        return m_u_n;
    }
    
    /** @brief Set gradient of u_n */
    void Set_grad_u_n(TPZFMatrix<REAL> & grad_u_n){
        m_grad_u_n = grad_u_n;
    }
    
    /** @brief Get gradient of u_n */
    TPZFMatrix<REAL> grad_u_n(){
        return m_grad_u_n;
    }
    
    /** @brief Set elastic strain at n */
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n){
        m_epsilon_e_n = epsilon_e_n;
    }
    
    /** @brief Get elastic strain at n */
    TPZFMatrix<REAL> epsilon_e_n(){
        return m_epsilon_e_n;
    }
    
    /** @brief Set plastic strain at n */
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n){
        m_epsilon_p_n = epsilon_p_n;
    }
    
    /** @brief Get plastic strain at n */
    TPZFMatrix<REAL> epsilon_p_n(){
        return m_epsilon_p_n;
    }
    
    void Write(TPZStream &buf, int withclassid)
    {
//                buf.Write(&m_Pressure_n);
//                buf.Write(&m_Pressure);
        DebugStop();
    }
    
    void Read(TPZStream &buf, void *context)
    {
        //        buf.Read(&m_Pressure_n);
        //        buf.Read(&m_Pressure);
        DebugStop();
    }
    
    void Print(std::ostream &out) const
    {
        //        out << m_Pressure_n;
        //        out << m_Pressure;
        DebugStop();
    }
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZPoroPermMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPoroPermMemory_h */
