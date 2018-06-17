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
#include "pzfmatrix.h"

class TPZPoroPermMemory {
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Basis functions
    
    /** @brief elliptic functions functions */
    TPZFMatrix<STATE> m_phi_u;
    
    /** @brief elliptic functions functions */
    TPZFMatrix<STATE> m_grad_phi_u;
    
    
    // initial state items
    
    /** @brief gradient of u_n at intial state*/
    TPZFNMatrix<9,REAL> m_grad_u_0;
    
    /** @brief sigma at intial state*/
    TPZFNMatrix<9,REAL> m_sigma_0;
    
    
    
    // last time state items
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> m_u;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> m_grad_u;
    
    
    
    // current time state items
    
    /** @brief displacements */
    TPZFNMatrix<3,REAL> m_u_n;
    
    /** @brief gradient of u_n */
    TPZFNMatrix<9,REAL> m_grad_u_n;
    
    
   /** @brief elastic strain at n */
    TPZFMatrix<REAL> m_epsilon_e_n;

   /** @brief plastic strain at n */
    TPZFMatrix<REAL> m_epsilon_p_n;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /** @brief weighted pressure at intial state */
    REAL m_p_0;
    
    /** @brief weighted pressure at the current timestep */
    REAL m_p_n;
    
    
public:
    
    /** @brief Default constructor */
    TPZPoroPermMemory();
    
    /** @brief Default destructor */
    ~TPZPoroPermMemory();
    
    TPZPoroPermMemory(const TPZPoroPermMemory &copy)
    {
        
        m_u = copy.m_u;
        m_u_n = copy.m_u_n;
        m_grad_u = copy.m_grad_u;
        m_grad_u_n = copy.m_grad_u_n;
        m_epsilon_e_n = copy.m_epsilon_e_n;
        m_epsilon_p_n = copy.m_epsilon_p_n;

    }
    
    TPZPoroPermMemory &operator=(const TPZPoroPermMemory &copy)
    {
        
        m_u = copy.m_u;
        m_u_n = copy.m_u_n;
        m_grad_u = copy.m_grad_u;
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
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // initial state items
    
    /** @brief set gradient of u_n at intial state */
    void Set_grad_u_0(TPZFMatrix<STATE> & grad_u_0)
    {
        m_grad_u_0 = grad_u_0;
    }
    
    /** @brief get gradient of u_n at intial state */
    TPZFMatrix<STATE> & grad_u_0()
    {
        return m_grad_u_0;
    }
    

    
    // current time state items
    
    /** @brief set displacements at current time */
    void Set_u(TPZFMatrix<STATE> & u)
    {
        m_u = u;
    }
    
    /** @brief get displacements at current time */
    TPZFMatrix<STATE> & u(){
        return m_u;
    }
    
    
    /** @brief set grad_u at current time */
    void Set_grad_u(TPZFMatrix<STATE> & grad_u)
    {
        m_grad_u = grad_u;
    }
    
    /** @brief get grad_u at current time */
    TPZFMatrix<STATE> & grad_u(){
        return m_grad_u;
    }
    
    
    
    // last time state items
    
    /** @brief set displacements at last time */
    void Set_u_n(TPZFMatrix<STATE> & u_n)
    {
        m_u_n = u_n;
    }
    
    /** @brief get displacements at last time */
    TPZFMatrix<STATE> & u_n(){
        return m_u_n;
    }
    
    /** @brief set grad_u at last time */
    void Set_grad_u_n(TPZFMatrix<STATE> & grad_u_n)
    {
        m_grad_u_n = grad_u_n;
    }
    
    /** @brief get grad_u at last time */
    TPZFMatrix<STATE> & grad_u_n()
    {
        return m_grad_u_n;
    }
    
    
    /** @brief Set elastic strain at last time */
    void Set_epsilon_e_n(TPZFMatrix<REAL> & epsilon_e_n)
    {
        m_epsilon_e_n = epsilon_e_n;
    }
    
    /** @brief Get elastic strain at last time */
    TPZFMatrix<REAL> epsilon_e_n()
    {
        return m_epsilon_e_n;
    }
    
    
    /** @brief Set plastic strain at last time */
    void Set_epsilon_p_n(TPZFMatrix<REAL> & epsilon_p_n)
    {
        m_epsilon_p_n = epsilon_p_n;
    }
    
    /** @brief Get plastic strain at last time */
    TPZFMatrix<REAL> epsilon_p_n(){
        return m_epsilon_p_n;
    }
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /** @brief set weighted pressure at intial state */
    void Set_p_0(STATE & p_0)
    {
        m_p_0 = p_0;
    }
    
    /** @brief get weighted pressure at intial state */
    STATE & p_0(){
        return m_p_0;
    }
    
    
    /** @brief set weighted pressure at current state */
    void Set_p_n(STATE & p_n)
    {
        m_p_n = p_n;
    }
    
    /** @brief get weighted pressure at current state */
    STATE & p_n(){
        return m_p_n;
    }
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Memory :
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    /**
     Write method not implemented.
     
     @param buf TPZStream buffer
     @param withclassid obsolete
     */
    void Write(TPZStream &buf, int withclassid) const
    {
//                buf.Write(&m_Pressure_n);
//                buf.Write(&m_Pressure);
        DebugStop();
    }
    
    /**
     Read method not implemented.
     
     @param buf TPZStream buffer
     @param context obsolete
     */
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

    
    /**
     Class unique identifier

     @return integer class id
     */
    virtual int ClassId() const;
    
};

inline std::ostream &operator<<(std::ostream &out,const TPZPoroPermMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* TPZPoroPermMemory_h */
