//
//  TRMPhaseInterfaceMemory.h
//  PZ
//
//  Created by Philippe Devloo on 7/5/15.
//
//

#ifndef __PZ__TRMPhaseInterfaceMemory__
#define __PZ__TRMPhaseInterfaceMemory__

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"


class TRMPhaseInterfaceMemory
{

    /** @brief contains the normal flux per surface area */
    REAL fun;
    
    /** @brief contains the volumetric left average pressure at last time step */
    REAL fp_avg_n_l;
    
    /** @brief contains the volumetric left saturation of alhpa at last time step */
    REAL fsa_n_l;
    
    /** @brief contains the volumetric left saturation of alhpa at last time step */
    REAL fsb_n_l;
    
    /** @brief contains the volumetric right average pressure at last time step */
    REAL fp_avg_n_r;
    
    /** @brief contains the volumetric right saturation of alhpa at last time step */
    REAL fsa_n_r;
    
    /** @brief contains the volumetric right saturation of alhpa at last time step */
    REAL fsb_n_r;

public:

    /** @brief Default constructor */
    TRMPhaseInterfaceMemory();
    
    /** @brief Default destructor */
    ~TRMPhaseInterfaceMemory();

    /** @brief Constructor based on a copy */
    TRMPhaseInterfaceMemory(const TRMPhaseInterfaceMemory &copy){
        fun      = copy.fun;
    }

    /** @brief Assignment operator */
    TRMPhaseInterfaceMemory &operator=(const TRMPhaseInterfaceMemory &other){
        if (this != & other) // prevent self-assignment
        {
            fun      = other.fun;
        }
        return *this;
    }
    
    
    /**
     * @defgroup Set and Get methods
     * @{
     */
    
    /** @brief Set the average normal flux */
    void Set_un(REAL un){
        fun = un;
    }
    
    /** @brief Set the average normal flux */
    REAL un(){
        return fun;
    }
    
    /** @brief Set the average weighted pressure at the previous timestep */
    void Set_p_avg_n_l(REAL p_avg_n_l){
        fp_avg_n_l = p_avg_n_l;
    }
    
    /** @brief Get the average weighted pressure at the previous timestep */
    REAL p_avg_n_l(){
        return fp_avg_n_l;
    }
    
    /** @brief Set the average weighted pressure at the previous timestep */
    void Set_p_avg_n_r(REAL p_avg_n_r){
        fp_avg_n_r = p_avg_n_r;
    }
    
    /** @brief Get the average weighted pressure at the previous timestep */
    REAL p_avg_n_r(){
        return fp_avg_n_r;
    }
    
    /** @brief Set alpha saturation at last step */
    void Set_sa_n_l(REAL sa_n_l){
        fsa_n_l = sa_n_l;
    }
    
    /** @brief Get alpha saturation at last step */
    REAL sa_n_l(){
        return fsa_n_l;
    }
    
    /** @brief Set alpha saturation at last step */
    void Set_sa_n_r(REAL sa_n_r){
        fsa_n_r = sa_n_r;
    }
    
    /** @brief Get alpha saturation at last step */
    REAL sa_n_r(){
        return fsa_n_r;
    }
    
    /** @brief Set alpha saturation at last step */
    void Set_sb_n_l(REAL sb_n_l){
        fsb_n_l = sb_n_l;
    }
    
    /** @brief Get alpha saturation at last step */
    REAL sb_n_l(){
        return fsb_n_l;
    }
    
    /** @brief Set alpha saturation at last step */
    void Set_sb_n_r(REAL sb_n_r){
        fsb_n_r = sb_n_r;
    }
    
    /** @brief Get alpha saturation at last step */
    REAL sb_n_r(){
        return fsb_n_r;
    }
    
    //@}

    void Write(TPZStream &buf, int withclassid) const
    {
        buf.Write(&fun);
    }
    
    void Read(TPZStream &buf, void *context)
    {
        buf.Read(&fun);
    }
    
    void Print(std::ostream &out) const
    {
        out << "TRMPhaseInterfaceMemory item, with values ";
        out << fun;
    }


};

inline std::ostream &operator<<(std::ostream &out,const TRMPhaseInterfaceMemory &mem)
{
    mem.Print(out);
    return out;
}


#endif /* defined(__PZ__TRMPhaseInterfaceMemory__) */
