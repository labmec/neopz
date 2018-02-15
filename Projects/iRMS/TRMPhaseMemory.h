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


class TRMPhaseMemory
{

private:
    
    /** @brief contains the normal flux per surface area */
    REAL fun;
    
    /** @brief contains the volumetric average pressure */
    REAL fp_avg;

    /** @brief contains the volumetric average pressure at last time step */
    REAL fp_avg_n;
    
    /** @brief contains the volumetric saturation of alpha phase */
    REAL fsa;
    
    /** @brief contains the volumetric saturation of alhpa at last time step */
    REAL fsa_n;
    
    /** @brief contains the volumetric saturation of beta phase */
    REAL fsb;

    /** @brief contains the volumetric saturation of alhpa at last time step */
    REAL fsb_n;

    
public:

    
    /** @brief Default constructor */
    TRMPhaseMemory();
    
    /** @brief Default destructor */
    ~TRMPhaseMemory();
    
    TRMPhaseMemory(const TRMPhaseMemory &copy)
    {
        fp_avg      = copy.fp_avg;
        fp_avg_n    = copy.fp_avg_n;
        fsa         = copy.fsa;
        fsa_n       = copy.fsa_n;
        fsb         = copy.fsb;
        fsb_n       = copy.fsb_n;
    }
    
    TRMPhaseMemory &operator=(const TRMPhaseMemory &other)
    {
        
        if (this != & other) // prevent self-assignment
        {
            fp_avg      = other.fp_avg;
            fp_avg_n    = other.fp_avg_n;
            fsa         = other.fsa;
            fsa_n       = other.fsa_n;
            fsb         = other.fsb;
            fsb_n       = other.fsb_n;
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
    
    

    void Write(TPZStream &buf, int withclassid) const{
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
