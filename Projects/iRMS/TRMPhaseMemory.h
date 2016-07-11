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
    
    /** @brief contains the volumetric average pressure */
    REAL fp;

    /** @brief contains the volumetric average pressure at last time step */
    REAL fp_n;
    
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
        fp      = copy.fp;
        fp_n    = copy.fp_n;
        fsa     = copy.fsa;
        fsa_n   = copy.fsa_n;
        fsb     = copy.fsb;
        fsb_n   = copy.fsb_n;
    }
    
    TRMPhaseMemory &operator=(const TRMPhaseMemory &other)
    {
        
        if (this != & other) // prevent self-assignment
        {
            fp      = other.fp;
            fp_n    = other.fp_n;
            fsa     = other.fsa;
            fsa_n   = other.fsa_n;
            fsb     = other.fsb;
            fsb_n   = other.fsb_n;
        }
        return *this;
        
    }

void Write(TPZStream &buf, int withclassid)
{
    buf.Write(&fp);
    buf.Write(&fp_n);
    buf.Write(&fsa);
    buf.Write(&fsa_n);
    buf.Write(&fsb);
    buf.Write(&fsb_n);
}

void Read(TPZStream &buf, void *context)
{
    buf.Read(&fp);
    buf.Read(&fp_n);
    buf.Read(&fsa);
    buf.Read(&fsa_n);
    buf.Read(&fsb);
    buf.Read(&fsb_n);
}

void Print(std::ostream &out) const
{
    
    out << "TRMPhaseMemory item, with values ";
    out << fp;
    out << fp_n;
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
