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
    // Store all the data required for the integration points
    // Store the saturation at n step
    // Also it can store the nonlinear part of the flux at n step
    // Store the xyz of the spatial properties
    /// Flux from left to right at the integration point
    STATE fNormalFlux;
    /// Saturation at the left
    STATE fLeftSaturation;
    
    /// Saturation at the right
    STATE fRightSaturation;

public:

    TRMPhaseInterfaceMemory()
    {
        DebugStop();
    }

    TRMPhaseInterfaceMemory(const TRMPhaseInterfaceMemory &copy)
    {
        DebugStop();
    }
    
    TRMPhaseInterfaceMemory &operator=(const TRMPhaseInterfaceMemory &copy)
    {
        DebugStop();
        return *this;
    }
    
// Note describe this class into the lyx doc
void UpdateSolutionMemory()
{
    //update saturation and pressure and total flux (un = unp1)
//    fPressureN = fPressureNp1;
}

void Write(TPZStream &buf, int withclassid)
{
//    buf.Write(&fPressureN);
//    buf.Write(&fPressureNp1);
}

void Read(TPZStream &buf, void *context)
{
//    buf.Read(&fPressureN);
//    buf.Read(&fPressureNp1);
}

void Print(std::ostream &out) const
{
//    out << fPressureN;
//    out << fPressureNp1;
}


};

inline std::ostream &operator<<(std::ostream &out,const TRMPhaseInterfaceMemory &mem)
{
    mem.Print(out);
    return out;
}



#endif /* defined(__PZ__TRMPhaseInterfaceMemory__) */
