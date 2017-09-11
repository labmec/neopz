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
    // Store all the data required for the integration points
    // Store the saturation at n step
    // Also it can store the nonlinear part of the flux at n step
    // Store the xyz of the spatial properties
    /// Pressure at the integration point
    STATE fPressure;
    /// Saturation at the previous time step
    STATE fSaturationN;    
    /// Saturation at the current time step
    STATE fSaturationNP1;

public:

// Note describe this class into the lyx doc
void UpdateSolutionMemory()
{
    //update saturation and pressure and total flux (un = unp1)
    fSaturationN = fSaturationNP1;
}

void Write(TPZStream &buf, int withclassid) const{
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

inline std::ostream &operator<<(std::ostream &out,const TRMPhaseMemory &mem)
{
    mem.Print(out);
    return out;
}



#endif /* defined(__PZ__TRMPhaseInterfaceMemory__) */
