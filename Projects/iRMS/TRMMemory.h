//
//  TRMMemory.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
//

#ifndef __PZ__TRMMemory__
#define __PZ__TRMMemory__

#include <stdio.h>
#include "pzreal.h"
#include "pzfilebuffer.h"


class TRMMemory {

// Store all the data required for the integration points
// Store the saturation at n step
// Also it can store the nonlinear part of the flux at n step
// Store the xyz of the spatial properties
    /// Pressure at the previous timestep
    STATE fPressureN;
    /// Pressure at the last iteration
    STATE fPressureNp1;
    
public:
    
// Note describe this class into the lyx doc
    void UpdateSolutionMemory()
    {
        //update saturation and pressure and total flux (un = unp1)
        fPressureN = fPressureNp1;
    }
  
    void Write(TPZStream &buf, int withclassid)
    {
        buf.Write(&fPressureN);
        buf.Write(&fPressureNp1);
    }

    void Read(TPZStream &buf, void *context)
    {
        buf.Read(&fPressureN);
        buf.Read(&fPressureNp1);
    }

    void Print(std::ostream &out) const
    {
        out << fPressureN;
        out << fPressureNp1;
    }
    
    
};

inline std::ostream &operator<<(std::ostream &out,const TRMMemory &mem)
{
    mem.Print(out);
    return out;
}



#endif /* defined(__PZ__TRMMemory__) */


// Optional:

// Gavitational Segregation flux values are storaged here!
