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


class TRMMemory {

// Store all the data required for the integration points
// Store the saturation at n step
// Also it can store the nonlinear part of the flux at n step
// Store the xyz of the spatial properties
    
// Note describe this class into the lyx doc
    void UpdateSolutionMemory(); //update saturation and pressure and total flux (un = unp1)
    
};





#endif /* defined(__PZ__TRMMemory__) */


// Optional:

// Gavitational Segregation flux values are storaged here!
