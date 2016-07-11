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
    

    void Write(TPZStream &buf, int withclassid)
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
