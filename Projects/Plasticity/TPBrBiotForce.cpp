//
//  TPBrBiotForce.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/13/14.
//
//

#include "TPBrBiotForce.h"
#include "pzmaterialid.h"


/** @brief Define the class id associated with the class */
/**
 * This id has to be unique for all classes
 * A non unique id is flagged at the startup of the program
 */
int TPBrBiotForce::ClassId() const
{
    return TPZBiotForceID;
}

/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
void TPBrBiotForce::Write(TPZStream &buf, int withclassid) 
{
    TPZFunction<STATE>::Write(buf,withclassid);
    buf.Write(&fRwell);
    buf.Write(&fRreservoir);
    buf.Write(&fPwell);
    buf.Write(&fPreservoir);
    buf.Write(&fBiot);
}

/** @brief read objects from the stream */
void TPBrBiotForce::Read(TPZStream &buf, void *context)
{
    TPZFunction<STATE>::Read(buf,context);
    buf.Read(&fRwell);
    buf.Read(&fRreservoir);
    buf.Read(&fPwell);
    buf.Read(&fPreservoir);
    buf.Read(&fBiot);
    ComputeConstant();
}

template class TPZRestoreClass<TPBrBiotForce,TPZBiotForceID>;
