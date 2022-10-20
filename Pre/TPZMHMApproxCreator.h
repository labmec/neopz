


#ifndef PZ_TPZMHMCREATOR_H
#define PZ_TPZMHMCREATOR_H

#include "pzvec.h"

class TPZMHMApproxCreator
{
private:
    /* data */
public:

    TPZVec<int64_t> fElementPartition;

    TPZMHMApproxCreator(/* args */);
    ~TPZMHMApproxCreator();
};


#endif