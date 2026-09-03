#ifndef PZ_TPZMHMCREATOR_H
#define PZ_TPZMHMCREATOR_H

#include "pzvec.h"

class TPZMHMApproxCreator
{
protected:
    TPZVec<int64_t> fElementPartition;
    
public:

    TPZMHMApproxCreator();
    TPZMHMApproxCreator(TPZVec<int64_t>& elementPartition);
    
    ~TPZMHMApproxCreator();
};


#endif
