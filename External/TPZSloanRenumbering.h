//---------------------------------------------------------------------------

#ifndef TPZSloanRenumberingH
#define TPZSloanRenumberingH

#include "pzrenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <list>
#include "pzstack.h"

class TPZSloanRenumbering : public TPZRenumbering {

  private:

    long W1() const{ return 1; }
    long W2() const{ return 2; }

    enum ENodeStatus
    {
      EInactive = 0,
      EPreActive = 1,
      EActive = 2,
      EPostActive = 3
    };

    long FindHighestPriority(const std::list<long> &Q,const TPZVec<long> &priority) const;



  public:


    virtual void Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter);

    TPZSloanRenumbering(long NElements, long NNodes);

    virtual ~TPZSloanRenumbering();

};
#endif
