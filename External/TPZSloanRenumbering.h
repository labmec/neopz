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

    int W1() const{ return 1; }
    int W2() const{ return 2; }

    enum ENodeStatus
    {
      EInactive = 0,
      EPreActive = 1,
      EActive = 2,
      EPostActive = 3
    };

    int FindHighestPriority(const std::list<int> &Q,const TPZVec<int> &priority) const;



  public:


    virtual void Resequence(TPZVec<int> &permGather, TPZVec<int> &permScatter);

    TPZSloanRenumbering(int NElements, int NNodes);

    virtual ~TPZSloanRenumbering();

};
#endif
