/**
 * @file
 * @brief Contains the TPZSloan class.
 */

#ifndef TPZSLOAN_H
#define TPZSLOAN_H

#include "TPZRenumbering.h"

#ifndef BORLAND
#include "sloan.h"
#else
#include "sloan\\sloan.h"
#endif

/**
 * @ingroup util
 * @brief Interface to sloan subrotines. \ref util "Utility"
 */
class TPZSloan : public TPZRenumbering {
 public:

	virtual void Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &iperm) override;
  
 private:

  int fMaxNodesElement;

public:
   
   TPZSloan(): TPZRenumbering(),
            fMaxNodesElement(27)
   {
   }
  TPZSloan(int NElements, int NNodes);

  virtual ~TPZSloan()
     {
     }

};

#endif

