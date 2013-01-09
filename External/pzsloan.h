/**
 * @file
 * @brief Contains the TPZSloan class.
 */

#ifndef TPZSLOAN_H
#define TPZSLOAN_H

#include "pzrenumbering.h"

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

	virtual void Resequence(TPZVec<long> &perm, TPZVec<long> &iperm);
  
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

