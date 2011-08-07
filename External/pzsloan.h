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
//  void Resequence(int * jj, int * jk, int n_nodes, int n_elements, int * nnn, int old_profile, int new_profile);
//  void Resequence(int n_nodes, int n_elements, int *nnn,int *npn, int *xnpn, int old_profile, int new_profile);
  virtual void Resequence(TPZVec<int> &perm, TPZVec<int> &iperm);
  
 private:

  int fMaxNodesElement;
/*
  int fNNodes;

  int fNElements;
	*/
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

