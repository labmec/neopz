//
// C++ Interface: tpzpermutation
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZPERMUTATION_H
#define TPZPERMUTATION_H

#include "pzmanvector.h"

/**
This class generates all permutations of n values

@author Philippe R. B. Devloo
*/
class TPZPermutation{
public:
    TPZPermutation(int n);

    ~TPZPermutation();
    
    /// Applies the current permutation on the vector in and produces the vector out
    void Permute(const TPZVec<int> &in, TPZVec<int> &out) const;
    
    void operator++();

protected:
  /// Variable which represents a counter for the permutations
    TPZManVector<int> fCounter;
    /// Variable which contains the current permutations
    TPZManVector<int> fOrder;
};

#endif
