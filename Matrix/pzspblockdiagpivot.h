//
// C++ Interface: %{MODULE}
//
// Description: 
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZSPBLOCKDIAGPIVOT_H
#define PZSPBLOCKDIAGPIVOT_H

#include <tpzsparseblockdiagonal.h>

/**
Derivation using decompose LU with pivot.

@author Philippe R. B. Devloo
*/
class TPZSpBlockDiagPivot : public TPZSparseBlockDiagonal
{
public:
    TPZSpBlockDiagPivot();

    ~TPZSpBlockDiagPivot();
    
    virtual int Decompose_LU();
    
    virtual int Substitution( TPZFMatrix * B ) const;
     
private:
  /** Attribute to store equation changes in LU decomposition.
   */
  TPZVec<int> fPivotIndices;
  
  int Substitution2( TPZFMatrix * B ) const;

};

#endif
