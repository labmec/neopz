/**
 * @file
 * @brief Contains TPZSpBlockDiagPivot class which does derivation using decompose LU with pivot.
 */
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

#include "tpzsparseblockdiagonal.h"

/**
 @brief Derivation using decompose LU with pivot. \ref matrix "Matrix"
 @ingroup matrix
 @author Philippe R. B. Devloo
 */
template<class TVar>
class TPZSpBlockDiagPivot : public TPZSparseBlockDiagonal<TVar>
{
public:
    TPZSpBlockDiagPivot();
	
    ~TPZSpBlockDiagPivot();
    
    virtual int Decompose_LU();
    
	virtual int Decompose_LU(std::list<int> &singular)
	{
        return TPZBlockDiagonal<TVar>::Decompose_LU(singular);
    }
	
	virtual int Substitution( TPZFMatrix<TVar> * B ) const;
	
private:
	/** @brief Attribute to store equation changes in LU decomposition.
	 */
	TPZVec<int> fPivotIndices;
	
	int Substitution2( TPZFMatrix<TVar> * B ) const;
	
};

#endif
