/**
 * @file
 * @brief Contains TPZSpBlockDiagPivot class which does derivation using decompose LU with pivot.
 */

#ifndef PZSPBLOCKDIAGPIVOT_H
#define PZSPBLOCKDIAGPIVOT_H

#include "tpzsparseblockdiagonal.h"

/**
 * @brief Derivation using decompose LU with pivot. \ref matrix "Matrix"
 * @ingroup matrix
 * @author Philippe R. B. Devloo
 */
template<class TVar>
class TPZSpBlockDiagPivot : public TPZSparseBlockDiagonal<TVar>
{
public:
    TPZSpBlockDiagPivot();
	
    ~TPZSpBlockDiagPivot();
    
    virtual int Decompose_LU();
    
	virtual int Decompose_LU(std::list<int64_t> &singular)
	{
        return TPZBlockDiagonal<TVar>::Decompose_LU(singular);
    }
	
	virtual int Substitution( TPZFMatrix<TVar> * B ) const;
	
private:
	/** @brief Attribute to store equation changes in LU decomposition. */
	TPZVec<int> fPivotIndices;
	
	int Substitution2( TPZFMatrix<TVar> * B ) const;
	
};

#endif
