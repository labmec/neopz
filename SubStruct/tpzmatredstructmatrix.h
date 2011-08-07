/**
 * @file
 * @brief Contains the TPZMatRedStructMatrix class. 
 */
/*
 *  tpzmatredstructmatrix.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 22/04/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

#ifndef TPZMATREDSTRUCTMATRIX
#define TPZMATREDSTRUCTMATRIX

#include "pzstrmatrix.h"
class TPZSubCompMesh;
class TPZFMatrix;

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub Structure"
 */
template< class TStructMatrix, class TSparseMatrix>
class TPZMatRedStructMatrix : TPZStructMatrix
{
public:
	
	TPZMatRedStructMatrix(TPZSubCompMesh *mesh);
	
	virtual ~TPZMatRedStructMatrix();
	
	TPZMatRedStructMatrix(const TPZMatRedStructMatrix &copy);
	
	virtual TPZStructMatrix *Clone();
	
	virtual TPZMatrix *Create();
	
private:
	
	int fInternalEqs;
	
};

#endif
