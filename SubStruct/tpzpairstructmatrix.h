/*
 *  tpzpairstructmatrix.h
 *  SubStruct
 *
 *  Created by Philippe Devloo on 20/04/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */
#ifndef PAIRSTRUCTMATRIX
#define PAIRSTRUCTMATRIX

#include "pzcmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzmatrix.h"

class TPZPairStructMatrix
{
	TPZCompMesh *fMesh;
	TPZVec<int> fPermuteScatter;
	
	void PermuteScatter(TPZVec<int> &index);
	
public:
	
	TPZPairStructMatrix(TPZCompMesh *mesh, TPZVec<int> &permutescatter)
	{
		fMesh = mesh;
		fPermuteScatter = permutescatter;
	}
	
	
	void Assemble(int mineq, int maxeq, TPZMatrix *first, TPZMatrix *second, TPZFMatrix &rhs);
	
};

#endif
