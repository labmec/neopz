/**
 * @file
 * @brief Contains the implementation of the TPZMatRedStructMatrix methods. 
 */
/*
 *  tpzmatredstructmatrix.cpp
 *  SubStruct
 *
 *  Created by Philippe Devloo on 22/04/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */

#include "tpzmatredstructmatrix.h"
#include "pzskylstrmatrix.h"
#include "tpzverysparsematrix.h"
#include "pzsubcmesh.h"
#include "pzmatred.h"

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::TPZMatRedStructMatrix(TPZSubCompMesh *mesh) : TPZStructMatrix(mesh)
{
	fInternalEqs = mesh->NumInternalEquations();
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::~TPZMatRedStructMatrix()
{
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::TPZMatRedStructMatrix(const TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix> &copy) : TPZStructMatrix(copy)
{
	fInternalEqs = copy.fInternalEqs;
}

template< class TStructMatrix, class TSparseMatrix>
TPZStructMatrix *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Clone()
{
	return new TPZMatRedStructMatrix(*this);
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatrix<REAL> *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Create()
{
	int neq = fMesh->NEquations();
	TStructMatrix strmat(fMesh);
	strmat.SetEquationRange(0,fInternalEqs);
	TPZMatrix<REAL> *InternalStiff = strmat.Create();
	TPZMatRed<REAL, TSparseMatrix> *matred = new TPZMatRed<REAL, TSparseMatrix>(neq,fInternalEqs);
	matred->SetK00(InternalStiff);
	return matred;
}

template class TPZMatRedStructMatrix<TPZSkylineStructMatrix,TPZVerySparseMatrix<REAL> >;
template class TPZMatRedStructMatrix<TPZSkylineStructMatrix,TPZFMatrix<REAL> >;
