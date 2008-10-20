//
// C++ Interface: tpzsparseblockdiagonalstructmatrix
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H
#define TPZSPARSEBLOCKDIAGONALSTRUCTMATRIX_H

#include "pzstrmatrix.h"

/**
This class will Build a sparse block diagonal preconditioner with a structure determined by the parameters passed to it

@author Philippe R. B. Devloo
*/
class TPZSparseBlockDiagonalStructMatrix : public TPZStructMatrix
{
public:
    TPZSparseBlockDiagonalStructMatrix(TPZCompMesh *mesh);

    ~TPZSparseBlockDiagonalStructMatrix();

  virtual TPZMatrix * Create();

    virtual TPZStructMatrix* Clone();
    int NumColors();

};

#endif
