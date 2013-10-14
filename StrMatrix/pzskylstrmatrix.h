/**
 * @file
 * @brief Contains the TPZSkylineStructMatrix class which implements SkyLine Structural Matrices.
 */

#ifndef TPZSKYLINESTRUCTMATRIX_H
#define TPZSKYLINESTRUCTMATRIX_H

#include "pzstrmatrix.h"

class TPZCompMesh;
template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

/**
 * @brief Implements SkyLine Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSkylineStructMatrix : public TPZStructMatrix {
protected:
    
    /** @brief the equations which should actually be assembled */
    TPZVec<long> fActiveEquations;
    
    /** @brief Equation destination */
    TPZVec<long> fEquationDestination;
    
    /** Returns the skyline matrix object */
    virtual TPZMatrix<STATE> * ReallyCreate(long neq, const TPZVec<long> &skyline);
    
public:    
	
	TPZSkylineStructMatrix(TPZCompMesh *);
    
    TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp);
    
    ~TPZSkylineStructMatrix();
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZStructMatrix * Clone();
    
public:
	
};

#endif //TPZSKYLINESTRUCTMATRIX_H
