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
    TPZVec<int> fActiveEquations;
    
    /** @brief Equation destination */
    TPZVec<int> fEquationDestination;
    
    /** Returns the skyline matrix object */
    virtual TPZMatrix<STATE> * ReallyCreate(int neq, const TPZVec<int> &skyline);
    
public:    
	
	TPZSkylineStructMatrix(TPZCompMesh *);
    
    TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp);
    
    ~TPZSkylineStructMatrix();
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
    
public:
	
};

#endif //TPZSKYLINESTRUCTMATRIX_H
