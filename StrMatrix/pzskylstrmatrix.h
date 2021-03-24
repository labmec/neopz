/**
 * @file
 * @brief Contains the TPZSkylineStructMatrix class which implements SkyLine Structural Matrices.
 */

#ifndef TPZSKYLINESTRUCTMATRIX_H
#define TPZSKYLINESTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzstrmatrix.h"

class TPZCompMesh;

/**
 * @brief Implements SkyLine Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSkylineStructMatrix : public TPZStructMatrix {
protected:
    
    /** @brief the equations which should actually be assembled */
    TPZVec<int64_t> fActiveEquations;
    
    /** @brief Equation destination */
    TPZVec<int64_t> fEquationDestination;
    
    /** Returns the skyline matrix object */
    virtual TPZMatrix<STATE> * ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline);
    
public:    
	
	TPZSkylineStructMatrix(TPZCompMesh *);
    
    TPZSkylineStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);
	
	TPZSkylineStructMatrix(const TPZSkylineStructMatrix &cp);
    
    ~TPZSkylineStructMatrix();
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZStructMatrix * Clone();

private :
    TPZSkylineStructMatrix();

    friend TPZPersistenceManager;
};

#endif //TPZSKYLINESTRUCTMATRIX_H
