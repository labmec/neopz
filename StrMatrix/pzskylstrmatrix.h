/**
 * @file
 * @brief Contains the TPZSkylineStructMatrix class which implements SkyLine Structural Matrices.
 */

#ifndef TPZSKYLINESTRUCTMATRIX_H
#define TPZSKYLINESTRUCTMATRIX_H

#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"
class TPZCompMesh;

/**
 * @brief Implements SkyLine Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSkylineStructMatrix : public TPZStructMatrix,
                               public TPZStructMatrixOR<STATE> {
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
	
    TPZMatrix<STATE> * Create() override;
	
    TPZStructMatrix * Clone() override;
    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
    friend TPZPersistenceManager;
};

#endif //TPZSKYLINESTRUCTMATRIX_H
