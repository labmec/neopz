/**
 * @file
 * @brief Contains the TPZSBandStructMatrix class which implements Symmetric Banded Structural Matrices.
 */

#ifndef TPZSBANDSTRUCTMATRIX_H
#define TPZSBANDSTRUCTMATRIX_H

#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"
/**
 * @brief Implements Symmetric Banded Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSBandStructMatrix : public TPZStructMatrix,
                             public TPZStructMatrixOR<STATE>
{
	
    TPZSBandStructMatrix(TPZCompMesh *);
	
	TPZSBandStructMatrix(TPZAutoPointer<TPZCompMesh>);
    
    TPZMatrix<STATE> * CreateAssemble(TPZBaseMatrix &rhs,TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
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

#endif //TPZSBANDSTRUCTMATRIX_H
