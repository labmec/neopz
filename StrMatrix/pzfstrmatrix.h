/**
 * @file
 * @brief Contains the TPZFStructMatrix class which implements Full Structural Matrices.
 */

#ifndef TPZFSTRUCTMATRIX_H
#define TPZFSTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"

class TPZCompMesh;

/**
 * @brief Implements Full Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZFStructMatrix : public TPZStructMatrix,
                         public TPZStructMatrixOR
{
public:    
	
    TPZFStructMatrix(TPZCompMesh *);
	
    TPZFStructMatrix(TPZAutoPointer<TPZCompMesh> );
    
    TPZMatrix<STATE>* Create() override;
	
    TPZStructMatrix * Clone() override;
    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
	
    TPZFStructMatrix() = default;
    
    friend TPZPersistenceManager;
};

#endif //TPZFSTRUCTMATRIX_H
