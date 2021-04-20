/**
 * @file
 * @brief Contains the TPZSymetricSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSymetricSpStructMatrix_H
#define TPZSymetricSpStructMatrix_H

#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"
#include "pzstack.h"
/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSymetricSpStructMatrix : public TPZStructMatrix,
                                  public TPZStructMatrixOR {
    
public:    
	
    TPZSymetricSpStructMatrix(TPZCompMesh *);
    TPZSymetricSpStructMatrix(TPZAutoPointer<TPZCompMesh>);
	
    TPZMatrix<STATE> * Create() override;
	TPZStructMatrix * Clone() override;
    
    virtual TPZMatrix<STATE> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
    
    friend TPZPersistenceManager;
};

#endif //TPZSymetricSpStructMatrix_H
