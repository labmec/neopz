/**
 * @file
 * @brief Contains the TPZSymetricSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSymetricSpStructMatrix_H
#define TPZSymetricSpStructMatrix_H

#include "TPZStructMatrix.h"
#include "pzstack.h"

#include "pzstrmatrixor.h"
/**
 * @brief Implements Sparse Symmetric Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSymetricSpStructMatrix : public TPZStructMatrix,
                                  public TPar {
    
public:    
	
    TPZSymetricSpStructMatrix(TPZCompMesh *);
    TPZSymetricSpStructMatrix(TPZAutoPointer<TPZCompMesh>);
	
    TPZMatrix<TVar> * Create() override;
	TPZStructMatrix * Clone() override;
    
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    
	virtual TPZMatrix<TVar> * CreateAssemble(TPZFMatrix<TVar> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
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
