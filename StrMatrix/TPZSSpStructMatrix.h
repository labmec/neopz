/**
 * @file
 * @brief Contains the TPZSSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSSpStructMatrix_H
#define TPZSSpStructMatrix_H

#include "TPZStructMatrixT.h"
#include "pzstack.h"

#include "pzstrmatrixor.h"
/**
 * @brief Implements Sparse Symmetric Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSSpStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar {
    
public:
	using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    TPZMatrix<TVar> * Create() override;
	TPZStructMatrix * Clone() override;
    
	virtual TPZMatrix<TVar> * CreateAssemble(TPZFMatrix<TVar> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
protected:
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
private :
    
    friend TPZPersistenceManager;
};

#endif //TPZSSpStructMatrix_H
