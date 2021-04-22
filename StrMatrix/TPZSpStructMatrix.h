/**
 * @file
 * @brief Contains the TPZSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSPSTRUCTMATRIX_H
#define TPZSPSTRUCTMATRIX_H

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"
#include "pzstack.h"
/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSpStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar{
public:
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    
    TPZMatrix<TVar> *Create() override;
	TPZStructMatrix * Clone() override;
    
	TPZMatrix<TVar> *CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;

    //@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
protected:
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
