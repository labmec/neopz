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
 * @brief Implements a sparse structural matrix using TPZFYsmpMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSpStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar{
public:
    using TPZStructMatrixT<TVar>::TPZStructMatrixT;
    
    TPZMatrix<TVar> *Create() override;
	TPZStructMatrix * Clone() override;

    void EndCreateAssemble(TPZBaseMatrix *) override;

    //@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
protected:
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
    friend TPZPersistenceManager;
};

#endif //TPZSPSTRUCTMATRIX_H
