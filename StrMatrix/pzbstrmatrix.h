/**
 * @file
 * @brief Contains the TPZBandStructMatrix class which implements Banded Structural Matrices.
 */

#ifndef TPZBANDSTRUCTMATRIX_H
#define TPZBANDSTRUCTMATRIX_H
#include "TPZStructMatrix.h"
#include "pzstrmatrixor.h"
/**
 * @brief Implements Banded Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar = STATE,
         class TPar = TPZStructMatrixOR<TVar>>
class TPZBandStructMatrix : public TPZStructMatrix,
                            public TPar{
public:    
	
    TPZBandStructMatrix(TPZCompMesh *);
    TPZBandStructMatrix(TPZAutoPointer<TPZCompMesh>);
	
    TPZBandStructMatrix(const TPZBandStructMatrix &copy) : TPZStructMatrix(copy)
    {
    }
	
    TPZMatrix<TVar> * Create() override;
	TPZStructMatrix * Clone() override;
    
    using TPZStructMatrix::CreateAssemble;
    virtual TPZMatrix<TVar> * CreateAssemble(TPZFMatrix<TVar> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);


	//@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
private :
    TPZBandStructMatrix() = default;
    
    friend TPZPersistenceManager;
};

#endif //TPZBANDSTRUCTMATRIX_H
