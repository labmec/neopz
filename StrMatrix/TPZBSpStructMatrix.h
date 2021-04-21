/**
 * @file
 * @brief Contains the TPBSpStructMatrix class which assembles on the pair equations.
 */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"

#include "pzstrmatrixor.h"
/** 
 * @ingroup structural
 * @brief Assembles only the pair equations. \ref structural "Structural Matrix"
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZBSpStructMatrix : public TPZSpStructMatrix<TVar,TPar> {
public:    
    TPZBSpStructMatrix(TPZCompMesh *);
    TPZBSpStructMatrix(TPZAutoPointer<TPZCompMesh>);
    
    int ClassId() const override;

    TPZMatrix<TVar> * Create() override;
	
    TPZStructMatrix * Clone() override;

    TPZMatrix<TVar> *CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
};

#endif //TPBSPSTRUCTMATRIX_H
