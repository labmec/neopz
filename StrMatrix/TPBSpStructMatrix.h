/**
 * @file
 * @brief Contains the TPBSpStructMatrix class which assembles on the pair equations.
 */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"
#include "pzelmat.h"
#include "pzcmesh.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

/** 
 * @ingroup structural
 * @brief Assembles only the pair equations. \ref structural "Structural Matrix"
 */
class TPBSpStructMatrix : public TPZSpStructMatrix {
public:    
    public:
int ClassId() const override;


    virtual TPZMatrix<STATE> * Create() override;
	
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
    virtual TPZStructMatrix * Clone() override;
	
    /** Used only for testing */
    static int main();
	
    TPBSpStructMatrix(TPZCompMesh *);
    
    TPBSpStructMatrix(const TPBSpStructMatrix &copy) : TPZRegisterClassId(&TPBSpStructMatrix::ClassId),
    TPZSpStructMatrix(copy)
    {
    }
	
};

#endif //TPBSPSTRUCTMATRIX_H
