/**
 * @file
 * @brief Contains the TPBSpStructMatrix class which assembles on the pair equations.
 */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"

/** 
 * @ingroup structural
 * @brief Assembles only the pair equations. \ref structural "Structural Matrix"
 */
class TPBSpStructMatrix : public TPZSpStructMatrix {
public:    
    public:
int ClassId() const override;


    virtual TPZBaseMatrix * Create() override;
	
    virtual TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
	
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
