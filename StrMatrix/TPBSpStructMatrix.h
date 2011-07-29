/* Generated by Together */

#ifndef TPBSPSTRUCTMATRIX_H
#define TPBSPSTRUCTMATRIX_H

#include "TPZSpStructMatrix.h"

struct TPZElementMatrix;
class TPZCompMesh;
class TPZFMatrix;
class TPZMatrix;
class TPZStructMatrix;

void UniformRefine(int num, TPZGeoMesh &m);

/** @ingroup structural
 * @brief This matrix assembles on the pair equations
 */
class TPBSpStructMatrix : public TPZSpStructMatrix {
public:    
	
    virtual TPZMatrix * Create();    
	
    virtual TPZMatrix * CreateAssemble(TPZFMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
    /** Used only for testing */
    static int main();
	
    TPBSpStructMatrix(TPZCompMesh *);
    
    TPBSpStructMatrix(const TPBSpStructMatrix &copy) : TPZSpStructMatrix(copy)
    {
    }
	
};

#endif //TPBSPSTRUCTMATRIX_H
