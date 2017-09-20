/**
 * @file
 * @brief Contains the TPZStructMatrixST class which responsible for a interface among Matrix and Finite Element classes.
 */

#ifndef TPZStructMatrixST_H
#define TPZStructMatrixST_H

#include <set>
#include <map>
#include <semaphore.h>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "pzequationfilter.h"
#include "TPZGuiInterface.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

class TPZStructMatrixST;
#include "TPZStructMatrixBase.h"


/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixST : public TPZStructMatrixBase {
    
public:
    
    TPZStructMatrixST(TPZCompMesh *);
    
    TPZStructMatrixST(TPZAutoPointer<TPZCompMesh> cmesh);
    
    TPZStructMatrixST(const TPZStructMatrixST &copy);
    
    virtual ~TPZStructMatrixST(){};
        
    /** @brief */
    virtual TPZMatrix<STATE> * Create();
    
    /** @brief */
    virtual TPZStructMatrixST * Clone();

    /** @brief */
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                                              unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual void Assemble(TPZMatrix<STATE> &mat, TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
                          unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }
    
    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
protected:
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void OnlyAssemble(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global right hand side */
    virtual void OnlyAssemble(TPZFMatrix<STATE> *rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief */
    virtual void ExecuteAssemble(TPZMatrix<STATE> *fGlobMatrix, TPZFMatrix<STATE> *fGlobRhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
};

#endif
