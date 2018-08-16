#ifndef TPZStructMatrixGCTP_H
#define TPZStructMatrixGCTP_H

#include <set>
#include <map>
#include <semaphore.h>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "TPZEquationFilter.h"
#include "TPZGuiInterface.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

class TPZStructMatrixGCTP;
#include "TPZStructMatrixBase.h"

#ifdef USING_BOOST
#include <boost/atomic.hpp>
#endif

/**
 * @brief It is responsible for a interface between Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * This class uses graph coloring and TPZThreadPool to assemble the matrix in parallel.
 * @ingroup structural
 */
class TPZStructMatrixGCTP : public TPZStructMatrixBase {
public:

    TPZStructMatrixGCTP() : TPZStructMatrixBase() {
    }

    TPZStructMatrixGCTP(TPZCompMesh *);

    TPZStructMatrixGCTP(TPZAutoPointer<TPZCompMesh> cmesh);

    TPZStructMatrixGCTP(const TPZStructMatrixGCTP &copy);

    virtual ~TPZStructMatrixGCTP() {
    };

    virtual TPZMatrix<STATE> * Create();

    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
            unsigned numthreads_assemble, unsigned numthreads_decompose) {
        SetNumThreads(numthreads_assemble);
        return CreateAssemble(rhs, guiInterface);
    }

    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    virtual TPZStructMatrixGCTP * Clone();

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    virtual void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface,
            unsigned numthreads_assemble, unsigned numthreads_decompose) {
        std::cout << "Nothing to do." << std::endl;
    }

    /** @brief Assemble the global right hand side */
    virtual void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

public:
    virtual int ClassId() const;
    void Read(TPZStream& buf, void* context);
    void Write(TPZStream& buf, int withclassid) const;

protected:

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

public:

    /** @brief Find the order to assemble the elements */
    static void OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder);

    /** @brief Establish whether the element should be computed */
    bool ShouldCompute(int matid) const {
        const size_t size = fMaterialIds.size();
        return size == 0 || fMaterialIds.find(matid) != fMaterialIds.end();
    }

    /** @brief Returns the material ids */
    const std::set<int> &MaterialIds() {
        return fMaterialIds;
    }

protected:

    TPZManVector<int64_t> fElementOrder;
    TPZVec<int64_t> fElementColors;
    int64_t fNColors;
};

#endif
