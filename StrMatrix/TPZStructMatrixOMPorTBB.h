//
// Created by victor on 29/08/2022.
//

/**
 * @file
 * @brief Contains the TPZStructMatrixOMPorTBB class which is a TPZStructMatrixBase using graph coloring to define the order
to process the elements and each color is processed and
synchronized.
 */

#ifndef PZ_TPZSTRUCTMATRIXOMPORTBB_H
#define PZ_TPZSTRUCTMATRIXOMPORTBB_H

#include <set>
#include <map>
//#include <semaphore>
#include <mutex>
#include <condition_variable>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzelmat.h"
#include "TPZSemaphore.h"
#include "TPZEquationFilter.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#include "TPZStrMatParInterface.h"

/**
 * @brief Refines geometrical mesh (all the elements) num times
 * @ingroup geometry
 */
//void UniformRefine(int num, TPZGeoMesh &m);

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
template<class TVar>
class TPZStructMatrixOMPorTBB : public virtual TPZStrMatParInterface{

public:

    TPZStructMatrixOMPorTBB();

    TPZStructMatrixOMPorTBB(const TPZStructMatrixOMPorTBB &copy);

    //! Move constructor
    TPZStructMatrixOMPorTBB(TPZStructMatrixOMPorTBB &&copy) = default;

    virtual ~TPZStructMatrixOMPorTBB() = default;

    //! Copy assignment operator
    TPZStructMatrixOMPorTBB& operator=(const TPZStructMatrixOMPorTBB &) = default;
    //! Move assignment operator
    TPZStructMatrixOMPorTBB& operator=(TPZStructMatrixOMPorTBB &&) = default;

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs) override;
    void Assemble(TPZBaseMatrix & rhs) override;

    /** @brief Sets buffer size to be preallocated for a TPZFMatrix provided
        to TPZElementMatrixT<TVar> for the dependency matrix, useful
    when dealing with large dependencies to avoid repeated dynamic allocation*/
    void BufferSizeForUserMatrix(const int sz){fUserMatSize=sz;}
public:
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;

protected:

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void Serial_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs);

    /** @brief Assemble the global right hand side */
    virtual void Serial_Assemble(TPZBaseMatrix & rhs);

    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZBaseMatrix & rhs);

    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs);

    void VerifyStiffnessSum(TPZBaseMatrix & mat);

    void AssemblingUsingTBBandColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs );



    void AssemblingUsingOMPandColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs );


    void AssemblingUsingTBBbutNotColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs );


    void AssemblingUsingOMPbutNotColoring(TPZBaseMatrix & mat, TPZBaseMatrix & rhs );

    void CalcStiffAndAssemble(TPZBaseMatrix & mat,TPZBaseMatrix & rhs,TPZCompEl *el,
                              TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);

public:
    void OrderElements();
    int  GetNumberColors();

    void SetTBBorOMP (bool type){
        fUsingTBB = type;
    }

    void SetShouldColor (bool type){
        fShouldColor = type;
    }

protected:
    TPZVec<int> fElVecColor;
    bool fShouldColor = false;
    bool fUsingTBB = false;

    /** @brief Vectors for mesh coloring */
    TPZVec<int64_t> fElBlocked, fElSequenceColor;

    /// vector of the size of the elements containing 0 or 1 if the element has been computed (in the order of computation sequence)
    TPZVec<int64_t> fElementsComputed;

    /// All elements below or equal this index have been computed
    int64_t fElementCompleted;

    /// variable indicating if a thread is sleeping
    int fSomeoneIsSleeping;

    std::atomic<int64_t> fCurrentIndex;

    /// user allocated mat for TPZElementMatrixT
    int fUserMatSize{0};
};

extern template class TPZStructMatrixOMPorTBB<STATE>;
extern template class TPZStructMatrixOMPorTBB<CSTATE>;

#endif //PZ_TPZSTRUCTMATRIXOMPORTBB_H
