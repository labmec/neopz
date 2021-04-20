/**
 * @file
 * @brief Contains the TPZStructMatrixTBBFlow class which responsible for a interface among Matrix and Finite Element classes using TBB library. Usage of this class without USING_TBB=ON when configuring the NeoPZ library will result in runtime errors.
 */

#ifndef TPZStructMatrixTBBFlow_H
#define TPZStructMatrixTBBFlow_H


#include "TPZStrMatParInterface.h"
#include "TPZSemaphore.h"
#include <mutex>

//forward declarations
class TPZElementMatrix;
class TPZBaseMatrix;
class TPZStructMatrix;
class TPZFlowGraph;

/**
 * @brief It is responsible for a interface among Matrix and Finite Element classes. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZStructMatrixTBBFlow : public virtual TPZStrMatParInterface
{
    
public:
    //! Default constructor.
    TPZStructMatrixTBBFlow();
    //! Copy constructor
    TPZStructMatrixTBBFlow(const TPZStructMatrixTBBFlow &copy) = default;
    //! Move constructor
    TPZStructMatrixTBBFlow(TPZStructMatrixTBBFlow &&copy) = default;
    //! Destructor
    ~TPZStructMatrixTBBFlow();
    //! Copy assignment operator
    TPZStructMatrixTBBFlow& operator=(const TPZStructMatrixTBBFlow &) = default;
    //! Move assignment operator
    TPZStructMatrixTBBFlow& operator=(TPZStructMatrixTBBFlow &&) = default;

        /**
     * Functions overriden from TPZStrMatParInterface
     */
    //@{

    void Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    void Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    
    TPZBaseMatrix * CreateAssemble(TPZBaseMatrix &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface) override;
    //@}
    
    //@{
    //!Read and Write methods
    int ClassId() const override;
    
    void Read(TPZStream &buf, void *context) override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    //@}

protected:

    
    /** @brief Assemble the global right hand side */
    virtual void MultiThread_Assemble(TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    /** @brief Assemble the global system of equations into the matrix which has already been created */
    virtual void MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    TPZAutoPointer<TPZFlowGraph> fFlowGraph;
};

#endif
