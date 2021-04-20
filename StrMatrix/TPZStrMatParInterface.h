/**
   \file TPZStrMatParInterface.h
   Describes the interface that should be implemented in a parallel
   scheme for \ref TPZStructMatrix.
*/


#ifndef TPZSTRUCTMATRIXPARINTERFACE_H
#define TPZSTRUCTMATRIXPARINTERFACE_H

#include "TPZSavable.h"
#include "TPZEquationFilter.h"

class TPZBaseMatrix;
class TPZGuiInterface;
template<class T>
class TPZAutoPointer;
template<class T>
class TPZVec;

/*!
  Describes the interface that should be implemented in a parallel
  scheme for \ref TPZStructMatrix
*/
class TPZStrMatParInterface : public virtual TPZSavable {
public:
    virtual ~TPZStrMatParInterface() = default;

    /**
     *  Methods to be derived in parallel schemes implementations.
     */
    
    //@{
    //! Assemble the global system of equations into a matrix that has already been created.
    virtual void Assemble(TPZBaseMatrix &stiffness, TPZBaseMatrix &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;
    //! Assemble the global right hand side vector.
    virtual void Assemble(TPZBaseMatrix &rhs,
                          TPZAutoPointer<TPZGuiInterface> guiInterface) = 0;
    //! NEED DOCU
    virtual TPZBaseMatrix *
    CreateAssemble(TPZBaseMatrix &rhs,
                   TPZAutoPointer<TPZGuiInterface> guiInterface);
    //@}

    //! Set number of threads to be used in the assembly.
    inline void SetNumThreads(int n) {
        this->fNumThreads = n;
    }
    //! Get number of threads to be used in the assembly.
    inline int GetNumThreads() const {
        return this->fNumThreads;
    }

    //@{
    //!Read and Write methods
    int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
protected:
    TPZStrMatParInterface() = default;
    TPZStrMatParInterface(const TPZStrMatParInterface&) = default;
    TPZStrMatParInterface(TPZStrMatParInterface&&) = default;
    TPZStrMatParInterface& operator=(const TPZStrMatParInterface&) = default;
    TPZStrMatParInterface& operator=(TPZStrMatParInterface&&) = default;
    //! Number of threads in Assemble process.
    int fNumThreads{0};
};

#endif