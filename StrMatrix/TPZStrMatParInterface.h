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
  scheme for TPZStructMatrix so anyone can implement a custom
  parallel interface. Write it such as:
  
      template<class TVar>
      class TPZStrMyParInterface: public virtual TPZStrMatParInterface{
        ...
      }

  where  `TVar` stands for `STATE` and `CSTATE`.
  @ingroup structural
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

    //! Operations to be performed at the beginning of CreateAssemble
    virtual void InitCreateAssemble(){}
    //! Operations to be performed at the end of CreateAssemble
    virtual void EndCreateAssemble(TPZBaseMatrix *stiff){}

    /*! Creates solver matrix and assembles it alongside global rhs.
     Avoid overriding it unless there are no other options*/
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
    //! Whether the rhs vector will be computed or not
    inline bool ComputeRhs() const;
    //! Set whether to compute the rhs vector
    inline void SetComputeRhs(bool rhs);
    
    //@{
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
    //! Whether the rhs will be computed
    bool fComputeRhs{true};
};


bool TPZStrMatParInterface::ComputeRhs() const
{
    return fComputeRhs;
}

void TPZStrMatParInterface::SetComputeRhs(bool rhs)
{
    fComputeRhs = rhs;
}
#endif