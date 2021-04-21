/**
    \file TPZStructMatrixT.h
    Template class defining common behaviour for the TPZStructMatrix hierarchy
*/

#ifndef TPZSTRUCTMATRIXT_H
#define TPZSTRUCTMATRIXT_H

#include "TPZStructMatrix.h"

/*!
  Describes the interface that should be implemented for a given 
  structural matrix, taking into account its type.

  It is expected that the child classes will be created as:

  template <
            typename TVar=STATE,
            typename TPar=TPZStructMatrixOR<TVar>
           >
  class TPZStrDerived :
                        public TPZStructMatrixT<TVar>,
                        public TPar
  {
    ...
  };

  Meaning that it will have one template parameter corresponding to its
  type and one template parameter corresponding to its Parallel Strategy,
  a class derived from \ref TPZStrMatParInterface.

*/
template<class TVar>
class TPZStructMatrixT : public TPZStructMatrix{
public:
    //Getting constructors from parent class
    using TPZStructMatrix::TPZStructMatrix;
    /**
     *  Methods to be implemented in child classes
     */
    //@{
    //! Clone method
    TPZStructMatrix *Clone() override = 0;
    //!Creates a matrix for assembling
    TPZMatrix<TVar> * Create() override = 0;
    

    int ClassId() const override{
        return Hash("TPZStructMatrixT") ^
            TPZStructMatrix::ClassId() << 1 ^
            ClassIdOrHash<TVar>() << 2;
    }

    /* 
     * The following are optional methods
     * that can be implemented if needed:
     * //! Creates solver matrix and assembles it alongside rhs
     * TPZMatrix<TVar> *CreateAssemble(TPZBaseMatrix &,
     *                 TPZAutoPointer<TPZGuiInterface>);
     *
     * //! Operations to be performed at the beginning of CreateAssemble
     * void InitCreateAssemble() override;
     *
     */
    
    //@}
};

extern template class TPZStructMatrixT<STATE>;
#endif
