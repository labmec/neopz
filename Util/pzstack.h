/**
 * @file pzstack.h
 * @brief A simple stack.
 */
// $Id: pzstack.h,v 1.4 2005-04-25 02:55:51 phil Exp $

#ifndef PZSTACK_H
#define PZSTACK_H

#include <iostream>
#include <cstdlib>
#include "pzerror.h"
#include "pzmanvector.h"


/**
 * @ingroup util
 *
 * This class implements a stack object.
 *
 * The TPZStack object will automatically increase its size as more
 * object are pushed. It inherits the memory management from the
 * TPZManVector class. The class T needs to implement the assignment
 * operator and copy constructor.
 */
template< class T, int NumExtAlloc = 100 >
class TPZStack : public TPZManVector< T, NumExtAlloc >
{
   public:
      /** @return A reference to the first element of the stack. */
      T & Peek() const;

      /**
       * Create the stack object, indicates the stack increments.
       *
       * The size of the stack is always equal zero upon creation.
       *
       * @param allocationsize Indicates the number of elements
       * allocated to store the elements which will be pushed onto the
       * stack.
       */
      TPZStack();

      /**
       * Pushes a copy of the object on the stack.
       *
       * @param object Element which will be copied onto the stack.
       */
      void Push(const T object);

      /**
       * Retrieve an object from the stack.
       *
       * If no objects exist on the stack, the method will return an
       * object created with the empty constructor. If the NODEBUG is
       * not defined and no * object exists on the stack, a warning
       * message is sent to PZError.
       */
      T Pop();

      /** Casting operator.
       *
       *  @return The fStore pointer.
       */
      operator T*() const { return this->fStore; }

      /**@shapeType DependencyLink*/
      /*#  TPZManVector<T> lnkUnnamed*/
};

//--| IMPLEMENTATION |----------------------------------------------------------

template<class T, int NumExtAlloc >
TPZStack<T, NumExtAlloc>::TPZStack() : TPZManVector<T, NumExtAlloc>(0) {
   this->Expand(NumExtAlloc);
}

// Puts an object on the stack
template<class T, int NumExtAlloc >
void TPZStack<T, NumExtAlloc>::Push(const T object) {
   Resize(this->NElements()+1);
   operator[](this->NElements()-1) = object;
}

// Retrieve an object from the stack
template<class T, int NumExtAlloc >
T TPZStack<T, NumExtAlloc>::Pop() {
   this->fNElements--;
   if(this->fNElements <0){
      this->fNElements = 0;
      PZError << "TPZStack popping beyond the stack object" << std::endl;
      PZError.flush();
      T temp;
      return temp;
   }
   return this->fStore[this->fNElements];
}

template <class T, int NumExtAlloc >
T & TPZStack<T, NumExtAlloc>::Peek() const {
   if(this->NElements() <= 0) {
      PZError << "TPZStack peek beyond the stack object" << std::endl;
      PZError.flush();
      exit(-1);
   }
   return operator[](this->NElements()-1);
}

#endif // PZSTACK_H

//--| PZ |----------------------------------------------------------------------
