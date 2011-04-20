
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tlink.hh
//
// Class:  TPZLink<class>
//
// Obs.:   Implementa uma lista ligada com elementos indefinidos.
//
// Versao: 12 / 1994.
//


#ifndef _TLINKHH_
#define _TLINKHH_

#include <stdio.h>


#ifdef WORKPOOL

class TPZWorkPool;

#endif


/**
 * @brief Implements a linked list of ElemType elements
 * @ingroup matrixutility
 */
template< class ElemType >
class TPZLink
{
  /**
     \struct Node
     Defines a Node structure that contains an element type index and a node index.
     @param elem Describes element type.
     @param next Describes next element.
   */
  struct Node
  {
    ElemType elem;
    Node     *next; 
  };

 public:
  /**
   * Simple constructor
   */
  TPZLink();
  /**
   * Simple destructor
   */
  ~TPZLink();

#ifdef WORKPOOL
  /**
   * Sets a memory block for auxiliar manipulations
   * @param *wp Available workpool
   */
  void SetWorkPool(TPZWorkPool *wp);
#endif

  /**
   * @name Manipulators
   * Those methods implement manipulations routine with the linked list
   */
  //@{
  /**
   * Inserts a element on the list
   * @param &elem Element being inserted
   */
  int Insert( ElemType &elem );
  /**
   * Appends an element to the list
   * @param &elem Element being appended
   */
  int Append( ElemType &elem );
  /**
   * Removes an element from the list
   */
  int Remove();
  /**
   * Updates the current list
   * @param &elem Updated element on the list
   */
  int Update( ElemType &elem );
  /**
   * Clears the entire list
   */
  int Clear();
  //@}

  TPZLink<ElemType> &operator=( TPZLink<ElemType> & );
  /**
   * Returns to the head of the list
   */
  int Head();
  /**
   * Moves to the next element on list
   */
  int Next();
  /**
   * Returns an element from the list
   * @param *pElem contains the returned element
   */
  int Get( ElemType *pElem );
  /**
   * Returns the node's element type
   */
  ElemType *GetNode();
  /**
   * Returns the last element on the list
   * @param *pELem contains the last element
   */
  int GetLast( ElemType *pElem );

 private:
  /**
   * Pointer to head of the list
   */
  Node *fHead;
  /**
   * Pointer to last element on list
   */
  Node *fLast;
  /**
   * Pointer to current element
   */
  Node *fThis;
  /**
   * What the hell
   */
  Node **fpBefore;
#ifdef WORKPOOL
  /**
   * Pointer to workpool
   */
  TPZWorkPool *fWp;
#endif

};



/*** Head ***/

template< class ElemType >
inline int
TPZLink<ElemType>::Head()
{
  fpBefore = &fHead;
  return( (fThis = fHead) != NULL );
}



/*** Next ***/

template< class ElemType >
inline int
TPZLink<ElemType>::Next()
{
  if ( fThis == NULL )
    return( 0 );

  fpBefore = &fThis->next;
  fThis    =  fThis->next;
  return( 1 );
}



/*** Get ***/

template< class ElemType >
inline int
TPZLink<ElemType>::Get( ElemType *pElem )
{
  if ( fThis == NULL )
    return( 0 );

  *pElem = fThis->elem;
  return( 1 );
}

/*** GetNode ***/

template< class ElemType >
inline ElemType *
TPZLink<ElemType>::GetNode()//TPZLink<ElemType>::GetNode()
{
  if(fThis) return &(fThis->elem);
  else return NULL;
}


/*** Get Last ***/

template< class ElemType >
inline int
TPZLink<ElemType>::GetLast( ElemType *pElem )
{
  if ( fLast == NULL )
    return( 0 );

  *pElem = fLast->elem;
  return( 1 );
}

//#include "tlink.h2"

#ifdef WORKPOOL
template< class ElemType >
inline void
TPZLink<ElemType>::SetWorkPool( TPZWorkPool *wp )
{
  fWp = wp;
}
#endif

#endif
