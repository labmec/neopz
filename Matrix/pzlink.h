/**
 * @file
 * @brief Contains TPZLink class which implements a linked list of ElemType elements.
 */
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
 * @brief Implements a linked list of ElemType elements. \ref matrixutility "Matrix utility"
 * @ingroup matrixutility
 */
template< class ElemType >
class TPZLink
{
	/**
     @brief Node structure
     @brief Defines a Node structure that contains an element type index and a node index.
     @param elem Describes element type.
     @param next Describes next element.
	 */
	struct Node
	{
		ElemType elem;
		Node     *next; 
	};
	
public:
	/** @brief Simple constructor */
	TPZLink();
	/** @brief Simple destructor */
	~TPZLink();
	
#ifdef WORKPOOL
	/**
	 * @brief Sets a memory block for auxiliar manipulations
	 * @param *wp Available workpool
	 */
	void SetWorkPool(TPZWorkPool *wp);
#endif
	
	/**
	 * @name Manipulators
	 * @brief Those methods implement manipulations routine with the linked list
	 */
	//@{
	/**
	 * @brief Inserts a element on the list
	 * @param &elem Element being inserted
	 */
	int Insert( ElemType &elem );
	/**
	 * @brief Appends an element to the list
	 * @param &elem Element being appended
	 */
	int Append( ElemType &elem );
	/**
	 * @brief Removes an element from the list
	 */
	int Remove();
	/**
	 * @brief Updates the current list
	 * @param &elem Updated element on the list
	 */
	int Update( ElemType &elem );
	/**
	 * @brief Clears the entire list
	 */
	int Clear();
	//@}
	
	TPZLink<ElemType> &operator=( TPZLink<ElemType> & );
	/**
	 * @brief Returns to the head of the list
	 */
	int Head();
	/**
	 * @brief Moves to the next element on list
	 */
	int Next();
	/**
	 * @brief Returns an element from the list
	 * @param *pElem contains the returned element
	 */
	int Get( ElemType *pElem );
	/**
	 * @brief Returns the node's element type
	 */
	ElemType *GetNode();
	/**
	 * @brief Returns the last element on the list
	 * @param *pElem contains the last element
	 */
	int GetLast( ElemType *pElem );
	
private:
	/**
	 * @brief Pointer to head of the list
	 */
	Node *fHead;
	/**
	 * @brief Pointer to last element on list
	 */
	Node *fLast;
	/**
	 * @brief Pointer to current element
	 */
	Node *fThis;
	/**
	 * What the hell
	 */
	Node **fpBefore;
#ifdef WORKPOOL
	/**
	 * @brief Pointer to workpool
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
