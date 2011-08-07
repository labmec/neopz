/**
 * @file
 * @brief Contains the implementation of the TPZLink methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tlink.cc
//
// Class:  TPZLink<class>
//
// Obs.:   Implementa uma lista ligada com elementos indefinidos.
//
// Versao: 12 / 1994.
//


#include "pzlink.h"
#include "pzespmat.h"


/*** Constructor ***/

template< class ElemType >
TPZLink<ElemType>::TPZLink()
{
	fHead = fThis = fLast = NULL;
	fpBefore = &fHead;
#ifdef WORKPOOL
	fWp = 0;
#endif
}

/*** Destructor ***/

template< class ElemType >
TPZLink<ElemType>::~TPZLink()
{
#ifdef WORKPOOL
	if(!fWp) {
#endif
		fThis = fHead;
		while ( fHead != NULL )
		{
			fHead = fHead->next;
			delete( fThis );
			fThis = fHead;
		}
#ifdef WORKPOOL
		//  } else {
		//	Clear();
	}
#endif
}



/*** Insert ***/

template< class ElemType >
int
TPZLink<ElemType>::Insert( ElemType &elem )
{
	Node *newNode;
#ifdef WORKPOOL
	if(fWp) newNode = (Node *) fWp->NewPointer(sizeof(Node));
	else newNode = new (Node);
#else
	newNode = new( Node );
#endif
	if ( newNode == NULL )
		return( 0 );
	
	newNode->elem = elem;
	if ( (newNode->next = fThis) == NULL )
		fLast = newNode;
	
	*fpBefore = fThis = newNode;
	
	return( 1 );
}



/*** Append ***/

template< class ElemType >
int
TPZLink<ElemType>::Append( ElemType &elem )
{
	// Vai para o final da lista.
	while ( fThis != NULL )
    {
		fpBefore = &fThis->next;
		fThis    =  fThis->next;
    }
	
	return( Insert( elem ) );
}



/*** Remove ***/

template< class ElemType >
int
TPZLink<ElemType>::Remove()
{
	if ( fThis == NULL )
		return( 0 );
	
	fThis = fThis->next;
#ifdef WORKPOOL
	if(!fWp) delete( *fpBefore );
	else fWp->Release(*fpBefore);
#else
	delete( *fpBefore );
#endif
	*fpBefore = fThis;
	
	// Se removeu o ultimo elemento, ajusta o ponteiro fLast.
	//
	if ( fThis == NULL )
    {
		Head();
		if ( fThis == NULL )
			fLast = NULL;
		else
		{
			while ( fThis->next != NULL )
				fThis = fThis->next;
			fLast = fThis;
			fThis = NULL;
		}
    }
	
	return( 1 );
}



/*** Update ***/

template< class ElemType >
int TPZLink<ElemType>::Update( ElemType &elem )
{
	if ( fThis == NULL )
		return( 0 );
	
	fThis->elem = elem;
	return( 1 );
}



/*** Clear ***/

template< class ElemType >
int
TPZLink<ElemType>::Clear()
{
	fThis = fHead;
	while ( fThis != NULL )
    {
		fThis = fThis->next;
#ifdef WORKPOOL
		if(fWp) fWp->Release(fHead);
		else delete fHead;
#else
		delete( fHead );
#endif
		fHead = fThis;
    }
	fpBefore = &fHead;
	fThis = fLast = NULL;
	return( 1 );
}



/*** Operator = ***/

template< class ElemType >
TPZLink<ElemType> &
TPZLink<ElemType>::operator=( TPZLink<ElemType> &src )
{
	// Guarda a posicao atual da lista 'src'.
	Node **oldBefore = src.fpBefore;
	
	Clear();
	
	// Copia os elementos da lista.
	src.Head();
	ElemType elem;
	while ( src.Get( &elem ) )
    {
		Append( elem );
		src.Next();
    }
	
	// Volta a lista 'src' aa posicao inicial.
	src.fpBefore = oldBefore;
	src.fThis    = *oldBefore;
	
	Head();
	return( *this );
}



/*********************** Tipos usados *************************/


template class TPZLink<TPZSpMatrix::TPZNode>;

