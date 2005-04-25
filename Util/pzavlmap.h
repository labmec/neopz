/**
 * @file pzavlmap.h
 * @brief Balanced AVL tree.
 */
// $Id: pzavlmap.h,v 1.2 2005-04-25 02:55:51 phil Exp $

#ifndef PZAVLMAP_H
#define PZAVLMAP_H

#include "pzpix.h"
#include "pzerror.h"


// Constants & inlines for maintaining balance & thread status in tree nodes.

#define AVLBALANCEMASK    3
#define AVLBALANCED       0
#define AVLLEFTHEAVY      1
#define AVLRIGHTHEAVY     2

#define LTHREADBIT        4
#define RTHREADBIT        8

//--| TPZAVLNode |--------------------------------------------------------------

/**
 * @ingroup util
 * @brief Implements a node of a binary tree of objects.
 *
 * This class is a templated extension of g++ library package.
 * It implements the node within a balanced AVL binary tree.
 */
template< class TIndex, class TObj >
class TPZAVLNode
{
   public:
      /** Index of the current item. */
      TIndex fIndex;

      /** Next item at the left of the current item. */
      TPZAVLNode< TIndex, TObj >* fLeftNode;

      /** Next item at the rigth of the current item. */
      TPZAVLNode< TIndex, TObj >* fRightNode;

      /** Current item. */
      TObj fContent;

      /** I have no idea what this is for. */
      char fStat;

      /** Constructor of the node. */
      TPZAVLNode( TIndex i, TObj &c,
		  TPZAVLNode< TIndex, TObj >* l = 0,
		  TPZAVLNode< TIndex, TObj >* r = 0);

      /** Destructor of the node. */
      ~TPZAVLNode();
};

template< class TIndex, class TObj >
inline TPZAVLNode< TIndex, TObj >::
TPZAVLNode( TIndex index, TObj &c,
	    TPZAVLNode< TIndex, TObj >* l,
	    TPZAVLNode< TIndex, TObj >* r )
   : fIndex( index ), fLeftNode( l ), fRightNode( r ),
     fContent( c ), fStat( 0 )
{
   // NOTHING TO DO HERE!
}

template< class TIndex, class TObj >
inline TPZAVLNode< TIndex, TObj >::~TPZAVLNode()
{
   // NOTHING TO DO HERE!
}

//--| TPZAVLMap |---------------------------------------------------------------

/**
 * @ingroup util
 * @brief Implements a binary tree of objects.
 *
 * This class is a templated extension of g++ library package.
 * It implements a balanced AVL binary tree.
 */
template< class TIndex, class TObj >
class TPZAVLMap
{
   protected:
      /** Number of elements in the tree. */
      int fNItems;

      /**
       * Default object to be returned in case of a search of a non
       * existing element.
       */
      TObj fDefault;

      /** This indicates that this class is not thread_safe? */
      static TPZAVLNode< TIndex, TObj >* _found_node;

      /**
       * Add/del item target @note this indicates that this class is
       * not thread safe?
       */
      static TIndex* _target_item;

   public:
      /** Root node of the map. */
      TPZAVLNode< TIndex, TObj >* fRoot;

      /** Kill node t at the map. */
      void _kill( TPZAVLNode< TIndex, TObj >* t );

      /** Adding at the map node t. */
      void _add( TPZAVLNode< TIndex, TObj >*& t );

      /** Delete node t into node p. */
      void _del( TPZAVLNode< TIndex, TObj >* p,
		 TPZAVLNode< TIndex, TObj >*& t );

   public:
      /** Constructor for the map passing a default object to it. */
      TPZAVLMap( TObj &dflt );

      /** Copy constructor. */
      TPZAVLMap( TPZAVLMap &map );

      /** Destructor. */
      ~TPZAVLMap();

      /** Current number of items. */
      int NItems();

      /** Return TPZPix of the first item or zero. */
      TPZPix First();

      /** Return TPZPix of the last item or zero. */
      TPZPix Last();

      /** Advanced to next item. */
      void Next( TPZPix& i );

      /** Turn to previous item. */
      void Previous( TPZPix& i );

      /** Return the index of the item associated with i. */
      TIndex& Key( TPZPix i );

      /** Access the contents of the i. */
      TObj& Contents( TPZPix i );

      /** Return pix of the item with index key. */
      TPZPix Seek( TIndex key );

      /** Verifies if an item already uses the key. */
      int Contains( TIndex  key ) { return Seek( key ) != 0; }

      /** Access object contented with index key. */
      TObj& operator [] ( TIndex key );

      /** Delete object with index key. */
      void Delete( TIndex  key );

      /** Delete all items of the map. */
      void CleanUp();

      /** Access to the default value. */
      TObj& Default();
};

template< class TIndex, class TObj >
inline TPZAVLMap< TIndex, TObj >::TPZAVLMap( TObj &dflt )
   : fNItems( 0 ), fDefault( dflt ), fRoot( 0 )
{
   // NOTHING TO DO HERE!
}

template< class TIndex, class TObj >
inline TPZAVLMap< TIndex, TObj >::~TPZAVLMap()
{
   _kill( fRoot );
}

template< class TIndex, class TObj >
inline int TPZAVLMap< TIndex, TObj >::NItems()
{
   return fNItems;
}

template< class TIndex, class TObj >
inline TObj &TPZAVLMap< TIndex, TObj >::Default()
{
   return fDefault;
}

template< class TIndex, class TObj >
inline TIndex &TPZAVLMap< TIndex, TObj >::Key( TPZPix i )
{
   if( i == 0 )
   {
      PZError << "Null TPZPix to TPZAVLMap::Key";
   }

   return ( static_cast< TPZAVLNode< TIndex, TObj > * >( i ) )->fIndex;
}

template< class TIndex, class TObj >
inline TObj& TPZAVLMap< TIndex, TObj >::Contents( TPZPix i )
{
   if( i == 0 )
   {
      PZError << "Null TPZPix to TPZAVLMap::Contents";
   }

   return ( static_cast< TPZAVLNode< TIndex, TObj > * >( i ) )->fContent;
}

template< class TIndex, class TObj >
inline void TPZAVLMap< TIndex, TObj >::CleanUp()
{
   _kill( fRoot );

   fNItems = 0;
   fRoot = 0;
}

typedef TPZAVLMap< int, void * > TPZVoidPtrMap;

//--| Out of class |------------------------------------------------------------

template<class TIndex,class TObj>
static inline int bf(TPZAVLNode<TIndex,TObj> *t) {
  return t->fStat &AVLBALANCEMASK;
}

template<class TIndex,class TObj>
static inline void set_bf(TPZAVLNode<TIndex,TObj> *t, char b) {
  t->fStat = (char)((t->fStat & ~AVLBALANCEMASK) | (b & AVLBALANCEMASK));
}

template<class TIndex,class TObj>
static inline int rthread(TPZAVLNode<TIndex,TObj> *t) {
  return t->fStat & RTHREADBIT;
}

template<class TIndex,class TObj>
static inline void set_rthread(TPZAVLNode<TIndex,TObj> *t, int b) {
  if(b)  t->fStat |= RTHREADBIT;
  else  t->fStat &= ~RTHREADBIT;
}

template<class TIndex,class TObj>
static inline int lthread(TPZAVLNode<TIndex,TObj> *t) {
  return t->fStat & LTHREADBIT;
}

template<class TIndex,class TObj>
static inline void set_lthread(TPZAVLNode<TIndex,TObj> *t, int b) {
  if(b)  t->fStat |= LTHREADBIT;
  else  t->fStat &= ~LTHREADBIT;
}

#endif //PZAVLMAP_H

//--| PZ |----------------------------------------------------------------------
