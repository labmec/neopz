// This may look like C code, but it is really -*- C++ -*-
/*
  Copyright (C) 1988 Free Software Foundation
  written by Doug Lea (dl@rocky.oswego.edu)

  This file is part of the GNU C++ Library.  This library is free
  software; you can redistribute it and/or modify it under the terms of
  the GNU Library General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your
  option) any later version.  This library is distributed in the hope
  that it will be useful, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the GNU Library General Public License for more details.
  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free Software
  Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/


#ifndef _TIndexTObjMap_h

#ifdef __GNUG__
#pragma interface
#endif
#define _TIndexTObjMap_h 1

#include "pzpix.h"
//#include "longdefs.h"
//#include "voiddefs.h"
//#include "gnudefs.h"

class TReceiveStorage;

template <class TIndex, class TObj>
class TMap
{
   protected:
      int  count;
      TObj def;

   public:
      /** Constructor for the map passing a default object to it. */
      TMap(TObj dflt);
      virtual ~TMap();

      int  length();                // current number of items
      int  empty();

      virtual int  contains(TIndex  key);   // is key mapped?
      virtual void clear();                 // delete all items

      virtual TObj& operator [] (TIndex  key) = 0; // access contents by key

      virtual void    del(TIndex  key) = 0;    // delete entry
      virtual TPZPix  first() = 0;             // TPZPix of first item or 0
      virtual void    next(TPZPix& i) = 0;     // advance to next or 0
      virtual TIndex& key(TPZPix i) = 0;       // access key at i
      virtual TObj&   contents(TPZPix i) = 0;  // access contents at i

      virtual int    owns(TPZPix i);             // is i a valid TPZPix  ?
      virtual TPZPix seek(TIndex  key);          // TPZPix of key

      TObj& dflt();    // access default val

      void  error(const char* msg);

      virtual int OK() = 0;   // rep invariant
};

template<class TIndex, class TObj>
inline TMap<TIndex,TObj>::~TMap() {}

template<class TIndex, class TObj>
inline int TMap<TIndex,TObj>::length()
{
   return count;
}

template<class TIndex, class TObj>
inline int TMap<TIndex,TObj>::empty()
{
   return count == 0;
}

template<class TIndex, class TObj>
inline TObj& TMap<TIndex,TObj>::dflt()
{
   return def;
}

template<class TIndex, class TObj>
inline TMap<TIndex,TObj>::TMap(TObj  dflt) :def(dflt)
{
   count = 0;
}

#endif

//--| PZ |----------------------------------------------------------------------
