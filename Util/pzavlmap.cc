
#include "pzreal.h"
#include <iostream>
using namespace std;
#include "pzavlmap.h"

/*
 traversal primitives
*/

template<class TIndex,class TObj>
TPZAVLMap<TIndex,TObj>::TPZAVLMap(TPZAVLMap<TIndex,TObj> &map) {
  fNItems = 0;
  fDefault = map.fDefault;
  fRoot = 0;
  for(TPZPix i = map.First(); i != 0; map.Next(i))
    (*this)[map.Key(i)] = map.Contents(i);
}

template<class TIndex,class TObj>
TPZPix TPZAVLMap<TIndex,TObj>::First() {
  TPZAVLNode<TIndex,TObj>* t = fRoot;
  if (t != 0)
    while(t->fLeftNode != 0)
      t = t->fLeftNode;

  return TPZPix(t);
}

template<class TIndex,class TObj>
TPZPix TPZAVLMap<TIndex,TObj>::Last() {
  TPZAVLNode<TIndex,TObj>* t = fRoot;
  if (t != 0)
    while(t->fRightNode != 0)
      t = t->fRightNode;

  return TPZPix(t);
}
//até todo feito
template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>::Next(TPZPix &i) {
  if(i != 0) {
    TPZAVLNode<TIndex,TObj>* r = ((TPZAVLNode<TIndex,TObj>*)i)->fRightNode;
    //if (!rthread(r))
    if (!rthread((TPZAVLNode<TIndex,TObj>*)i))//corrigido por Cedric 27/10/98 : 18:30
      while (!lthread(r))
	r = r->fLeftNode;
    i = TPZPix(r);
  }
}
/*
inline void LongIntAVLMap::prev(Pix& i)
{
  if (i != 0) i = Pix(pred((LongIntAVLNode*)i));
}

*/
template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>::Previous(TPZPix &i) {
  if(i != 0) {
    TPZAVLNode<TIndex,TObj>* l = ((TPZAVLNode<TIndex,TObj>*)i)->fLeftNode;
    //if(!lthread(l))
    if(!lthread((TPZAVLNode<TIndex,TObj>*)i))//corrigido por Cedric 27/10/98 : 18:30
      while (!rthread(l))
	l = l->fRightNode;
    i = TPZPix(l);
  }
}

template<class TIndex,class TObj>
TPZPix TPZAVLMap<TIndex,TObj>::Seek(TIndex  key) {
  TPZAVLNode<TIndex,TObj>* t = fRoot;
  if (t == 0)  return 0;
  for (;;) {
    //		int cmp = (key < t->fIndex) ? -1 : ((key > t->fIndex) ? 1 : 0);
    if(key == t->fIndex)  return TPZPix(t);
    else if(key < t->fIndex) {
      if(lthread(t))  return 0;
      else  t = t->fLeftNode;
    }
    else {
      if(rthread(t)) return 0;
      else  t = t->fRightNode;
    }
  }
}


/*
 The combination of threads and AVL bits make adding & deleting
 interesting, but very awkward.

 We use the following statics to avoid passing them around recursively
*/

static int _need_rebalancing;   // to send back balance info from rec. calls

template<class TIndex,class TObj>
TIndex *TPZAVLMap<TIndex,TObj>::_target_item; // add/del_item target


template<class TIndex,class TObj>
TPZAVLNode<TIndex,TObj> *TPZAVLMap<TIndex,TObj>::_found_node; // returned added/deleted node

static int    _already_found;   // for deletion subcases


template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>:: _add(TPZAVLNode<TIndex,TObj> *&t) {
  if(*_target_item == t->fIndex) {
    _found_node = t;
    return;
  }
  else if(*_target_item < t->fIndex) {
    if(lthread(t)) {
      ++fNItems;
      _found_node = new TPZAVLNode<TIndex,TObj>(*_target_item,fDefault);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->fLeftNode = t->fLeftNode;
      _found_node->fRightNode = t;
      t->fLeftNode = _found_node;
      set_lthread(t, 0);
      _need_rebalancing = 1;
    }
    else  _add(t->fLeftNode);

    if(_need_rebalancing) {
      switch(bf(t)) {
      case AVLRIGHTHEAVY:
	set_bf(t, AVLBALANCED);
	_need_rebalancing = 0;
	return;
      case AVLBALANCED:
	set_bf(t, AVLLEFTHEAVY);
	return;
      case AVLLEFTHEAVY:
	{
	  TPZAVLNode<TIndex,TObj>* l = t->fLeftNode;
	  if(bf(l) == AVLLEFTHEAVY) {
	    if(rthread(l))  t->fLeftNode = l;
	    else  t->fLeftNode = l->fRightNode;
	    set_lthread(t, rthread(l));
	    l->fRightNode = t;
	    set_rthread(l, 0);
	    set_bf(t, AVLBALANCED);
	    set_bf(l, AVLBALANCED);
	    t = l;
	    _need_rebalancing = 0;
	  }
	  else {
	    TPZAVLNode<TIndex,TObj>* r = l->fRightNode;
	    set_rthread(l, lthread(r));

	    if(lthread(r)) l->fRightNode = r;
	    else  l->fRightNode = r->fLeftNode;
	    r->fLeftNode = l;
	    set_lthread(r, 0);
	    set_lthread(t, rthread(r));
	    if(rthread(r))  t->fLeftNode = r;
	    else  t->fLeftNode = r->fRightNode;
	    r->fRightNode = t;
	    set_rthread(r, 0);

	    if(bf(r) == AVLLEFTHEAVY)  set_bf(t, AVLRIGHTHEAVY);
	    else  set_bf(t, AVLBALANCED);
	    if(bf(r) == AVLRIGHTHEAVY)  set_bf(l, AVLLEFTHEAVY);
	    else  set_bf(l, AVLBALANCED);

	    set_bf(r, AVLBALANCED);
	    t = r;
	    _need_rebalancing = 0;
	    return;
	  }
	}
      }
    }
  }
  else {
    if (rthread(t)) {
      ++fNItems;
      _found_node = new TPZAVLNode<TIndex,TObj>(*_target_item,fDefault);
      set_rthread(t, 0);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->fLeftNode = t;
      _found_node->fRightNode = t->fRightNode;
      t->fRightNode = _found_node;
      _need_rebalancing = 1;
    }
    else  _add(t->fRightNode);

    if(_need_rebalancing) {
      switch(bf(t)) {
      case AVLLEFTHEAVY:
	set_bf(t, AVLBALANCED);
	_need_rebalancing = 0;
	return;
      case AVLBALANCED:
	set_bf(t, AVLRIGHTHEAVY);
	return;
      case AVLRIGHTHEAVY:
	{
	  TPZAVLNode<TIndex,TObj>* r = t->fRightNode;
	  if(bf(r) == AVLRIGHTHEAVY) {
	    if (lthread(r))  t->fRightNode = r;
	    else  t->fRightNode = r->fLeftNode;
	    set_rthread(t, lthread(r));
	    r->fLeftNode = t;
	    set_lthread(r, 0);
	    set_bf(t, AVLBALANCED);
	    set_bf(r, AVLBALANCED);
	    t = r;
	    _need_rebalancing = 0;
	  }
	  else {
	    TPZAVLNode<TIndex,TObj>* l = r->fLeftNode;
	    set_lthread(r, rthread(l));
	    if(rthread(l))  r->fLeftNode = l;
	    else  r->fLeftNode = l->fRightNode;
	    l->fRightNode = r;
	    set_rthread(l, 0);
	    set_rthread(t, lthread(l));
	    if(lthread(l))  t->fRightNode = l;
	    else  t->fRightNode = l->fLeftNode;
	    l->fLeftNode = t;
	    set_lthread(l, 0);

	    if(bf(l) == AVLRIGHTHEAVY)  set_bf(t, AVLLEFTHEAVY);
	    else  set_bf(t, AVLBALANCED);
	    if(bf(l) == AVLLEFTHEAVY) set_bf(r, AVLRIGHTHEAVY);
	    else  set_bf(r, AVLBALANCED);
	    set_bf(l, AVLBALANCED);
	    t = l;
	    _need_rebalancing = 0;
	    return;
	  }
	}
      }
    }
  }
}

template<class TIndex,class TObj>
TObj &TPZAVLMap<TIndex,TObj>::operator [] (TIndex item) {
  if (fRoot == 0) {
    ++fNItems;
    fRoot = new TPZAVLNode<TIndex,TObj>(item,fDefault);
    set_rthread(fRoot, 1);
    set_lthread(fRoot, 1);
    return fRoot->fContent;
  }
  else {
    _target_item = &item;
    _need_rebalancing = 0;
    _add(fRoot);
    return _found_node->fContent;
  }
}

template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>::_del(TPZAVLNode<TIndex,TObj> *par,TPZAVLNode<TIndex,TObj> *&t) {
  int COMP;
  if(_already_found) {
    if(rthread(t))
      COMP = 0;
    else
      COMP = 1;
  }
  else
    COMP = (*_target_item < t->fIndex) ? -1 : ((*_target_item > t->fIndex) ? 1 : 0);
  if (COMP == 0) {
    if (lthread(t) && rthread(t)) {
      _found_node = t;
      if(t == par->fLeftNode) {
        set_lthread(par, 1);
        par->fLeftNode = t->fLeftNode;
      }
      else  {
        set_rthread(par, 1);
        par->fRightNode = t->fRightNode;
      }
      _need_rebalancing = 1;
      return;
    }
    else if (lthread(t)) {
      _found_node = t;

      TPZAVLNode<TIndex,TObj> *s = t->fRightNode;
      if(!rthread(t)) while(!lthread(s)) s = s->fLeftNode;

      if(s != 0 && lthread(s))  s->fLeftNode = t->fLeftNode;
      t = t->fRightNode;
      _need_rebalancing = 1;
      return;
    }
    else if(rthread(t))  {
      _found_node = t;

      TPZAVLNode<TIndex,TObj>* p = t->fLeftNode;
      if(!lthread(t)) while(!rthread(p)) p = p->fRightNode;

      if(p != 0 && rthread(p))  p->fRightNode = t->fRightNode;
      t = t->fLeftNode;
      _need_rebalancing = 1;
      return;
    }
    else {                       // replace item & find someone deletable
      TPZAVLNode<TIndex,TObj>* p = t->fLeftNode;
      if(!lthread(t)) while(!rthread(p)) p = p->fRightNode;

      t->fIndex = p->fIndex;
      t->fContent = p->fContent;
      _already_found = 1;
      COMP = -1;                // fall through below to left
    }
  }

  if(COMP < 0) {
    if(lthread(t))  return;
    _del(t, t->fLeftNode);
    if (!_need_rebalancing) return;
    switch (bf(t)) {
    case AVLLEFTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLRIGHTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLRIGHTHEAVY:
      {
      	TPZAVLNode<TIndex,TObj>* r = t->fRightNode;
	switch (bf(r))  {
	case AVLBALANCED:
	  if(lthread(r))  t->fRightNode = r;
	  else  t->fRightNode = r->fLeftNode;
	  set_rthread(t, lthread(r));
	  r->fLeftNode = t;
	  set_lthread(r, 0);
	  set_bf(t, AVLRIGHTHEAVY);
	  set_bf(r, AVLLEFTHEAVY);
	  _need_rebalancing = 0;
	  t = r;
	  return;
	case AVLRIGHTHEAVY:
	  if(lthread(r))  t->fRightNode = r;
	  else  t->fRightNode = r->fLeftNode;
	  set_rthread(t, lthread(r));
	  r->fLeftNode = t;
	  set_lthread(r, 0);
	  set_bf(t, AVLBALANCED);
	  set_bf(r, AVLBALANCED);
	  t = r;
	  return;
	case AVLLEFTHEAVY:
	  {
	    TPZAVLNode<TIndex,TObj>* l = r->fLeftNode;
	    set_lthread(r, rthread(l));
	    if(rthread(l))  r->fLeftNode = l;
	    else  r->fLeftNode = l->fRightNode;
	    l->fRightNode = r;
	    set_rthread(l, 0);
	    set_rthread(t, lthread(l));
	    if(lthread(l)) t->fRightNode = l;
	    else  t->fRightNode = l->fLeftNode;
	    l->fLeftNode = t;
	    set_lthread(l, 0);
	    if(bf(l) == AVLRIGHTHEAVY) set_bf(t, AVLLEFTHEAVY);
	    else  set_bf(t, AVLBALANCED);
	    if(bf(l) == AVLLEFTHEAVY) set_bf(r, AVLRIGHTHEAVY);
	    else set_bf(r, AVLBALANCED);
	    set_bf(l, AVLBALANCED);
	    t = l;
	    return;
	  }
	}
      }
    }
  }
  else {
    if(rthread(t))  return;
    _del(t, t->fRightNode);
    if(!_need_rebalancing) return;

    switch(bf(t)) {
    case AVLRIGHTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLLEFTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLLEFTHEAVY:
      {
      	TPZAVLNode<TIndex,TObj>* l = t->fLeftNode;
      	switch (bf(l))
	  {
	  case AVLBALANCED:
	    if(rthread(l))  t->fLeftNode = l;
	    else  t->fLeftNode = l->fRightNode;
	    set_lthread(t,rthread(l));
	    l->fRightNode = t;
	    set_rthread(l, 0);
	    set_bf(t, AVLLEFTHEAVY);
	    set_bf(l, AVLRIGHTHEAVY);
	    _need_rebalancing = 0;
	    t = l;
	    return;
	  case AVLLEFTHEAVY:
	    if(rthread(l)) t->fLeftNode = l;
	    else  t->fLeftNode = l->fRightNode;
	    set_lthread(t, rthread(l));
	    l->fRightNode = t;
	    set_rthread(l, 0);
	    set_bf(t, AVLBALANCED);
	    set_bf(l, AVLBALANCED);
	    t = l;
	    return;
	  case AVLRIGHTHEAVY:
            {
	      TPZAVLNode<TIndex,TObj>* r = l->fRightNode;
	      set_rthread(l, lthread(r));
	      if(lthread(r)) l->fRightNode = r;
	      else  l->fRightNode = r->fLeftNode;
	      r->fLeftNode = l;
	      set_lthread(r, 0);
	      set_lthread(t, rthread(r));
	      if(rthread(r))  t->fLeftNode = r;
	      else t->fLeftNode = r->fRightNode;
	      r->fRightNode = t;
	      set_rthread(r, 0);
	      if(bf(r) == AVLLEFTHEAVY) set_bf(t, AVLRIGHTHEAVY);
	      else  set_bf(t, AVLBALANCED);
	      if(bf(r) == AVLRIGHTHEAVY) set_bf(l, AVLLEFTHEAVY);
	      else  set_bf(l, AVLBALANCED);
	      set_bf(r, AVLBALANCED);
	      t = r;
	      return;
	    }
	  }
      }
    }
  }
}

template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>::Delete(TIndex item) {
  if(fRoot == 0) return;
  _need_rebalancing = 0;
  _already_found = 0;
  _found_node = 0;
  _target_item = &item;
  _del(fRoot, fRoot);
  if(_found_node) {
    delete(_found_node);
    if(--fNItems == 0)  fRoot = 0;
  }
}

template<class TIndex,class TObj>
void TPZAVLMap<TIndex,TObj>::_kill(TPZAVLNode<TIndex,TObj> *t) {
  if(t != 0) {
    if(!lthread(t)) _kill(t->fLeftNode);
    if (!rthread(t)) _kill(t->fRightNode);
    delete t;
  }
}


template class TPZAVLMap<int, void *>;

template class TPZAVLMap<int,int>;

template class TPZAVLMap<int,int *>;

template class TPZAVLMap<int,REAL>;

template class TPZAVLMap<int,REAL *>;

class TPZGeoEl;
template class TPZAVLMap<int,TPZGeoEl *>;
template class TPZAVLMap<TPZGeoEl *, TPZGeoEl *>;
