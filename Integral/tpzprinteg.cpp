//
// C++ Implementation: tpzprinteg
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzprinteg.h"
#include "tpzint1point.h"
#include "pzquad.h"

template<class TFather>
TPZPrInteg<TFather>::~TPZPrInteg()
{
}

template<class TFather>
    inline TPZIntPoints *TPZPrInteg<TFather>::PrismExtend(int order)
{
  std::cout << "Please implement me " << __PRETTY_FUNCTION__;
  return 0;
}

template<>
    inline TPZIntPoints *TPZPrInteg<TPZInt1Point>::PrismExtend(int order)
{
  return new TPZPrInteg<TPZPrInteg<TPZInt1Point> >(order);
}

template<>
    inline TPZIntPoints *TPZPrInteg< TPZPrInteg<TPZInt1Point> >::PrismExtend(int order)
{
  return new TPZPrInteg<TPZPrInteg<TPZPrInteg<TPZInt1Point> > >(order);
}

template class TPZPrInteg<TPZInt1d>;
template class TPZPrInteg<TPZIntTriang>;
template class TPZPrInteg<TPZIntQuad>;
template class TPZPrInteg<TPZIntCube3D>;
template class TPZPrInteg<TPZIntTetra3D>;
template class TPZPrInteg<TPZIntPrism3D>;
template class TPZPrInteg<TPZIntPyram3D>;

