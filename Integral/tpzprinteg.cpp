/**
 * @file
 * @brief Contains the implementation of the TPZPrInteg methods. 
 */

#include "tpzprinteg.h"
#include "pzquad.h"

template<class TFather>
TPZPrInteg<TFather>::~TPZPrInteg()
{
}

template<class TFather>
TPZIntPoints *TPZPrInteg<TFather>::PrismExtend(int /*order*/)
{
	std::cout << "Please implement me " << __PRETTY_FUNCTION__;
	return 0;
}

#ifndef BORLAND
template<>
TPZIntPoints *TPZPrInteg<TPZInt1Point>::PrismExtend(int order)
{
	return new TPZPrInteg<TPZPrInteg<TPZInt1Point> >(order);
}

template<>
TPZIntPoints *TPZPrInteg< TPZPrInteg<TPZInt1Point> >::PrismExtend(int order)
{
	return new TPZPrInteg<TPZPrInteg<TPZPrInteg<TPZInt1Point> > >(order);
}
#endif

template class TPZPrInteg<TPZInt1Point>;
template class TPZPrInteg<TPZInt1d>;
template class TPZPrInteg<TPZIntQuad>;
template class TPZPrInteg<TPZIntTriang>;
template class TPZPrInteg<TPZIntCube3D>;
template class TPZPrInteg<TPZIntTetra3D>;
template class TPZPrInteg<TPZIntPyram3D>;
template class TPZPrInteg<TPZIntPrism3D>;

