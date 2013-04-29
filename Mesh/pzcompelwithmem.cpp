/**
 * @file
 * @brief Contains the declaration of the TPZCompElWithMem class, it is as TPZCompEl with enable material memory feature.
 */

#include "pzcompelwithmem.h"

template<class TBASE>
TPZCompElWithMem<TBASE>::~TPZCompElWithMem() {
	SetFreeIntPtIndices();
}


template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePoint> >, TPZCOMPELWITHMEMPOINTID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeLinear> >, TPZCOMPELWITHMEMLINEARID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTriang> >, TPZCOMPELWITHMEMTRIANGID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeQuad> >, TPZCOMPELWITHMEMQUADID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeCube> >, TPZCOMPELWITHMEMCUBEID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapeTetra> >, TPZCOMPELWITHMEMTETRAID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePrism> >, TPZCOMPELWITHMEMPRISMID>;
template class TPZRestoreClass<TPZCompElWithMem<TPZIntelGen<pzshape::TPZShapePiram> >, TPZCOMPELWITHMEMPIRAMID>;



