/** @file pzavlmapgen.cc */
// $Id: pzavlmapgen.cc,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#include "pzreal.h"
#include "pzavlmap.h"

class TPZGeoEl;

template class TPZAVLMap< int, void* >;
template class TPZAVLMap< int, int >;
template class TPZAVLMap< int, int* >;
template class TPZAVLMap< int, REAL >;
template class TPZAVLMap< int, REAL* >;
template class TPZAVLMap< int, TPZGeoEl* >;
template class TPZAVLMap< TPZGeoEl*, TPZGeoEl* >;

//--| PZ |----------------------------------------------------------------------
