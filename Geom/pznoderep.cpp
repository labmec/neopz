/**
 * @file
 * @brief Creates TPZNodeRep classes for several master elements. 
 */

#ifdef BORLAND
#include "pznoderep.h"
#endif
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpzpyramid.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

using namespace pztopology;
namespace pzgeom{

template<int N, class Topology>
bool TPZNodeRep<N,Topology>::IsLinearMapping() const
{
    return true;
}


template class TPZNodeRep<1,TPZPoint>;
template class TPZNodeRep<2,TPZLine>;
template class TPZNodeRep<3,TPZTriangle>;
template class TPZNodeRep<4,TPZQuadrilateral>;
template class TPZNodeRep<5,TPZPyramid>;
template class TPZNodeRep<4,TPZTetrahedron>;
template class TPZNodeRep<6,TPZPrism>;
template class TPZNodeRep<8,TPZCube>;

} // end of namespace pzgeom



