//
//  pzgeom_utility.c
//  AcoplamentoH1Hdiv
//
//  Created by Philippe Devloo on 31/08/19.
//

#include "pzgeom_utility.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeom_util");
#endif


using namespace pzgeom;


using namespace pztopology;

template bool IsInSideParametricDomain<TPZPoint>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZLine>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZQuadrilateral>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZTriangle>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZPyramid>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZTetrahedron>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZPrism>(int side, const TPZVec<REAL> &pt, REAL tol);
template bool IsInSideParametricDomain<TPZCube>(int side, const TPZVec<REAL> &pt, REAL tol);


template void GetSideShapeFunction<TPZPoint>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZLine>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZTriangle>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZQuadrilateral>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZPyramid>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZTetrahedron>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZPrism>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void GetSideShapeFunction<TPZCube>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);

#include "fad.h"
//template<class T=REAL>
//class Fad;
template void GetSideShapeFunction<TPZPoint>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZLine>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZTriangle>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZQuadrilateral>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZPyramid>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZTetrahedron>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZPrism>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void GetSideShapeFunction<TPZCube>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);

