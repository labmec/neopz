#include "TPZTopologyUtils.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzpyramid.h"
#include "tpzprism.h"

namespace pztopology{
    REAL GetTolerance(){return gTolerance;}

    void SetTolerance(const REAL &tol){
        if(tol > 0) gTolerance = tol;
        else {
            typedef std::numeric_limits< REAL > dbl;

            std::cout.precision(dbl::max_digits10);
            std::cout<<"Invalid tolerance parameter for topologies. Trying to set: "<<tol<<std::endl;
            std::cout<<"This value will be ignored. Tolerance is set at: "<<gTolerance<<std::endl;
        }
    }

    template<class Topology>
    void GetPermutation(const int& permute, TPZVec<int> &permutation){
        #ifdef PZDEBUG
        if(permute < 0 || permute >= Topology::NPermutations){
            PZError<<"GetPermutation: invalid parameter: permute = "<<permute<<std::endl;
            PZError<<"Aborting..."<<std::endl;
            DebugStop();
        }
        #endif
        permutation.Resize(Topology::NSides,-1);
        for(int i = 0; i < Topology::NSides; i++) permutation[i] = Topology::fPermutations[permute][i];
    }

}

template void pztopology::GetPermutation<pztopology::TPZPoint>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZLine>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZTriangle>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZQuadrilateral>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZTetrahedron>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZCube>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZPrism>(const int& permute, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZPyramid>(const int& permute, TPZVec<int> &permutation);



#define IMPLEMENTPERMUTATION(TGEO,TTOPOL) \
\
template<>\
void pztopology::GetPermutation<TGEO>(const int& permute, TPZVec<int> &permutation){\
    return pztopology::GetPermutation<TTOPOL>(permute,permutation);\
};


#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"
//linear elements
IMPLEMENTPERMUTATION(pzgeom::TPZGeoPoint,pztopology::TPZPoint)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoLinear,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoTriangle,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoQuad,pztopology::TPZQuadrilateral)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoTetrahedra,pztopology::TPZTetrahedron)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoCube,pztopology::TPZCube)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoPrism,pztopology::TPZPrism)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoPyramid,pztopology::TPZPyramid)

#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"
//quadratic elements
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticLine,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticTrig,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticQuad,pztopology::TPZQuadrilateral)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticTetra,pztopology::TPZTetrahedron)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticCube,pztopology::TPZCube)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticPrism,pztopology::TPZPrism)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadraticPyramid,pztopology::TPZPyramid)

#include "TPZWavyLine.h"
#include "tpzarc3d.h"
#include "tpzellipse3d.h"
#include "TPZTriangleSphere.h"
#include "TPZQuadSphere.h"
IMPLEMENTPERMUTATION(pzgeom::TPZWavyLine,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZArc3D,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZEllipse3D,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle>,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>,pztopology::TPZQuadrilateral)
#include "tpzgeoblend.h"

//blend elements (linear)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoPoint>,pztopology::TPZPoint)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoLinear>,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>,pztopology::TPZQuadrilateral)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>,pztopology::TPZTetrahedron)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>,pztopology::TPZCube)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism>,pztopology::TPZPrism)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZGeoPyramid>,pztopology::TPZPyramid)
//blend elements (quadratic)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticLine>,pztopology::TPZLine)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticTrig>,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticQuad>,pztopology::TPZQuadrilateral)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticTetra>,pztopology::TPZTetrahedron)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticCube>,pztopology::TPZCube)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticPrism>,pztopology::TPZPrism)
IMPLEMENTPERMUTATION(pzgeom::TPZGeoBlend<pzgeom::TPZQuadraticPyramid>,pztopology::TPZPyramid)

IMPLEMENTPERMUTATION(pzgeom::TPZTriangleSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>,pztopology::TPZTriangle)
IMPLEMENTPERMUTATION(pzgeom::TPZQuadSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>,pztopology::TPZQuadrilateral)

#undef IMPLEMENTPERMUTATION
