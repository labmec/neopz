//
//  TPZQuadSphere.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZQuadSphere.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeoquad");
#endif

TPZFMatrix<REAL> TensorProd(TPZFMatrix<REAL> &mat1, TPZFMatrix<REAL> &mat2);

namespace pzgeom {
	
	  template<class GeomQuad>
    int TPZQuadSphere<GeomQuad>::ClassId() const
    {
      return Hash("TPZQuadSphere") ^
        ClassIdOrHash<GeomQuad>() << 1;
    }
    template<class GeomQuad>
    void TPZQuadSphere<GeomQuad>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &sz)
    {
        TPZManVector<REAL,3> center(lowercorner);
        REAL radius = 1.;
        center[0] += radius/2.;
        center[1] += radius/2.;
        center[2] += radius/2.;
        for (int i=0; i<3; i++) {
            sz[i] = 2.*radius;
        }
        REAL coords[4][3] = {
            {-1,-1,-0.1},
            { 1,-1,-0.1},
            { 1, 1,-0.1},
            {-1, 1,-0.1}
        };
        for (int i=0; i<4; i++) {
            REAL norm = sqrt(coords[i][0]*coords[i][0]+coords[i][1]*coords[i][1]+coords[i][2]*coords[i][2]);
            for(int j=0; j<3; j++) coords[i][j] *= radius/norm;
        }
        TPZManVector<int64_t,4> indices(4);
        for (int i=0; i<4; i++) {
            indices[i] = gmesh.NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3);
            for (int j=0; j<3; j++) {
                xco[j] = coords[i][j]+center[j];
            }
            gmesh.NodeVec()[indices[i]].Initialize(xco, gmesh);
        }
        TPZGeoElRefPattern<TPZQuadSphere<> > *gel = new TPZGeoElRefPattern<TPZQuadSphere<> >(indices,matid,gmesh);
        gel->Geom().SetData(radius, center);
    }

// 	/**
// 	 * Creates a geometric element according to the type of the father element
// 	 */
// 	/** @brief Creates a geometric element according to the type of the father element */
//     template<class GeomQuad>
// TPZGeoEl *TPZQuadSphere<GeomQuad>::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
// 																						TPZVec<int64_t>& nodeindexes,
// 																						int matid,
// 																						int64_t& index)
	
// 	{
// 		return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
// 	}
	
    /** @brief declare geometry as blended element */
    template<class GeomQuad>
    bool TPZQuadSphere<GeomQuad>::IsGeoBlendEl() const
    {
        return false;
    }
	
    /** @brief declare geometry as blended element */
    template<>
    bool TPZQuadSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >::IsGeoBlendEl() const
    {
        return true;
    }
}

template class pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>;

template class pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad > >;

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad> >>;

template class TPZGeoElRefLess<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad> >;

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > >>;

template class TPZGeoElRefLess<pzgeom::TPZQuadSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > >;
