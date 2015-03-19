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

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoquad"));
#endif

TPZFMatrix<REAL> TensorProd(TPZFMatrix<REAL> &mat1, TPZFMatrix<REAL> &mat2);

namespace pzgeom {
	
    template<class GeomQuad>
	TPZGeoEl *TPZQuadSphere<GeomQuad>::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
	{
    
		int ns = orig->NSideNodes(side);
		TPZManVector<long> nodeindices(ns);
		int in;
		for(in=0; in<ns; in++)
		{
			nodeindices[in] = orig->SideNodeIndex(side,in);
		}
		long index;
		
		TPZGeoMesh *mesh = orig->Mesh();
		MElementType type = orig->Type(side);
		
		TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
		TPZGeoElSide me(orig,side);
		TPZGeoElSide newelside(newel,newel->NSides()-1);
		
		newelside.InsertConnectivity(me);
		newel->Initialize();
		
		return newel;
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	/** @brief Creates a geometric element according to the type of the father element */
    template<class GeomQuad>
TPZGeoEl *TPZQuadSphere<GeomQuad>::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
																						TPZVec<long>& nodeindexes,
																						int matid,
																						long& index)
	
	{
		return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
	}
	
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




/**
 * @ingroup geometry
 * @brief Id for three dimensional arc element
 */

template<>
int TPZGeoElRefPattern<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad> >::ClassId() const {
	return TPZGEOELEMENTQUADSPHEREID;
}

template class pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad>;

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad> >, TPZGEOELEMENTQUADSPHEREID>;


template class TPZGeoElRefLess<pzgeom::TPZQuadSphere<pzgeom::TPZGeoQuad> >;


template<>
int TPZGeoElRefPattern<pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > >::ClassId() const {
    return TPZGEOELEMENTQUADSPHEREBLENDID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadSphere< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > >, TPZGEOELEMENTQUADSPHEREBLENDID>;


template class TPZGeoElRefLess<pzgeom::TPZQuadSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > >;
