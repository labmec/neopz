//
//  TPZQuadTorus.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZQuadTorus.h"
#include "tpzgeomid.h"

namespace pzgeom {

TPZGeoEl *TPZQuadTorus::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
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

#include "tpzgeoelmapped.h"
    
    /**
     * Creates a geometric element according to the type of the father element
     */
    /** @brief Creates a geometric element according to the type of the father element */
    TPZGeoEl *TPZQuadTorus::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                      TPZVec<long>& nodeindexes,
                                      int matid,
                                      long& index)

    {
        return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    }

    

}

#include "tpzgeoelrefpattern.h"

/**
 * @ingroup geometry
 * @brief Id for three dimensional arc element
 */

template<>
int TPZGeoElRefPattern<pzgeom::TPZQuadTorus>::ClassId() const {
	return TPZGEOELEMENTQUADTORUSID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadTorus>, TPZGEOELEMENTQUADTORUSID>;


template class TPZGeoElRefLess<pzgeom::TPZQuadTorus>;

