//
//  TPZTriangleTorus.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZTriangleTorus.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

namespace pzgeom {

TPZGeoEl *TPZTriangleTorus::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
{
    
        int ns = orig->NSideNodes(side);
        TPZManVector<int64_t> nodeindices(ns);
        int in;
        for(in=0; in<ns; in++)
        {
            nodeindices[in] = orig->SideNodeIndex(side,in);
        }
        int64_t index;
        
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
    TPZGeoEl *TPZTriangleTorus::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                      TPZVec<int64_t>& nodeindexes,
                                      int matid,
                                      int64_t& index)

    {
        return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    }

    int TPZTriangleTorus::ClassId() const{
        return Hash("TPZTriangleTorus") ^ TPZGeoTriangle::ClassId() << 1;
    }

}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZTriangleTorus>>;


template class TPZGeoElRefLess<pzgeom::TPZTriangleTorus>;

