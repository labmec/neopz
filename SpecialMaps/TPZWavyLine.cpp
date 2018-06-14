//
//  TPZWavyLine.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZWavyLine.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

namespace pzgeom {

TPZGeoEl *TPZWavyLine::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
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
    TPZGeoEl *TPZWavyLine::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                      TPZVec<int64_t>& nodeindexes,
                                      int matid,
                                      int64_t& index)

    {
        return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    }


    int TPZWavyLine::ClassId() const{
        return Hash("TPZWavyLine") ^ TPZGeoLinear::ClassId() << 1;
    }
    
    void TPZWavyLine::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> x0(lowercorner), x1(lowercorner), wavedir(3,0.);
        x1[0] = x0[0]+2.;
        wavedir[1] = 1.;
        int numwaves = 4;
        TPZManVector<int64_t> nodind(2);
        nodind[0] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodind[0]].Initialize(x0, gmesh);
        nodind[1] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodind[1]].Initialize(x1, gmesh);
        
        TPZGeoElRefPattern<TPZWavyLine> *gel = new TPZGeoElRefPattern<TPZWavyLine> (nodind,matid,gmesh);
        size[0] = 2.;
        size[1] = 2.;
        gel->Geom().SetData(wavedir, numwaves);
       
    }

}



template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZWavyLine>>;

template class TPZGeoElRefLess<pzgeom::TPZWavyLine>;

