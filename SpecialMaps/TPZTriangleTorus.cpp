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


    int TPZTriangleTorus::ClassId() const{
        return Hash("TPZTriangleTorus") ^ TPZGeoTriangle::ClassId() << 1;
    }

    void pzgeom::TPZTriangleTorus::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        REAL R = 1., r = 0.8;
        size[0] = 2.;
        size[1] = 2.;
        size[2] = 1.;
        TPZManVector<REAL,3> origin(lowercorner);
        origin[0] += 1.;
        
        TPZFNMatrix<12,REAL> phitheta(2,3,0.);
        phitheta(0,0) = 0.;
        phitheta(1,0) = 0;
        
        phitheta(0,1) = M_PI/3;
        phitheta(1,1) = 0;
        
        phitheta(0,2) = 0;
        phitheta(1,2) = M_PI;
        
        TPZManVector<int64_t,4> indexes(3);
        constexpr REAL coords[3][2] = {
            {0,0},{1,0},{0,1}
        };

        TPZTriangleTorus torus;
        torus.SetOrigin(origin);
        torus.SetDataRadius(R, r);
        torus.SetDataPhiTheta(phitheta);
        for (int i=0; i<3; i++) {
            indexes[i] = gmesh.NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3);
            TPZManVector<REAL, 2>loc(2);
            loc[0] = coords[i][0];
            loc[1] = coords[i][1];
            torus.X(phitheta, loc, xco);
            gmesh.NodeVec()[indexes[i]].Initialize(xco, gmesh);
        }
        TPZGeoElRefPattern<pzgeom::TPZTriangleTorus> *gel = new TPZGeoElRefPattern<pzgeom::TPZTriangleTorus>(indexes,matid,gmesh);
        gel->Geom().SetOrigin(origin);
        gel->Geom().SetDataRadius(R, r);
        gel->Geom().SetDataPhiTheta(phitheta);
    }

}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZTriangleTorus>>;


template class TPZGeoElRefLess<pzgeom::TPZTriangleTorus>;

