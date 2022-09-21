//
//  TPZCylinder.cpp
//  pz
//
//  Created by Philippe Devloo on 04/05/18.
//

#include "TPZCylinderMap.h"
#include "pzvec_extras.h"
#include "tpzgeoelmapped.h"
#include "pzgeoquad.h"

namespace pzgeom {
    
    
    /// axis direction with the vertical axis
    template<class TGeo>
    void TPZCylinderMap<TGeo>::SetCylinderAxis(const TPZVec<REAL> &axis)
    {
        // build two normal vectors
        std::map<REAL,int> vals;
        for(int i=0; i<3; i++) vals[fabs(axis[i])] = i;
        TPZManVector<REAL,3> gen(3,0.), orto1(3,0.), orto2(3,0.);
        auto it = vals.rbegin();
        gen[it->second] = 1.;
        Cross(axis, gen, orto1);
        REAL normax = sqrt(Norm(axis));
        REAL normort1 = sqrt(Norm(orto1));
        for (int i=0; i<3; i++) {
            axis[i] /= normax;
            orto1[i] /= normort1;
        }
        Cross(orto1, axis, orto2);
        for (int i=0; i<3; i++) {
            fRotation(0,i) = orto1[i];
            fRotation(1,i) = orto2[i];
            fRotation(2,i) = axis[i];
        }
    }
    
    /// compute the corner coordinates of the corner nodes
    template<class TGeo>
    void TPZCylinderMap<TGeo>::ComputeCornerCoordinates(TPZGeoMesh &gmesh)
    {
        int nnodes = TGeo::NNodes;
        for (int in=0; in<nnodes; in++) {
            int64_t nodeindex = TGeo::fNodeIndexes[in];
            TPZManVector<REAL,3> co(3);
            gmesh.NodeVec()[nodeindex].GetCoordinates(co);
            TPZManVector<REAL,3> localco(3);
            for (int i=0; i<3; i++) {
                localco[i] = 0.;
                for (int j=0; j<3; j++) {
                    localco[i] += fRotation(i,j)*co[i]-fOrigin[i];
                }
            }
            REAL radius = sqrt(localco[0]*localco[0]+localco[1]*localco[1]);
            if(fabs(radius-fRadius) > 1.e-6) DebugStop();
            REAL theta = atan2(localco[1],localco[0]);
            REAL z = localco[2];
            fCornerCo(0,in) = theta;
            fCornerCo(1,in) = z;
        }
    }
    
    
    template<class TGeo>
    void TPZCylinderMap<TGeo>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size){
        TPZManVector<int64_t, TGeo::NNodes> nodeind(TGeo::NNodes,-1);
        TPZManVector<REAL,3> x(TGeo::Dimension,0);

        auto CalcX = [&lowercorner](int i, TPZVec<REAL> &x){
            switch(i){
            case 0:
                x = {1,0,0};
                break;
            case 1:
                x = {0,1,0};
                break;
            case 2:
                x = {0,1,1};
                break;
            case 3://only for quads
                x = {1,0,1};
                break;
            default:
                PZError<<__PRETTY_FUNCTION__
                       <<"\n invalid number of nodes!\n";
                DebugStop();
            }
            for(int i = 0; i < 3; i++) {x[i] += lowercorner[i];}
        };

        
        for(auto in = 0; in < TGeo::NNodes; in++){
            CalcX(in,x);
            nodeind[in] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeind[in]].Initialize(x,gmesh);

        }
        auto *gel = new TPZGeoElRefPattern<TPZCylinderMap<TGeo>> (nodeind,matid,gmesh);

        constexpr REAL radius{1};
        gel->Geom().SetOrigin({0,0,0}, radius);
        gel->Geom().SetCylinderAxis({0,0,1});
        gel->Geom().ComputeCornerCoordinates(gmesh);
        
        if constexpr (std::is_same_v<TGeo,TPZGeoTriangle>){
            size[0] = size[1] = 1;
        }else if constexpr (std::is_same_v<TGeo,TPZGeoQuad>){
            size[0] = size[1] = 2;
        }else{
            PZError<<__PRETTY_FUNCTION__
                   <<"\ninvalid element type!Aborting...\n";
            DebugStop();
        }
    }


};

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"

#define IMPLEMENTCYLINDERMAP(TGEO)                                    \
    template class pzgeom::TPZCylinderMap<TGEO>;                      \
    template class TPZGeoElRefLess<pzgeom::TPZCylinderMap<TGEO> >;    \
    template class TPZGeoElRefPattern<pzgeom::TPZCylinderMap<TGEO> >;

IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoTriangle)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoQuad)

#undef IMPLEMENTCYLINDERMAP
