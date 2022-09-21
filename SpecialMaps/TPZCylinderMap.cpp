//
//  TPZCylinder.cpp
//  pz
//
//  Created by Philippe Devloo on 04/05/18.
//

#include "TPZCylinderMap.h"
#include "pzfmatrix.h"
#include "pzvec_extras.h"
#include "tpzgeoelmapped.h"
#include "pzgeoquad.h"

namespace pzgeom {
    
    constexpr REAL tol = std::numeric_limits<REAL>::epsilon()*1000;
    /// axis direction with the vertical axis
    template<class TGeo>
    void TPZCylinderMap<TGeo>::SetCylinderAxis(const TPZVec<REAL> &axis)
    {
        //master cylinder has axis in this direction
        TPZManVector<REAL,3> orig_axis = {0,0,1};

        //let us normalize axis
        REAL normaxis = 0;
        for(int ix = 0; ix < axis.size(); ix++){normaxis += axis[ix]*axis[ix];}
        normaxis = sqrt(normaxis);
        /*
          We know three things about this rotation matrix:
          Let us assume that x is the original axis and y is the new axis.
          1. Rx = y
          2. R^t y = x
          3. R(x \times y) = x \times y
         */
        const auto &x = orig_axis;
        const auto &y = axis;
        TPZManVector<REAL,3> z(3,0.);
        //Let us compute the cross product between x and y
        Cross(x,y,z);

        if(fabs(z[0]) < tol && fabs(z[1]) < tol && fabs(z[2]) < tol){
            //x and y are aligned
            //let us compute the inner product
            REAL inner{0};
            for(int ix = 0; ix < 3; ix++){inner += x[ix] * y[ix];}
            const REAL sign = inner > 0 ? 1 : - 1;
            for(int i = 0; i< 3; i++){
                fRotation(i,i) = sign;
            }
        }else{
            //now we build a 9x9 matrix to find the rotation matrix
            TPZFNMatrix<81,REAL> M(9,9,0.), b(9,1,0.);
            for(int c = 0; c < 3; c++){
                for(int i = 0; i < 3; i++){
                    for(int j =0 ; j < 3; j++){
                        M(3*c+i,3*c+j) = x[j];
                        M(3*c+i,3*j) = y[j];
                        M(3*c+2+i,3*c+j) = z[j];
                    }
                }
                b(3*c+0,0) = y[c];
                b(3*c+1,0) = x[c];
                b(3*c+2,0) = z[c];
            }

            M.Solve_Cholesky(&b);
            for(int i = 0; i < 9; i++){
                const auto r = i/3;
                const auto c = i%3;
                fRotation(r,c) = b(i,0);
            }
        }
        fRotation.Print(std::cout);
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
                    localco[i] += fRotation(j,i)*(co[i]-fOrigin[i]);
                }
            }
            REAL radius = sqrt(localco[0]*localco[0]+localco[1]*localco[1]);
            if(fabs(radius-fRadius) > tol) {
                PZError<<__PRETTY_FUNCTION__
                       <<"\nError:"
                       <<"coordinates "<< co
                       <<"\nlocal coordinates: "<<localco
                       <<"\nradius: "<<fRadius
                       <<"\ncomputed radius: "<<radius
                       <<"Aborting..."<<std::endl;
                DebugStop();
            }
            REAL theta = atan2(localco[1],localco[0]);
            REAL z = localco[2];
            fCornerCo(0,in) = theta;
            fCornerCo(1,in) = z;
        }
    }
    
    
    template<class TGeo>
    void TPZCylinderMap<TGeo>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size){
        TPZManVector<int64_t, TGeo::NNodes> nodeind(TGeo::NNodes,-1);
        TPZManVector<REAL,3> x(3,0);

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
        gel->Geom().SetOrigin(lowercorner, radius);
        gel->Geom().SetCylinderAxis({0,0,1});
        gel->Geom().ComputeCornerCoordinates(gmesh);
        
        size[0] = size[1] = 1;
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
