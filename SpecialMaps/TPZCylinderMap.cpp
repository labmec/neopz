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
    template<class TGeo>
    void TPZCylinderMap<TGeo>::SetCylinderAxis(const TPZVec<REAL> &orig_axis)
    {

        TPZManVector<REAL,3> reference_axis = {0,0,1};

        //let us normalize axis
        REAL normaxis = 0;
        for(const auto &ax : orig_axis){normaxis += ax*ax;}
        normaxis = sqrt(normaxis);
        const TPZManVector<REAL,3> axis = {orig_axis[0]/normaxis,
                                           orig_axis[1]/normaxis,
                                           orig_axis[2]/normaxis};
        const auto &x = reference_axis;
        const auto &y = axis;
        TPZManVector<REAL,3> orth1(3,0.);
        Cross(x,y,orth1);

        //if they are aligned, then it is just a matter of sign
        if(fabs(orth1[0]) < tol && fabs(orth1[1]) < tol && fabs(orth1[2]) < tol){
            //x and y are aligned
            //let us compute the inner product
            REAL inner{0};
            for(int ix = 0; ix < 3; ix++){inner += x[ix] * y[ix];}
            const REAL sign = inner > 0 ? 1 : - 1;
            for(int i = 0; i< 3; i++){
                fRotation(i,i) = sign;
            }
        }else{
            /**if they are not aligned, orth1 is already an orth vector,
               we need to normalise it*/
            {
                REAL normorth1{0};
                for(auto &xx : orth1) {normorth1 += xx*xx;}
                normorth1 = sqrt(normorth1);
                for(auto &xx : orth1) {xx /= normorth1;}
            }

            TPZManVector<REAL,3> orth2(3,0.);
            Cross(y,orth1,orth2);
            for(int i = 0; i < 3; i++){
                fRotation(i,0) = orth1[i];
                fRotation(i,1) = orth2[i];
                fRotation(i,2) = y[i];
        
            }
        }
    }
    
    template<class TGeo>
    void TPZCylinderMap<TGeo>::SetRotationMatrix(const TPZFMatrix<REAL> &orig_axis)
    {
        /**master cylinder has axis in the z-direction.
         therefore, a rotation matrix to convert from reference axis
        to cylinder's axis will have the cylinder's axis as the 3 column
        (for [0,0,1] is transformed to the cylinder's axis).
        we now need to find two orthonormal vectors to compose our matrix*/
        fRotation = orig_axis;
#ifdef PZDEBUG
        for(int i = 0; i<3; i++) {
            for(int j= 0; j<3; j++) {
                REAL inner = 0.;
                for(int k = 0; k<3; k++) {
                    inner += orig_axis(i,k)*orig_axis(j,k);
                }
                if(i==j && fabs(inner-1.) > 1.e-8) {
                    DebugStop();
                } else if(i!=j && fabs(inner) > 1.e-8) {
                    DebugStop();
                }
            }
        }
#endif
    }
    
    /// compute the corner coordinates of the corner nodes
    template<class TGeo>
    void TPZCylinderMap<TGeo>::ComputeCornerCoordinates(TPZGeoMesh &gmesh)
    {
        constexpr int nnodes = TGeo::NNodes;
        //atan2 returns in the range -pi,pi
        REAL mintheta{2*M_PI}, maxtheta{-2*M_PI};
        for (int in=0; in<nnodes; in++) {
            const int64_t nodeindex = TGeo::fNodeIndexes[in];
            TPZManVector<REAL,3> co(3);
            gmesh.NodeVec()[nodeindex].GetCoordinates(co);
            TPZManVector<REAL,3> localco(3);
            for (int i=0; i<3; i++) {
                localco[i] = 0.;
                for (int j=0; j<3; j++) {
                    localco[i] += fRotation(j,i)*(co[j]-fOrigin[j]);
                }
            }
            const REAL radius = sqrt(localco[0]*localco[0]+localco[1]*localco[1]);
            const REAL theta = atan2(localco[1],localco[0]);
            const REAL z = localco[2];
            fCylindricalCo(0,in) = radius;
            fCylindricalCo(1,in) = theta;
            fCylindricalCo(2,in) = z;
            if(theta > maxtheta) {maxtheta = theta;}
            if(theta < mintheta) {mintheta = theta;}
        }


        /*
          now we need choose a proper range for theta. all angles
          should be no further spaced than pi.
          
          examples of possible mistakes:
          range: -pi, pi (the one from atan2)
          angles: -3pi/4, 3pi/4
          computed average: 0
          desired average: -pi or pi

          range: 0, 2pi
          angles: pi/4, 7pi/4
          computed average: pi
          desired average: 0 or 2pi
         */
        if(maxtheta - mintheta > M_PI){
            mintheta = 2*M_PI;
            maxtheta = -mintheta;
            for(auto in = 0 ; in < nnodes ; in++){
                const REAL theta = std::fmod((fCylindricalCo(1,in)+2*M_PI),2*M_PI); 
                fCylindricalCo(1,in) = theta;
                if(theta > maxtheta) {maxtheta = theta;}
                if(theta < mintheta) {mintheta = theta;}
            }
            if(maxtheta - mintheta > M_PI){
                PZError<<__PRETTY_FUNCTION__
                       <<"\nUnable to find suitable range for converting "
                       <<"this element's corner nodes to cylindrical coordinates.\n"
                       <<"Computed coordinates were:\n";
                TPZManVector<REAL,3> co(3);
                for(auto in = 0 ; in < nnodes ; in++){
                    const int64_t nodeindex = TGeo::fNodeIndexes[in];
                    gmesh.NodeVec()[nodeindex].GetCoordinates(co);
                    PZError<<"\tnode "<<in
                           <<"\n\t\tcartesian:";
                    for(int ix = 0; ix < 3; ix++){PZError<<' '<<co[ix];}
                    PZError<<"\n\t\tcylindrical:";
                    for(int ix = 0; ix < 3; ix++){PZError<<' '<<fCylindricalCo(ix,in);}
                }
                PZError <<"\nAborting...";
                DebugStop();
            }
        }
    }

    template<class TGeo>
    void TPZCylinderMap<TGeo>::Read(TPZStream& buf, void* context)
    {
        TGeo::Read(buf,0);
        fCylindricalCo.Read(buf,0);
        buf.Read(fOrigin);
        fRotation.Read(buf, 0);
    }
        
    template<class TGeo>
    void TPZCylinderMap<TGeo>::Write(TPZStream &buf, int withclassid) const
    {
        TGeo::Write(buf, withclassid);
        fCylindricalCo.Write(buf,0);
        buf.Write(fOrigin);
        fRotation.Write(buf, withclassid);
    }
    
    
    template<class TGeo>
    void TPZCylinderMap<TGeo>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &elsize){
        TPZManVector<int64_t, TGeo::NNodes> nodeind(TGeo::NNodes,-1);
        TPZFNMatrix<24,REAL> x = {
            {0.5,0,0},
            {1,0,0},
            {0,1,0},
            {0,0.5,0},
            {0.5,0,0.5},
            {1,0,0.5},
            {0,1,0.5},
            {0,0.5,0.5}
        };
        x.Transpose();
        for(int i=0; i<8; i++) for(int jc = 0; jc<3; jc++) {
            x(jc,i) += lowercorner[jc];
        }
        
        // inserting the nodes in the mesh
        const int64_t firstnode = gmesh.NodeVec().NElements();
        gmesh.NodeVec().Resize(firstnode + x.Cols());
        for(int i = 0 ; i < x.Cols() ; i++) {
            TPZManVector<REAL,3> coor(3,0.);
            for(int in = 0 ; in < 3 ; in++){
                coor[in] = x(in,i);
            }
            gmesh.NodeVec()[i+firstnode].Initialize(coor, gmesh);
        }

        TPZManVector<int64_t,8> nodind(TGeo::NCornerNodes);
        switch(TGeo::Type()) {
            case EOned:
                nodeind = {0,1};
                break;
            case ETriangle:
                nodeind = {0,1,2};
                break;
            case EQuadrilateral:
                nodeind = {0,1,2,3};
                break;
            case EPrisma:
                nodeind = {0,1,2,4,5,6};
                break;
            case EPiramide:
                nodeind = {0,1,2,3,6};
                break;
            case ETetraedro:
                nodeind = {0,1,2,4};
                break;
            case ECube:
                nodeind = {0,1,2,3,4,5,6,7};
                break;
            default:
                DebugStop();
        }
        for(auto& it : nodeind) it += firstnode;                
        
        auto *gel = new TPZGeoElRefPattern<TPZCylinderMap<TGeo>> (nodeind,matid,gmesh);

        constexpr REAL radius{1};
        gel->Geom().SetOrigin(lowercorner);
        TPZFNMatrix<9,REAL> axis(3,3);
        axis.Identity();
        gel->Geom().SetRotationMatrix(axis);
        gel->Geom().ComputeCornerCoordinates(gmesh);
        
        elsize[0] = elsize[1] = elsize[2] = 1;
    }


};

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"

#define IMPLEMENTCYLINDERMAP(TGEO)                                    \
    template class pzgeom::TPZCylinderMap<TGEO>;                      \
    template class TPZGeoElRefLess<pzgeom::TPZCylinderMap<TGEO> >;    \
    template class TPZGeoElRefPattern<pzgeom::TPZCylinderMap<TGEO> >;

IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoLinear)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoTriangle)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoQuad)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoTetrahedra)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoCube)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoPrism)
IMPLEMENTCYLINDERMAP(pzgeom::TPZGeoPyramid)

#undef IMPLEMENTCYLINDERMAP
