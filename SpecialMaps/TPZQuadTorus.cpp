//
//  TPZQuadTorus.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZQuadTorus.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef PZ_LOG
static PZLogger logger("pz.geom.pzgeoquad0");
#endif

namespace pzgeom {

// TPZGeoEl *TPZQuadTorus::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
// {
    
//         int ns = orig->NSideNodes(side);
//         TPZManVector<int64_t> nodeindices(ns);
//         int in;
//         for(in=0; in<ns; in++)
//         {
//             nodeindices[in] = orig->SideNodeIndex(side,in);
//         }
//         int64_t index;
        
//         TPZGeoMesh *mesh = orig->Mesh();
//         MElementType type = orig->Type(side);
        
//         TPZGeoEl *newel = mesh->CreateGeoBlendElement(type, nodeindices, bc, index);
//         TPZGeoElSide me(orig,side);
//         TPZGeoElSide newelside(newel,newel->NSides()-1);
        
//         newelside.InsertConnectivity(me);
//         newel->Initialize();
        
//         return newel;
// }


    // /**
    //  * Creates a geometric element according to the type of the father element
    //  */
    // /** @brief Creates a geometric element according to the type of the father element */
    // TPZGeoEl *TPZQuadTorus::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
    //                                   TPZVec<int64_t>& nodeindexes,
    //                                   int matid,
    //                                   int64_t& index)

    // {
    //     return ::CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    // }

    

    
//    void TPZQuadTorus::Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//    {
//        
//        TPZFNMatrix<9,REAL> GradPhi(3,3,0.);
//        DebugStop();
//        //TPZGeoQuad::Jacobian(fPhiTheta, param, jacobian, axes, detjac, jacinv);
//        TPZFNMatrix<6> axest(3,2);
//        axes.Transpose(&axest);
//        axest.Multiply(jacobian, GradPhi);
//        TPZFNMatrix<6,REAL> DxDphi(3,3,0.);
//        TPZManVector<REAL,3> ft(3,0.);
//        TPZGeoQuad::X(fPhiTheta,param,ft);
//        DxDphi(0,0) = -cos(ft[1]) * sin(ft[0]);
//        DxDphi(0,1) = -(3. + cos(ft[0])) * sin(ft[1]);
//        DxDphi(1,0) = -sin(ft[1]) * sin(ft[0]);
//        DxDphi(1,1) = cos(ft[1]) * (3. + cos(ft[0]));
//        DxDphi(2,0) = cos(ft[0]);
//        DxDphi(2,1) = 0.;
//        TPZFMatrix<REAL> VecMatrix;
//        DxDphi.Multiply(GradPhi, VecMatrix);
//
//        TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
//        
//        int spacedim = fPhiTheta.Rows();
//        
//        for (int j=0; j<spacedim; j++) {
//            minx[j] = fPhiTheta.GetVal(j,0);
//            maxx[j] = fPhiTheta.GetVal(j,0);
//        }
//
//        for(int i = 0; i < 4; i++) {
//            for(int j = 0; j < spacedim; j++) {
//                minx[j] = minx[j] < fPhiTheta.GetVal(j,i) ? minx[j]:fPhiTheta.GetVal(j,i);
//                maxx[j] = maxx[j] > fPhiTheta.GetVal(j,i) ? maxx[j]:fPhiTheta.GetVal(j,i);
//            }
//        }
//        REAL delx = 0.;
//        for (int j=0; j<spacedim; j++) {
//            delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
//        }
//        VecMatrix *= 1./delx;
//        
//        VecMatrix.GramSchmidt(axest,jacobian);
//        axest.Transpose(&axes);
//        detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
//        
//        if(IsZero(detjac))
//        {
//#ifdef PZDEBUG
//            std::stringstream sout;
//            sout << "Singular Jacobian " << detjac;
//            LOGPZ_ERROR(logger, sout.str())
//#endif
//            detjac = ZeroTolerance();
//        }
//        
//        jacinv(0,0) =  jacobian(1,1)/detjac;
//        jacinv(1,1) =  jacobian(0,0)/detjac;
//        jacinv(0,1) = -jacobian(0,1)/detjac;
//        jacinv(1,0) = -jacobian(1,0)/detjac;
//        
//        jacobian *= delx;
//        jacinv *= 1./delx;
//        detjac *= (delx*delx);
//        
//    }
        
    int TPZQuadTorus::ClassId() const {
        return Hash("TPZQuadTorus") ^ TPZGeoQuad::ClassId() << 1;
    }

    void pzgeom::TPZQuadTorus::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        REAL R = 1., r = 0.8;
        size[0] = 2.;
        size[1] = 2.;
        size[2] = 1.;
        TPZManVector<REAL,3> origin(lowercorner);
        origin[0] += 1.;
        TPZQuadTorus torus;
        torus.SetOrigin(origin);
        torus.SetDataRadius(R, r);
        TPZFNMatrix<12,REAL> phitheta(2,4,0.);
        phitheta(0,0) = 0.;
        phitheta(0,1) = M_PI/5.;
        phitheta(0,2) = M_PI;
        phitheta(0,3) = M_PI+M_PI/3.;
        
        phitheta(1,0) = M_PI;
        phitheta(1,1) = 2*M_PI;
        phitheta(1,2) = 2*M_PI-M_PI/5.;
        phitheta(1,3) = M_PI;
        torus.SetDataPhiTheta(phitheta);
        TPZManVector<int64_t,4> indexes(4);
        REAL coords[4][2] = {
            {-1,-1},{1,-1},{1,1},{-1,1}
        };
        for (int i=0; i<4; i++) {
            indexes[i] = gmesh.NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3), loc(2);
            loc[0] = coords[i][0];
            loc[1] = coords[i][1];
            
            torus.X(phitheta, loc, xco);
            gmesh.NodeVec()[indexes[i]].Initialize(xco, gmesh);
        }
        TPZGeoElRefPattern<pzgeom::TPZQuadTorus> *gel = new TPZGeoElRefPattern<pzgeom::TPZQuadTorus>(indexes,matid,gmesh);
        gel->Geom().SetOrigin(origin);
        gel->Geom().SetDataRadius(R, r);
        gel->Geom().SetDataPhiTheta(phitheta);
    }
    

}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadTorus>>;

template class TPZGeoElRefLess<pzgeom::TPZQuadTorus>;
