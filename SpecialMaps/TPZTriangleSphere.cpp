//
//  TPZTriangleSphere.cpp
//  PZ
//
//  Created by Philippe Devloo on 3/21/14.
//
//

#include "TPZTriangleSphere.h"
#include "tpzgeomid.h"
#include "tpzgeoelmapped.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeotriangle");
#endif

namespace pzgeom {

    
    template<class GeomTriang>
    void TPZTriangleSphere<GeomTriang>::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &sz)
    {
        TPZManVector<REAL,3> center(lowercorner);
        REAL radius = 1.;
        center[0] += radius/2.;
        center[1] += radius/2.;
        center[2] += radius/2.;
        for (int i=0; i<3; i++) {
            sz[i] = 2.*radius;
        }
        REAL coords[3][3] = {
            {-1,-1,-0.1},
            { 1,-1,-0.1},
            { 1, 1,-0.1}
        };
        for (int i=0; i<3; i++) {
            REAL norm = sqrt(coords[i][0]*coords[i][0]+coords[i][1]*coords[i][1]+coords[i][2]*coords[i][2]);
            for(int j=0; j<3; j++) coords[i][j] *= radius/norm;
        }
        TPZManVector<int64_t,3> indices(3);
        for (int i=0; i<3; i++) {
            indices[i] = gmesh.NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3);
            for (int j=0; j<3; j++) {
                xco[j] = coords[i][j]+center[j];
            }
            gmesh.NodeVec()[indices[i]].Initialize(xco, gmesh);
        }
        TPZGeoElRefPattern<TPZTriangleSphere<> > *gel = new TPZGeoElRefPattern<TPZTriangleSphere<> >(indices,matid,gmesh);
        gel->Geom().SetData(radius, center);
    }

    // /**
    //  * Creates a geometric element according to the type of the father element
    //  */
    // /** @brief Creates a geometric element according to the type of the father element */
    // template< class GeomTriang>
    // TPZGeoEl *TPZTriangleSphere<GeomTriang>::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
    //                                   TPZVec<int64_t>& nodeindexes,
    //                                   int matid,
    //                                   int64_t& index)

    // {
    //     return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    // }

    template< class GeomTriang>
    bool TPZTriangleSphere<GeomTriang>::IsGeoBlendEl() const
    {
        return false;
    }
	
    /** @brief declare geometry as blended element */
    template<>
    bool TPZTriangleSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> >::IsGeoBlendEl() const
    {
        return true;
    }

	
	/* @brief Computes the jacobian of the map between the master element and deformed element */
//	void TPZTriangleSphere::Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
//	{
//		TPZVec<REAL> XTriangle(3,0.0);
//		TPZFNMatrix<9,REAL> GradOneoverNorm(3,1,0.0);
//		TPZFNMatrix<9,REAL> TensorXtGradX(3,2,0.0);	
//		TPZFNMatrix<9,REAL> GradXtScaled(3,3,0.0);	
//		
//		TPZFNMatrix<3*NNodes> coord(3,NNodes);
//		CornerCoordinates(gel, coord);
//		TPZGeoTriangle::X(coord,param,XTriangle);
//		TPZFMatrix<REAL> XtminusXc(3,1,0.0);
//		TPZVec<REAL> Xc= fXc;
//		
//		for (int i = 0; i < XTriangle.size(); i++) { XtminusXc(i,0)= XTriangle[i] - Xc[i]; }
//		REAL NormValue = Norm(XtminusXc);		
//		
//		TPZGeoTriangle::Jacobian(gel, param, jacobian , axes, detjac, jacinv);	
//		TPZFNMatrix<6> axest(3,2);
//		axes.Transpose(&axest);
//		
//		TPZFNMatrix<9,REAL> GradXt;
//		axest.Multiply(jacobian, GradXt);					
//		
//		GradOneoverNorm(0,0) = XtminusXc(0,0)*GradXt(0,0)+XtminusXc(1,0)*GradXt(1,0)+XtminusXc(2,0)*GradXt(2,0);
//		GradOneoverNorm(1,0) = XtminusXc(0,0)*GradXt(0,1)+XtminusXc(1,0)*GradXt(1,1)+XtminusXc(2,0)*GradXt(2,1);
//		GradOneoverNorm = -(fR/(NormValue*NormValue*NormValue))*GradOneoverNorm;
//		
//		TensorXtGradX(0,0)= XtminusXc(0,0)*GradOneoverNorm(0,0);
//		TensorXtGradX(0,1)= XtminusXc(0,0)*GradOneoverNorm(1,0);
//		TensorXtGradX(1,0)= XtminusXc(1,0)*GradOneoverNorm(0,0);
//		TensorXtGradX(1,1)= XtminusXc(1,0)*GradOneoverNorm(1,0);			
//		TensorXtGradX(2,0)= XtminusXc(2,0)*GradOneoverNorm(0,0);
//		TensorXtGradX(2,1)= XtminusXc(2,0)*GradOneoverNorm(1,0);
//		
//		TPZFNMatrix<9,REAL> VecMatrix;
//		VecMatrix=((fR/(NormValue))*GradXt)+TensorXtGradX;		
//			
//			TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
//			
//			int spacedim = coord.Rows();
//			
//			for (int j=0; j<spacedim; j++) { minx[j] = coord.GetVal(j,0); maxx[j] = coord.GetVal(j,0);}
//			
//			for(int i = 0; i < 3; i++) {
//				for(int j = 0; j < spacedim; j++) {
//					minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
//					maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
//				}
//			}
//			REAL delx = 0.;
//			for (int j=0; j<spacedim; j++) {
//				delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
//			}
//		
//			VecMatrix *= 1./delx;	
//			VecMatrix.GramSchmidt(axest,jacobian);
//			axest.Transpose(&axes);
//			detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
//			
//			if(IsZero(detjac))
//			{
//#ifdef PZDEBUG
//				std::stringstream sout;
//				sout << "Singular Jacobian " << detjac;
//				LOGPZ_ERROR(logger, sout.str())
//#endif
//				detjac = ZeroTolerance();
//			}
//
//			jacinv(0,0) =  jacobian(1,1)/detjac;
//			jacinv(1,1) =  jacobian(0,0)/detjac;
//			jacinv(0,1) = -jacobian(0,1)/detjac;
//			jacinv(1,0) = -jacobian(1,0)/detjac;			
//			jacobian *= delx;
//			jacinv *= 1./delx;
//			detjac *= (delx*delx);			
//			
//		
//	}
	
//	void TPZTriangleSphere::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
//	{
//		
//		TPZGeoTriangle::X(nodes,loc,result);
//		TPZFMatrix<REAL> XtminusXc(3,1,0.0);
//		TPZVec<REAL> Xc= fXc;
//		
//		XtminusXc(0,0)= result[0] - Xc[0];
//		XtminusXc(1,0)= result[1] - Xc[1];
//		XtminusXc(2,0)= result[2] - Xc[2];
//		REAL NormValue = Norm(XtminusXc);
//		
//		TPZVec<REAL> Xsphere(3,0.0);
//		Xsphere[0] = (fR/NormValue)*(XtminusXc(0,0)+Xc[0]);
//		Xsphere[1] = (fR/NormValue)*(XtminusXc(1,0)+Xc[1]);
//		Xsphere[2] = (fR/NormValue)*(XtminusXc(2,0)+Xc[2]);		
//		result=Xsphere;
//		
//	}	
	
}

template class pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle >;

template class pzgeom::TPZTriangleSphere< pzgeom::TPZGeoBlend< pzgeom::TPZGeoTriangle > >;

template class TPZGeoElRefLess<pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > >;

template class TPZRestoreClass< TPZGeoElRefPattern< pzgeom::TPZTriangleSphere<pzgeom::TPZGeoTriangle > >>;

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZTriangleSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> > >>;

template class TPZGeoElRefLess<pzgeom::TPZTriangleSphere<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle> > >;



