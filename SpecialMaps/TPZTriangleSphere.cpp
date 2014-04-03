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
#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeotriangle"));
#endif

namespace pzgeom {

TPZGeoEl *TPZTriangleSphere::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
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
    TPZGeoEl *TPZTriangleSphere::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                      TPZVec<long>& nodeindexes,
                                      int matid,
                                      long& index)

    {
        return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    }

    
	
	/* @brief Computes the jacobian of the map between the master element and deformed element */
	void TPZTriangleSphere::Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
	{
		TPZVec<REAL> XTriangle(3,0.0);
		TPZFNMatrix<9,REAL> GradOneoverNorm(3,1,0.0);
		TPZFNMatrix<9,REAL> TensorXtGradX(3,2,0.0);	
		TPZFNMatrix<9,REAL> GradXtScaled(3,3,0.0);	
		
		TPZFNMatrix<3*NNodes> coord(3,NNodes);
		CornerCoordinates(gel, coord);
		TPZGeoTriangle::X(coord,param,XTriangle);
		TPZFMatrix<REAL> XtminusXc(3,1,0.0);
		TPZVec<REAL> Xc= fXc;
		
		for (int i = 0; i < XTriangle.size(); i++) { XtminusXc(i,0)= XTriangle[i] - Xc[i]; }
		REAL NormValue = Norm(XtminusXc);		
		
		TPZGeoTriangle::Jacobian(gel, param, jacobian , axes, detjac, jacinv);	
		TPZFNMatrix<6> axest(3,2);
		axes.Transpose(&axest);
		
		TPZFNMatrix<9,REAL> GradXt;
		axest.Multiply(jacobian, GradXt);					
		
		GradOneoverNorm(0,0) = XtminusXc(0,0)*GradXt(0,0)+XtminusXc(1,0)*GradXt(1,0)+XtminusXc(2,0)*GradXt(2,0);
		GradOneoverNorm(1,0) = XtminusXc(0,0)*GradXt(0,1)+XtminusXc(1,0)*GradXt(1,1)+XtminusXc(2,0)*GradXt(2,1);
		GradOneoverNorm = -(fR/(NormValue*NormValue*NormValue))*GradOneoverNorm;
		
		TensorXtGradX(0,0)= XtminusXc(0,0)*GradOneoverNorm(0,0);
		TensorXtGradX(0,1)= XtminusXc(0,0)*GradOneoverNorm(1,0);
		TensorXtGradX(1,0)= XtminusXc(1,0)*GradOneoverNorm(0,0);
		TensorXtGradX(1,1)= XtminusXc(1,0)*GradOneoverNorm(1,0);			
		TensorXtGradX(2,0)= XtminusXc(2,0)*GradOneoverNorm(0,0);
		TensorXtGradX(2,1)= XtminusXc(2,0)*GradOneoverNorm(1,0);
		
		TPZFNMatrix<9,REAL> VecMatrix;
		VecMatrix=((fR/(NormValue))*GradXt)+TensorXtGradX;		
			
			TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
			
			int spacedim = coord.Rows();
			
			for (int j=0; j<spacedim; j++) { minx[j] = coord.GetVal(j,0); maxx[j] = coord.GetVal(j,0);}
			
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < spacedim; j++) {
					minx[j] = minx[j] < coord.GetVal(j,i) ? minx[j]:coord.GetVal(j,i);
					maxx[j] = maxx[j] > coord.GetVal(j,i) ? maxx[j]:coord.GetVal(j,i);
				}
			}
			REAL delx = 0.;
			for (int j=0; j<spacedim; j++) {
				delx = delx > (maxx[j]-minx[j]) ? delx : (maxx[j]-minx[j]);
			}
		
			VecMatrix *= 1./delx;	
			VecMatrix.GramSchmidt(axest,jacobian);
			axest.Transpose(&axes);
			detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
			
			if(IsZero(detjac))
			{
#ifdef DEBUG
				std::stringstream sout;
				sout << "Singular Jacobian " << detjac;
				LOGPZ_ERROR(logger, sout.str())
#endif
				detjac = ZeroTolerance();
			}

			jacinv(0,0) =  jacobian(1,1)/detjac;
			jacinv(1,1) =  jacobian(0,0)/detjac;
			jacinv(0,1) = -jacobian(0,1)/detjac;
			jacinv(1,0) = -jacobian(1,0)/detjac;			
			jacobian *= delx;
			jacinv *= 1./delx;
			detjac *= (delx*delx);			
			
		
	}
	
	void TPZTriangleSphere::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
	{
		
		TPZGeoTriangle::X(nodes,loc,result);
		TPZFMatrix<REAL> XtminusXc(3,1,0.0);
		TPZVec<REAL> Xc= fXc;
		
		XtminusXc(0,0)= result[0] - Xc[0];
		XtminusXc(1,0)= result[1] - Xc[1];
		XtminusXc(2,0)= result[2] - Xc[2];
		REAL NormValue = Norm(XtminusXc);
		
		TPZVec<REAL> Xsphere(3,0.0);
		Xsphere[0] = (fR/NormValue)*(XtminusXc(0,0)+Xc[0]);
		Xsphere[1] = (fR/NormValue)*(XtminusXc(1,0)+Xc[1]);
		Xsphere[2] = (fR/NormValue)*(XtminusXc(2,0)+Xc[2]);		
		result=Xsphere;
		
	}	
	

}

/**
 * @ingroup geometry
 * @brief Id for three dimensional arc element
 */

template<>
int TPZGeoElRefPattern<pzgeom::TPZTriangleSphere>::ClassId() const {
	return TPZGEOELEMENTTRIANGLESPHEREID;
}

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZTriangleSphere>, TPZGEOELEMENTTRIANGLESPHEREID>;


template class TPZGeoElRefLess<pzgeom::TPZTriangleSphere>;

