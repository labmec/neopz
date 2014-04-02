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

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoquad0"));
#endif

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


    /**
     * Creates a geometric element according to the type of the father element
     */
    /** @brief Creates a geometric element according to the type of the father element */
    TPZGeoEl *TPZQuadTorus::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                      TPZVec<long>& nodeindexes,
                                      int matid,
                                      long& index)

    {
        return ::CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
    }

	void TPZQuadTorus::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
	{
		
		TPZGeoQuad::X(this->fPhiTheta,loc,result);
		TPZVec <REAL> toro(3,0.0);
		
		toro[0] = (fR + fr*cos(result[0]))*cos(result[1]);
		toro[1] = (fR + fr*cos(result[0]))*sin(result[1]);		
		toro[2] = fr*sin(result[0]);		
		result=toro;

	}
    
	void TPZQuadTorus::Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
	{
		
		TPZFNMatrix<9,REAL> GradPhi(3,3,0.);
		TPZGeoQuad::Jacobian(fPhiTheta, param, jacobian, axes, detjac, jacinv);
		TPZFNMatrix<6> axest(3,2);
		axes.Transpose(&axest);
		axest.Multiply(jacobian, GradPhi);
		TPZFNMatrix<6,REAL> DxDphi(3,3,0.);
		TPZManVector<REAL,3> ft(3,0.);
		TPZGeoQuad::X(fPhiTheta,param,ft);
		DxDphi(0,0) = -cos(ft[1]) * sin(ft[0]);
		DxDphi(0,1) = -(3. + cos(ft[0])) * sin(ft[1]);
		DxDphi(1,0) = -sin(ft[1]) * sin(ft[0]);
		DxDphi(1,1) = cos(ft[1]) * (3. + cos(ft[0]));
		DxDphi(2,0) = cos(ft[0]);
		DxDphi(2,1) = 0.;
		TPZFMatrix<REAL> VecMatrix;
		DxDphi.Multiply(GradPhi, VecMatrix);

		TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
		
		int spacedim = fPhiTheta.Rows();
        
        for (int j=0; j<spacedim; j++) {
            minx[j] = fPhiTheta.GetVal(j,0);
            maxx[j] = fPhiTheta.GetVal(j,0);
        }

		for(int i = 0; i < 4; i++) {
			for(int j = 0; j < spacedim; j++) {
                minx[j] = minx[j] < fPhiTheta.GetVal(j,i) ? minx[j]:fPhiTheta.GetVal(j,i);
                maxx[j] = maxx[j] > fPhiTheta.GetVal(j,i) ? maxx[j]:fPhiTheta.GetVal(j,i);
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

}




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

