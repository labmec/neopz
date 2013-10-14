/**
 * @file
 * @brief Contains the implementation of the TPZGeoTetrahedra methods. 
 */

#include "pzgeotetrahedra.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "tpzgeoelrefpattern.h"

using namespace pzshape;
using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeotetrahedra"));
#endif

namespace pzgeom {
	
	const double tol = pzgeom_TPZNodeRep_tol;
	
	void TPZGeoTetrahedra::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
		phi(1,0)  = pt[0];
		phi(2,0)  = pt[1];
		phi(3,0)  = pt[2];
		
		dphi(0,0) = -1.0;
		dphi(1,0) = -1.0;
		dphi(2,0) = -1.0;
		dphi(0,1) =  1.0;
		dphi(1,1) =  0.0;
		dphi(2,1) =  0.0;
		dphi(0,2) =  0.0;
		dphi(1,2) =  1.0;
		dphi(2,2) =  0.0;
		dphi(0,3) =  0.0;
		dphi(1,3) =  0.0;
		dphi(2,3) =  1.0;
	}
	
	void TPZGeoTetrahedra::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
        jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
		REAL spacephi[10];
		TPZFMatrix<REAL> phi(4,1,spacephi,10);
		REAL spacedphi[20];
		TPZFMatrix<REAL> dphi(3,4,spacedphi,20);
		Shape(param,phi,dphi);
		jacobian.Zero();
		
		int i,j;
		for(i=0;i<4;i++) {
			for(j=0;j<3;j++) {
				jacobian(j,0) += coord.GetVal(j,i)*dphi(0,i);
				jacobian(j,1) += coord.GetVal(j,i)*dphi(1,i);
				jacobian(j,2) += coord.GetVal(j,i)*dphi(2,i);
			}
		}
		
		detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0);//-a02 a11 a20
		detjac += jacobian(0,1)*jacobian(1,2)*jacobian(2,0);//+ a01 a12 a20
		detjac += jacobian(0,2)*jacobian(1,0)*jacobian(2,1);//+ a02 a10 a21
		detjac -= jacobian(0,0)*jacobian(1,2)*jacobian(2,1);//- a00 a12 a21
		detjac -= jacobian(0,1)*jacobian(1,0)*jacobian(2,2);//- a01 a10 a22
		detjac += jacobian(0,0)*jacobian(1,1)*jacobian(2,2);//+ a00 a11 a22
		
		if(IsZero(detjac))
		{
#ifdef DEBUG
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
			LOGPZ_ERROR(logger, sout.str())
#endif
			detjac = ZeroTolerance();
		}
		
		jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
		jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;//a02 a21 - a01 a22
		jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
		jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;//a12 a20 - a10 a22
		jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
		jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;//a02 a10 - a00 a12
		jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
		jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;//a01 a20 - a00 a21
		jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11
		
		axes.Zero();
		axes(0,0) = 1.;
		axes(1,1) = 1.;
		axes(2,2) = 1.;
	}
	
	void TPZGeoTetrahedra::X(const TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
		REAL spacephi[10],spacedphi[20];
		int i,j;
		TPZFMatrix<REAL> phi(4,1,spacephi,10);
		TPZFMatrix<REAL> dphi(3,4,spacedphi,20);
		Shape(loc,phi,dphi);
		for(j=0;j<3;j++) {
			result[j] = 0.0;
			for(i=0;i<4;i++) 
				result[j] += coord.GetVal(j,i)*phi(i,0);
		}
	}
	
	TPZGeoEl *TPZGeoTetrahedra::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		if(side<0 || side>14){
			cout << "TPZGeoTetrahedra::CreateBCCompEl with bad side = " << side << "not implemented\n";	
			return 0;
		}
		
		if(side==14) {
			cout << "TPZGeoTetrahedra::CreateBCCompEl with side = 14 not implemented\n";
			return 0;
		}
		if(side<4) {
			TPZManVector<long> nodeindexes(1);
			//		TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} else if (side > 3 && side < 10) {//side =4 a 9 : lados
			TPZManVector<long> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			//		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} else if (side > 9) {//side = 10 a 13 : faces
			TPZManVector<long> nodes(3);
			int in;
			for (in=0;in<3;in++){
				nodes[in] = orig->SideNodeIndex(side,in);
			}
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
			//		TPZGeoElT2d *gel = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
			for (in=0;in<6;in++){
				TPZGeoElSide(gel,in).SetConnectivity(TPZGeoElSide(orig,TPZShapeTetra::ContainedSideLocId(side,in)));
			}
			TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} 
		else PZError << "TPZGeoTetrahedra::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	void TPZGeoTetrahedra::FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
	{
		ChangedPoint.Resize(OriginalPoint.NElements(),0.);
		ChangedPoint = OriginalPoint;
		
		switch(side)
		{
			case 4:
			{
				if(OriginalPoint[1] == 1. || OriginalPoint[2] == 1.)
				{
					REAL den = 0.25 + OriginalPoint[1]*OriginalPoint[1] + OriginalPoint[2]*OriginalPoint[2];
					ChangedPoint[0] = (0.5*tol)/sqrt(den);
					ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1]*tol) / sqrt(den);
					ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2]*tol)/sqrt(den);
				}
				break;
			}
				
			case 5:
			{
				if(OriginalPoint[0] + OriginalPoint[1] == 0.)
				{
					REAL den = 0.5 + OriginalPoint[2]*OriginalPoint[2];
					ChangedPoint[0] = (0.5*tol)/sqrt(den);
					ChangedPoint[1] = (0.5*tol)/sqrt(den);
					ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2]*tol)/sqrt(den);
				}
				break;
			}
				
			case 6:
			{
				if(OriginalPoint[0] == 1. || OriginalPoint[2] == 1.)
				{
					REAL den = 0.25 + OriginalPoint[0]*OriginalPoint[0] + OriginalPoint[2]*OriginalPoint[2];
					ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0]*tol)/sqrt(den);
					ChangedPoint[1] = (0.5*tol)/sqrt(den);
					ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2]*tol)/sqrt(den);
				}
				break;
			}
				
			case 7:
			{
				if(OriginalPoint[0]== 1. || OriginalPoint[1] == 1.)
				{
					REAL den = 0.25 + OriginalPoint[0]*OriginalPoint[0] + OriginalPoint[1]*OriginalPoint[1];
					ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0]*tol)/sqrt(den);
					ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1]*tol)/sqrt(den);
					ChangedPoint[2] = (0.5*tol)/sqrt(den);
				}
				break;
			}
				
			case 8:
			{
				if(OriginalPoint[0] + OriginalPoint[2] == 0.)
				{
					REAL den = 0.5 + OriginalPoint[1]*OriginalPoint[1];
					ChangedPoint[0] = (0.5*tol)/sqrt(den);
					ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1]*tol)/sqrt(den);
					ChangedPoint[2] = (0.5*tol)/sqrt(den);
				}
				break;
			}
				
			case 9:
			{
				if(OriginalPoint[1] + OriginalPoint[2] == 0.)
				{
					REAL den = 0.5 + OriginalPoint[0]*OriginalPoint[0];
					ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0]*tol)/sqrt(den);
					ChangedPoint[1] = (0.5*tol)/sqrt(den);
					ChangedPoint[2] = (0.5*tol)/sqrt(den);
				}
				break;
			}
				
			case 10:
			{
				if(OriginalPoint[2] == 1.)
				{
					ChangedPoint[0] = tol/sqrt(11.);
					ChangedPoint[1] = tol/sqrt(11.);
					ChangedPoint[2] = 1. - (3.*tol)/sqrt(11.);
				}
				break;
			}
				
			case 11:
			{
				if(OriginalPoint[1] == 1.)
				{
					ChangedPoint[0] = tol/sqrt(11.);
					ChangedPoint[1] = 1. - (3.*tol)/sqrt(11.);
					ChangedPoint[2] = tol/sqrt(11.);
				}
				break;
			}
				
			case 12:
			{
				if(OriginalPoint[0] + OriginalPoint[1] + OriginalPoint[2] == 0.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = tol;
					ChangedPoint[2] = tol;
				}
				break;
			}
				
			case 13:
			{
				if(OriginalPoint[0] == 1.)
				{
					ChangedPoint[0] = 1. - (3.*tol)/sqrt(11.);
					ChangedPoint[1] = tol/sqrt(11.);
					ChangedPoint[2] = tol/sqrt(11.);
				}
				break;
			}
		}
	}
	
	/** Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZGeoTetrahedra::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
												 TPZVec<long>& nodeindexes,
												 int matid,
												 long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
};
