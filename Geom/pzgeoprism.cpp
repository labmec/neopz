/**
 * @file
 * @brief Contains the implementation of the TPZGeoPrism methods. 
 */

#include "pzgeoprism.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapeprism.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoprism"));
#endif
using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	static const double tol = pzgeom_TPZNodeRep_tol;
	
	void TPZGeoPrism::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

		phi(0,0)  = .5*(1.-pt[0]-pt[1])*(1.-pt[2]);
		phi(1,0)  = .5*pt[0]*(1.-pt[2]);
		phi(2,0)  = .5*pt[1]*(1.-pt[2]);
		phi(3,0)  = .5*(1.-pt[0]-pt[1])*(1.+pt[2]);
		phi(4,0)  = .5*pt[0]*(1.+pt[2]);
		phi(5,0)  = .5*pt[1]*(1.+pt[2]);
		
		dphi(0,0) = -.5*(1.-pt[2]);
		dphi(1,0) = -.5*(1.-pt[2]);
		dphi(2,0) = -.5*(1.-pt[0]-pt[1]);
		
		dphi(0,1) =  .5*(1.-pt[2]);
		dphi(1,1) =  .0;
		dphi(2,1) = -.5*pt[0];
		
		dphi(0,2) =  .0;
		dphi(1,2) =  .5*(1.-pt[2]);
		dphi(2,2) = -.5*pt[1];
		
		dphi(0,3) = -.5*(1.+pt[2]);
		dphi(1,3) = -.5*(1.+pt[2]);
		dphi(2,3) =  .5*(1.-pt[0]-pt[1]);
		
		dphi(0,4) =  .5*(1.+pt[2]);
		dphi(1,4) =  .0;
		dphi(2,4) =  .5*pt[0];
		
		dphi(0,5) =  .0;
		dphi(1,5) =  .5*(1.+pt[2]);
		dphi(2,5) =  .5*pt[1];
		
	}
	
	
	void TPZGeoPrism::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
		jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
		REAL spacephi[6];
		TPZFMatrix<REAL> phi(6,1,spacephi,6);
		REAL spacedphi[18];
		TPZFMatrix<REAL> dphi(3,6,spacedphi,18);
		Shape(param,phi,dphi);
		jacobian.Zero();
		//  TPZGeoNode *np;
		int i,j;
		for(i=0;i<6;i++) {
			//    np = NodePtr(i);
			for(j=0;j<3;j++) {
				jacobian(j,0) += coord.GetVal(j,i)*dphi(0,i);//np->coord.GetVal(j)*dphi(0,i);
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
	
	void TPZGeoPrism::X(const TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
		REAL spacephi[6],spacedphi[18];
		int i,j;
		TPZFMatrix<REAL> phi(6,1,spacephi,5);
		TPZFMatrix<REAL> dphi(3,6,spacedphi,18);
		Shape(loc,phi,dphi);
		for(j=0;j<3;j++) {
			result[j] = 0.0;
			for(i=0;i<6;i++) result[j] += coord.GetVal(j,i)*phi(i,0);
		}
	}
	
	
	TPZGeoEl *TPZGeoPrism::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		TPZGeoEl *result = 0;
		if(side<0 || side>20) {
			cout << "TPZGeoPrism::CreateBCCompEl bad side = " << side << " not implemented\n";		
			return result;
		}
		
		if(side==20) {
			cout << "TPZGeoPrism::CreateBCCompEl with side = 20 not implemented\n";
			return result;
		}
		
		long index;
		if(side<6) {
			TPZManVector<long> nodeindexes(1);
			TPZGeoEl *gel;
			//		int nodestore [4];
			nodeindexes[0] = orig->NodeIndex(side);
			gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			result = gel;
		} else if (side > 5 && side < 15) {//side = 6 a 14 : arestas
			TPZManVector<long> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);//(TPZCompElPr3d::SideNodes[s][0]);
			nodes[1] = orig->SideNodeIndex(side,1);//NodeIndex(TPZCompElPr3d::SideNodes[s][1]);
			//TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapePrism::ContainedSideLocId(side,0)));
			//(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][0]));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapePrism::ContainedSideLocId(side,1)));
			//(TPZGeoElSide(this,TPZCompElPr3d::SideNodes[s][1]));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));//(TPZGeoElSide(this,side));
			result = gel;//->CreateCompEl(cmesh,index);
		} 
		else if (side > 14) {//side = 15 a 19 : faces
			TPZManVector<long> nodes(4);//4o = -1 para face triangular
			int iside;
			for (iside=0;iside<4;iside++){
				nodes[iside] = orig->SideNodeIndex(side,iside);
			}
			TPZGeoEl *gelt;
			TPZGeoEl *gelq;
			if(side>15 && side<19) {
				gelq = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
				//gelq = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
				
				for (iside=0;iside<8;iside++){
					TPZGeoElSide(gelq,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePrism::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gelq,8).SetConnectivity(TPZGeoElSide(orig,side));
				result = gelq;
			} else {
				nodes.Resize(3);
				gelt = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
				//			gelt = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
				int iside;
				for (iside=0;iside<6;iside++){
					TPZGeoElSide(gelt,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePrism::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gelt,6).SetConnectivity(TPZGeoElSide(orig,side));
				result = gelt;
			}
		} else {
			PZError << "TPZGeoPrism::CreateBCGeoEl. Side = " << side << endl;
			return 0;
		}
		//	cout << "\n\nBoundary element " << bc << "  created for prism side " << side << endl;
		//	result->Print(cout);
		//	cout << "\n\nPrism Element\n";
		
		
		return result;
	}
	
	void TPZGeoPrism::FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
	{
		ChangedPoint.Resize(OriginalPoint.NElements(),0.);
		ChangedPoint = OriginalPoint;
		
		switch(side)
		{
			case 6:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 1.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = 1. - 2.*tol;
				}
				break;
			}
				
			case 7:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = tol;
				}
				break;
			}
				
			case 8:
			{
				if(OriginalPoint[0] == 1. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = 1.-tol;
					ChangedPoint[1] = tol/2.;
				}
				break;
			}
				
			case 12:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 1.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = 1. - 2.*tol;
				}
				break;
			}
				
			case 13:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = tol;
				}
				break;
			}
				
			case 14:
			{
				if(OriginalPoint[0] == 1. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = 1.-tol;
					ChangedPoint[1] = tol/2.;
				}
				break;
			}
				
			case 16:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 1.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = 1. - 2.*tol;
				}
				break;
			}
				
			case 17:
			{
				if(OriginalPoint[0] == 0. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = tol;
					ChangedPoint[1] = tol;
				}
				break;
			}
				
			case 18:
			{
				if(OriginalPoint[0] == 1. && OriginalPoint[1] == 0.)
				{
					ChangedPoint[0] = 1.-tol;
					ChangedPoint[1] = tol/2.;
				}
				break;
			}
		}
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoPrism::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											TPZVec<long>& nodeindexes,
											int matid,
											long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
};
