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
	
	void TPZGeoPrism::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {

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
	
	
	void TPZGeoPrism::Jacobian(TPZFMatrix & coord, TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
		
		
#ifdef DEBUG
		int nnodes = NNodes;
		if (nnodes != 6) {
			PZError << "TPZGeoPrism.jacobian only implemented for"
			" 6 nodes, NumberOfNodes = " << nnodes << "\n";
		}
		const REAL tol = 1.e-3;
		if( param[0] < 0.-tol || param[0] > 1.+tol || param.NElements() != 3 ||
		   param[1] < 0.-tol || param[1] > 1.+tol || param[2] < -1.-tol || param[2] > 1.+tol) {
			PZError << "TPZGeoPrism::jacobian. param out of range : "
			" param.NElements() = " << param.NElements() <<
			"\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
			return;
		}
#endif
		jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
		REAL spacephi[6];
		TPZFMatrix phi(6,1,spacephi,6);
		REAL spacedphi[18];
		TPZFMatrix dphi(3,6,spacedphi,18);
		Shape(param,phi,dphi);
		jacobian.Zero();
		//  TPZGeoNode *np;
		int i,j;
		for(i=0;i<6;i++) {
			//    np = NodePtr(i);
			for(j=0;j<3;j++) {
				jacobian(j,0) += coord(j,i)*dphi(0,i);//np->Coord(j)*dphi(0,i);
				jacobian(j,1) += coord(j,i)*dphi(1,i);
				jacobian(j,2) += coord(j,i)*dphi(2,i);
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
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
			LOGPZ_ERROR(logger, sout.str())
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
	
	void TPZGeoPrism::X(TPZFMatrix & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
		REAL spacephi[6],spacedphi[18];
		int i,j;
		TPZFMatrix phi(6,1,spacephi,5);
		TPZFMatrix dphi(3,6,spacedphi,18);
		Shape(loc,phi,dphi);
		for(j=0;j<3;j++) {
			result[j] = 0.0;
			for(i=0;i<6;i++) result[j] += coord(j,i)*phi(i,0);
		}
	}
	
	bool TPZGeoPrism::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
        double zero = 1.E-5;
		
		REAL qsi = InternalPar[0]; 
		REAL eta = InternalPar[1]; 
		REAL zeta = InternalPar[2];
		if((qsi + eta) > 1. || qsi < 0. || eta < 0. || zeta < -1. || zeta > 1.)
		{
			cout << "Point (qsi,eta,zeta) = (" << qsi << "," << eta << "," << zeta << ") is out of TPZGeoPrism Master Element Range!\n";
			cout << "See TPZGeoPrism::MapToSide() method!\n";
			DebugStop();
		}
		bool regularmap = true;
		switch(side)
		{
			case 6://1D
				SidePar.Resize(1); 
				JacToSide.Resize(1,3);
				if(fabs(eta-1.) < zero)
				{
					SidePar[0] = 0.;
					JacToSide(0,0) = 1.; 
					JacToSide(0,1) = 1.; 
					JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
					SidePar[0] = 2.*qsi/(1.-eta) - 1.;
					JacToSide(0,0) = 2./(1.-eta); JacToSide(0,1) = 2.*qsi/((1.-eta)*(1.-eta)); JacToSide(0,2) = 0.;
				}
				break;
				
			case 7://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(qsi+eta) < zero)
				{
					SidePar[0] = 0.;
                    JacToSide(0,0) = 1.; 
					JacToSide(0,1) = 1.; 
					JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*qsi/(qsi + eta);
                    JacToSide(0,0) = -2.*eta/((qsi+eta)*(qsi+eta)); JacToSide(0,1) = 2.*qsi/((qsi+eta)*(qsi+eta)); JacToSide(0,2) = 0.;
				}
				break;
				
			case 8://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(qsi-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*eta/(1.-qsi);
                    JacToSide(0,0) = -2.*eta/((1.-qsi)*(1.-qsi)); JacToSide(0,1) = -2./(1.-qsi); JacToSide(0,2) = 0.;
				}
				break;
				
			case 9://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				SidePar[0] = zeta;
				JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 1.;
				break;
				
			case 10://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				SidePar[0] = zeta;
				JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 1.;
				break;
				
			case 11://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				SidePar[0] = zeta;
				JacToSide(0,0) = 0.; JacToSide(0,1) = 0.; JacToSide(0,2) = 1.;
				break;
				
			case 12://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(eta-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 2.*qsi/(1.-eta) - 1.;
                    JacToSide(0,0) = 2./(1.-eta); JacToSide(0,1) = 2.*qsi/((1.-eta)*(1.-eta)); JacToSide(0,2) = 0.;
				}
				break;
				
			case 13://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(qsi+eta) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*qsi/(qsi + eta);
                    JacToSide(0,0) = -2.*eta/((qsi+eta)*(qsi+eta)); JacToSide(0,1) = 2.*qsi/((qsi+eta)*(qsi+eta)); JacToSide(0,2) = 0.;
				}
				break;
				
			case 14://1D
				SidePar.Resize(1); JacToSide.Resize(1,3);
				if(fabs(qsi-1.) < zero)
				{
                    SidePar[0] = 0.;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*eta/(1.-qsi);
                    JacToSide(0,0) = -2.*eta/((1.-qsi)*(1.-qsi)); JacToSide(0,1) = -2./(1.-qsi); JacToSide(0,2) = 0.;
				}
				break;
				
			case 15://2D - triangle
				SidePar.Resize(2); JacToSide.Resize(2,3);
				SidePar[0] = qsi; SidePar[1] = eta;
				JacToSide(0,0) = 1.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
				JacToSide(1,0) = 0.; JacToSide(1,1) = 1.; JacToSide(1,2) = 0.;
				break;
				
			case 16://2D - quadrilateral
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(eta-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = zeta;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 1.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 2.*qsi/(1.-eta) - 1.; SidePar[1] = zeta;
                    JacToSide(0,0) = 2./(1.-eta); JacToSide(0,1) = 2.*qsi/((1.-eta)*(1.-eta)); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(0,2) = 1.;
				}
				break;
				
			case 17://2D - quadrilateral
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(qsi+eta) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = zeta;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 1.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*qsi/(qsi + eta); SidePar[1] = zeta;
                    JacToSide(0,0) = -2.*eta/((qsi+eta)*(qsi+eta)); JacToSide(0,1) = 2.*qsi/((qsi+eta)*(qsi+eta)); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(0,2) = 1.;
				}
				break;
				
			case 18://2D - quadrilateral
				SidePar.Resize(2); JacToSide.Resize(2,3);
				if(fabs(qsi-1.) < zero)
				{
                    SidePar[0] = 0.; SidePar[1] = zeta;
                    JacToSide(0,0) = 1.; JacToSide(0,1) = 1.; JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(1,2) = 1.;
					regularmap = false;
				}
				else
				{
                    SidePar[0] = 1. - 2.*eta/(1.-qsi); SidePar[1] = zeta;
                    JacToSide(0,0) = -2.*eta/((1.-qsi)*(1.-qsi)); JacToSide(0,1) = -2./(1.-qsi); JacToSide(0,2) = 0.;
                    JacToSide(1,0) = 0.; JacToSide(1,1) = 0.; JacToSide(0,2) = 1.;
				}
				break;
				
			case 19://2D - triangle
				SidePar.Resize(2); JacToSide.Resize(2,3);
				SidePar[0] = qsi; SidePar[1] = eta;
				JacToSide(0,0) = 1.; JacToSide(0,1) = 0.; JacToSide(0,2) = 0.;
				JacToSide(1,0) = 0.; JacToSide(1,1) = 1.; JacToSide(1,2) = 0.;
				break;
		}
		if(side < 6 || side > 19)
		{
			cout << "Cant compute MapToSide method in TPZGeoPrism class!\nParameter (SIDE) must be between 6 and 19!\nMethod Aborted!\n"; 
			DebugStop();
		}
		return regularmap;
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
		
		int index;
		if(side<6) {
			TPZManVector<int> nodeindexes(1);
			TPZGeoEl *gel;
			//		int nodestore [4];
			nodeindexes[0] = orig->NodeIndex(side);
			gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			result = gel;
		} else if (side > 5 && side < 15) {//side = 6 a 14 : arestas
			TPZManVector<int> nodes(2);
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
			TPZManVector<int> nodes(4);//4o = -1 para face triangular
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
											TPZVec<int>& nodeindexes,
											int matid,
											int& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
};
