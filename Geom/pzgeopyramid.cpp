/**
 * @file
 * @brief Contains the implementation of the TPZGeoPyramid methods. 
 */

#include "pzgeopyramid.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeopyramid"));
#endif

#include <cmath>

using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	const double tol = pzgeom_TPZNodeRep_tol;
	
	void TPZGeoPyramid::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		if(fabs(pt[0])<1.e-10 && fabs(pt[1])<1.e-10 && pt[2]==1.) {
			//para testes com transforma�es geometricas-->>Que  o que faz o RefPattern!!
			//(0,0,1) nunca �um ponto de integra�o
			phi(0,0)  = 0.;
			phi(1,0)  = 0.;
			phi(2,0)  = 0.;
			phi(3,0)  = 0.;
			phi(4,0)  = 1.;
			dphi(0,0)  = -0.25;
			dphi(1,0)  = -0.25;
			dphi(2,0)  = -0.25;
			dphi(0,1)  = 0.25;
			dphi(1,1)  = -0.25;
			dphi(2,1)  = -0.25;
			dphi(0,2)  = 0.25;
			dphi(1,2)  = 0.25;
			dphi(2,2)  = -0.25;
			dphi(0,3)  = -0.25;
			dphi(1,3)  = 0.25;
			dphi(2,3)  = -0.25;
			dphi(0,4)  = 0;
			dphi(1,4)  = 0;
			dphi(2,4)  = 1.;
			
			
			
			return;
		}
		
		REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
		REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
		REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
		REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
		REAL lmez = (1.-pt[2]);
		phi(0,0)  = T0xz*T0yz*lmez;
		phi(1,0)  = T1xz*T0yz*lmez;
		phi(2,0)  = T1xz*T1yz*lmez;
		phi(3,0)  = T0xz*T1yz*lmez;
		phi(4,0)  = pt[2];
		REAL lmexmez = 1.-pt[0]-pt[2];
		REAL lmeymez = 1.-pt[1]-pt[2];
		REAL lmaxmez = 1.+pt[0]-pt[2];
		REAL lmaymez = 1.+pt[1]-pt[2];
		dphi(0,0) = -.25*lmeymez / lmez;
		dphi(1,0) = -.25*lmexmez / lmez;
		dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
		
		dphi(0,1) =  .25*lmeymez / lmez;
		dphi(1,1) = -.25*lmaxmez / lmez;
		dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
		
		dphi(0,2) =  .25*lmaymez / lmez;
		dphi(1,2) =  .25*lmaxmez / lmez;
		dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
		
		dphi(0,3) = -.25*lmaymez / lmez;
		dphi(1,3) =  .25*lmexmez / lmez;
		dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
		
		dphi(0,4) =  0.0;
		dphi(1,4) =  0.0;
		dphi(2,4) =  1.0;
	}
	
	void TPZGeoPyramid::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
        jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
		REAL spacephi[5];
		TPZFMatrix<REAL> phi(5,1,spacephi,5);
		REAL spacedphi[15];
		TPZFMatrix<REAL> dphi(3,5,spacedphi,15);
		Shape(param,phi,dphi);
		jacobian.Zero();
		
		int i,j;
		for(i=0;i<5;i++) {
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
	
	void TPZGeoPyramid::X(const TPZFMatrix<REAL> & coord, TPZVec<REAL> & loc,TPZVec<REAL> &result){
		REAL spacephi[10],spacedphi[20];
		int i,j;
		TPZFMatrix<REAL> phi(5,1,spacephi,10);
		TPZFMatrix<REAL> dphi(3,5,spacedphi,20);
		Shape(loc,phi,dphi);
		for(j=0;j<3;j++) {
			result[j] = 0.0;
			for(i=0;i<5;i++) result[j] += coord.GetVal(j,i)*phi(i,0);
		}
	}
	
	
	TPZGeoEl *TPZGeoPyramid::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		if(side<0 || side>18) {
			cout << "TPZGeoPyramid::CreateBCGeoEl Bad parameter side = " 
			<< side << "not implemented\n";
			return 0;
		}
		
		if(side==18) {
			cout << "TPZGeoPyramid::CreateBCCompEl with side = 18 not implemented\n";
			return 0;
		}
		
		if(side<5) {
			TPZManVector<long> nodeindexes(1);
			//		TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide origside(orig,side);
			TPZGeoElSide(gel,0).SetConnectivity(origside);
			return gel;
		} 
		else if (side > 4 && side < 13) {//side =5 a 12 : lados
			TPZManVector<long> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			//		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if (side > 12) {//side = 13 a 17 : faces
			TPZManVector<long> nodes(4);//4o = -1 para face triangular
			int iside;
			for (iside=0;iside<4;iside++){
				nodes[iside] = orig->SideNodeIndex(side,iside);
			}
			if(side==13) {
				long index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
				//      		gelq = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
				for (iside=0; iside<8; iside++){
					TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			} 
			else {
				nodes.Resize(3);
				long index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
				//			gelt = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
				for (iside=0; iside<6; iside++){
					TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			}
		} 
		else 
			PZError << "TPZGeoPyramid::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	void TPZGeoPyramid::FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
	{
		ChangedPoint.Resize(OriginalPoint.NElements(),0.);
		ChangedPoint = OriginalPoint;
		
		switch(side)
		{
			case 5:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 6:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 7:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 8:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 9:
			{
				if( ChangedPoint[0] == 1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == 1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 10:
			{
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 11:
			{
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] ==  1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 12:
			{
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] ==  1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 13:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 14:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  1.) )
				{
					ChangedPoint[1] = 1. - tol;
				}
				break;
			}
				
			case 15:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  -1.) )
				{
					ChangedPoint[0] = -1. + tol;
				}
				break;
			}
				
			case 16:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] == -1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] == -1.) )
				{
					ChangedPoint[1] = -1. + tol;
				}
				break;
			}
				
			case 17:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  -1.) )
				{
					ChangedPoint[0] = 1. - tol;
				}
				break;
			}
		}
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoPyramid::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											  TPZVec<long>& nodeindexes,
											  int matid,
											  long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
};
