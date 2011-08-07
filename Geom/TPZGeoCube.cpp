/**
 * @file
 * @brief Contains the implementation of the TPZGeoCube methods. 
 */

#include "TPZGeoCube.h"
#include "pzgeoel.h"
#include "pzshapecube.h"
#include "pzquad.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger = Logger::getLogger("pz.geom.tpzgeocube");
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	void TPZGeoCube::X(TPZFMatrix &nodes,TPZVec<REAL> & loc,TPZVec<REAL> &result){
		
		int nrow = nodes.Rows();
		int ncol = nodes.Cols();
		if(nrow != 3 || ncol != 8){//8x3 n� por linhas
			cout << "TPZGeoCube::X nodes matrix error size\n";
		}
		int i,j;
		TPZFNMatrix<12> phi(8,1);
		TPZFNMatrix<30> dphi(3,8);
		Shape(loc,phi,dphi);
		for(j=0;j<3;j++) {
			result[j] = 0.0;
			for(i=0;i<8;i++) result[j] += nodes(j,i)*phi(i,0);
		}
	}
	
	bool TPZGeoCube::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide) {
		TPZTransform Transf = pztopology::TPZCube::SideToSideTransform(TPZGeoCube::NSides - 1, side);
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		int R = Transf.Mult().Rows();
		int C = Transf.Mult().Cols();
		
		JacToSide.Resize(R,C);
		for(int i = 0; i < R; i++)
		{
			for(int j = 0; j < C; j++) JacToSide(i,j) = Transf.Mult()(i,j);
		}
		return true;
	}
	
	void TPZGeoCube::Shape(TPZVec<REAL> &pt,TPZFMatrix &phi,TPZFMatrix &dphi) {
		
		REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
		x[0] = (1.-pt[0])/2.;
		x[1] = (1.+pt[0])/2.;
		dx[0] = -0.5;
		dx[1] = 0.5;
		y[0] = (1.-pt[1])/2.;
		y[1] = (1.+pt[1])/2.;
		dy[0] = -0.5;
		dy[1] = 0.5;
		z[0] = (1.-pt[2])/2.;
		z[1] = (1.+pt[2])/2.;
		dz[0] = -0.5;
		dz[1] = 0.5;
		
		phi(0,0) = x[0]*y[0]*z[0];
		phi(1,0) = x[1]*y[0]*z[0];
		phi(2,0) = x[1]*y[1]*z[0];
		phi(3,0) = x[0]*y[1]*z[0];
		phi(4,0) = x[0]*y[0]*z[1];
		phi(5,0) = x[1]*y[0]*z[1];
		phi(6,0) = x[1]*y[1]*z[1];
		phi(7,0) = x[0]*y[1]*z[1];
		dphi(0,0) = dx[0]*y[0]*z[0];
		dphi(1,0) = x[0]*dy[0]*z[0];
		dphi(2,0) = x[0]*y[0]*dz[0];
		dphi(0,1) = dx[1]*y[0]*z[0];
		dphi(1,1) = x[1]*dy[0]*z[0];
		dphi(2,1) = x[1]*y[0]*dz[0];
		dphi(0,2) = dx[1]*y[1]*z[0];
		dphi(1,2) = x[1]*dy[1]*z[0];
		dphi(2,2) = x[1]*y[1]*dz[0];
		dphi(0,3) = dx[0]*y[1]*z[0];
		dphi(1,3) = x[0]*dy[1]*z[0];
		dphi(2,3) = x[0]*y[1]*dz[0];
		dphi(0,4) = dx[0]*y[0]*z[1];
		dphi(1,4) = x[0]*dy[0]*z[1];
		dphi(2,4) = x[0]*y[0]*dz[1];
		dphi(0,5) = dx[1]*y[0]*z[1];
		dphi(1,5) = x[1]*dy[0]*z[1];
		dphi(2,5) = x[1]*y[0]*dz[1];
		dphi(0,6) = dx[1]*y[1]*z[1];
		dphi(1,6) = x[1]*dy[1]*z[1];
		dphi(2,6) = x[1]*y[1]*dz[1];
		dphi(0,7) = dx[0]*y[1]*z[1];
		dphi(1,7) = x[0]*dy[1]*z[1];
		dphi(2,7) = x[0]*y[1]*dz[1];
	}
	
	void TPZGeoCube::Jacobian(TPZFMatrix &nodes,TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv){
		
#ifdef DEBUG
		//  if (NNodes != 8) {
		//    PZError << "TPZGeoCube.jacobian only implemented for"
		//      " 8 nodes, NumberOfNodes = " << NNodes << "\n";
		//  }
		if(param.NElements() != 3 || param[0] < -1.001 || param[0] > 1.001 ||
		   param[1] < -1.001 || param[1] > 1.001 || param[2] < -1.001 || param[2] > 1.001) {
			PZError << "TPZGeoCube.jacobian. param out of range : "
			" param.NElements() = " << param.NElements() <<
			"\nparam[0] = " << param[0] << " param[1] = " << param[1] << " param[2] = " << param[2] << "\n";
			return;
		}
#endif
		jacobian.Resize(3,3); axes.Resize(3,3); jacinv.Resize(3,3);
		int nrow = nodes.Rows();
		int ncol = nodes.Cols();
		if(nrow != 3 || ncol != 8){//8x3 n� por linhas
			cout << "TPZGeoCube::X nodes matrix error size\n";
		}
		
		TPZFNMatrix<12> phi(8,1);
		TPZFNMatrix<30> dphi(3,8);
		Shape(param,phi,dphi);
		jacobian.Zero();
		REAL coor;
		int i,j;
		for(i=0;i<8;i++) {
			for(j=0;j<3;j++) {
				coor = nodes(j,i); // ??
				jacobian(j,0) += coor*dphi(0,i);
				jacobian(j,1) += coor*dphi(1,i);
				jacobian(j,2) += coor*dphi(2,i);
			}
		}
		
		detjac = -jacobian(0,2)*jacobian(1,1)*jacobian(2,0)
        + jacobian(0,1)*jacobian(1,2)*jacobian(2,0)
        + jacobian(0,2)*jacobian(1,0)*jacobian(2,1)
        - jacobian(0,0)*jacobian(1,2)*jacobian(2,1)
        - jacobian(0,1)*jacobian(1,0)*jacobian(2,2)
		+ jacobian(0,0)*jacobian(1,1)*jacobian(2,2);
		
		if(IsZero(detjac))
		{
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
#ifdef LOG4CXX
			LOGPZ_ERROR(logger , sout.str())
#endif
			detjac = ZeroTolerance();
		}
		jacinv(0,0) = (-jacobian(1,2)*jacobian(2,1)+jacobian(1,1)*jacobian(2,2))/detjac;//-a12 a21 + a11 a22
		jacinv(0,1) = ( jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))/detjac;// a02 a21 - a01 a22
		jacinv(0,2) = (-jacobian(0,2)*jacobian(1,1)+jacobian(0,1)*jacobian(1,2))/detjac;//-a02 a11 + a01 a12
		jacinv(1,0) = ( jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))/detjac;// a12 a20 - a10 a22
		jacinv(1,1) = (-jacobian(0,2)*jacobian(2,0)+jacobian(0,0)*jacobian(2,2))/detjac;//-a02 a20 + a00 a22
		jacinv(1,2) = ( jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))/detjac;// a02 a10 - a00 a12
		jacinv(2,0) = (-jacobian(1,1)*jacobian(2,0)+jacobian(1,0)*jacobian(2,1))/detjac;//-a11 a20 + a10 a21
		jacinv(2,1) = ( jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))/detjac;// a01 a20 - a00 a21
		jacinv(2,2) = (-jacobian(0,1)*jacobian(1,0)+jacobian(0,0)*jacobian(1,1))/detjac;//-a01 a10 + a00 a11
		
		axes.Zero();
		axes(0,0) = 1.;
		axes(1,1) = 1.;
		axes(2,2) = 1.;
	}
	
	TPZGeoEl *TPZGeoCube::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc) {
		
		if(side<0 || side>26) {
			cout <<  "TPZGeoCube::CreateBCGeoEl unexpected side = " << side << endl;
			return 0;
		}
		if(side==26) {
			cout << "TPZGeoCube::CreateBCGeoEl with side = 26 not implemented\n";
			return 0;
		}
		
		if(side<8) {
			TPZManVector<int> nodeindexes(1);
			//    TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			int index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//    gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} 
		else 
			if(side > 7 && side < 20) {//side = 8 a 19 : arestas
				TPZManVector<int> nodes(2);
				nodes[0] = orig->SideNodeIndex(side,0);
				nodes[1] = orig->SideNodeIndex(side,1);
				//      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
				int index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
				TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,0)));
				TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,1)));
				TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			} 
			else 
				if(side > 19 && side < 26) {//side = 20 a 25 : faces
					TPZManVector<int> nodes(4);
					int in;
					for (in=0;in<4;in++){
						nodes[in] = orig->SideNodeIndex(side,in);
					}
					int index;
					TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
					//		TPZGeoElQ2d *gel = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
					for (in=0;in<8;in++){
						TPZGeoElSide(gel,in).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,in)));
					}
					TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
					return gel;
				} 
				else 
					PZError << "TPZGeoCube::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoCube::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										   TPZVec<int>& nodeindexes,
										   int matid,
										   int& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
};
