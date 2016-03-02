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
	
	void TPZGeoCube::Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
		
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
				coor = nodes.GetVal(j,i); // ??
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
#ifdef PZDEBUG
			std::stringstream sout;
			sout << "Singular Jacobian " << detjac;
			LOGPZ_ERROR(logger, sout.str())
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
			TPZManVector<long> nodeindexes(1);
			//    TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//    gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} 
		else 
			if(side > 7 && side < 20) {//side = 8 a 19 : arestas
				TPZManVector<long> nodes(2);
				nodes[0] = orig->SideNodeIndex(side,0);
				nodes[1] = orig->SideNodeIndex(side,1);
				//      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
				long index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
				TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,0)));
				TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,1)));
				TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			} 
			else 
				if(side > 19 && side < 26) {//side = 20 a 25 : faces
					TPZManVector<long> nodes(4);
					int in;
					for (in=0;in<4;in++){
						nodes[in] = orig->SideNodeIndex(side,in);
					}
					long index;
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
										   TPZVec<long>& nodeindexes,
										   int matid,
										   long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
};
