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
	
    void TPZGeoCube::Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv){
        
        DebugStop();
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
