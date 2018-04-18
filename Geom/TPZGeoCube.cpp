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
			TPZManVector<int64_t> nodeindexes(1);
			//    TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			int64_t index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//    gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		} 
		else 
			if(side > 7 && side < 20) {//side = 8 a 19 : arestas
				TPZManVector<int64_t> nodes(2);
				nodes[0] = orig->SideNodeIndex(side,0);
				nodes[1] = orig->SideNodeIndex(side,1);
				//      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
				int64_t index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
				TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,0)));
				TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeCube::ContainedSideLocId(side,1)));
				TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			} 
			else 
				if(side > 19 && side < 26) {//side = 20 a 25 : faces
					TPZManVector<int64_t> nodes(4);
					int in;
					for (in=0;in<4;in++){
						nodes[in] = orig->SideNodeIndex(side,in);
					}
					int64_t index;
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
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoCube::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<int64_t,3> nodeindexes(8);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i=0; i<NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j=0; j<3; j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        CreateGeoElement(gmesh, ECube, nodeindexes, matid, index);
    }
    

	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoCube::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										   TPZVec<int64_t>& nodeindexes,
										   int matid,
										   int64_t& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
        
        int TPZGeoCube::ClassId() const{
            return Hash("TPZGeoCube") ^ TPZNodeRep<8, pztopology::TPZCube>::ClassId() << 1;
        }
        
        void TPZGeoCube::Read(TPZStream& buf, void* context) {
            TPZNodeRep<8, pztopology::TPZCube>::Read(buf,context);
        }

        void TPZGeoCube::Write(TPZStream& buf, int withclassid) const {
            TPZNodeRep<8, pztopology::TPZCube>::Write(buf,withclassid);
        }


};
