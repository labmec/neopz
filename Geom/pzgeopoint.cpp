/**
 * @file
 * @brief Contains the implementation of the TPZGeoPoint methods. 
 */

#include "pzgeopoint.h"
#include <ostream>               // for operator<<, basic_ostream, endl
#include "pzerror.h"             // for PZError
#include "pzgeoel.h"             // for TPZGeoEl
#include "pzgeoelside.h"         // for TPZGeoElSide
#include "pzgmesh.h"             // for TPZGeoMesh
#include "pzmanvector.h"         // for TPZManVector
#include "pzvec.h"               // for TPZVec
#include "tpzgeoelrefpattern.h"  // for CreateGeoElementPattern

using namespace std;

namespace pzgeom {
	
	// TPZGeoEl *TPZGeoPoint::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
	// 	if(side==0) {
	// 		TPZManVector<int64_t> nodeindexes(1);
	// 		nodeindexes[0] = orig->NodeIndex(0);
	// 		int64_t index;
	// 		TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			
	// 		TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,0));
	// 		return gel;
	// 	}
	// 	else PZError << "TPZGeoPoint::CreateBCGeoEl. Side = " << side << endl;
	// 	return 0;
	// }
	
	// /** Creates a geometric element according to the type of the father element */
	// TPZGeoEl *TPZGeoPoint::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
	// 										TPZVec<int64_t>& nodeindexes,
	// 										int matid,
	// 										int64_t& index)
	// {
	// 	return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	// }
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoPoint::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<int64_t,1> nodeindexes(1);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i = 0; i < NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            co.Resize(3, 0.);
            for (int j=0; j< co.size(); j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        gmesh.CreateGeoElement(EPoint, nodeindexes, matid, index);
    }
    
    int TPZGeoPoint::ClassId() const{
        return Hash("TPZGeoPoint") ^ TPZNodeRep<1, pztopology::TPZPoint>::ClassId() << 1;
    }
        
    void TPZGeoPoint::Read(TPZStream& buf, void* context) {
        TPZNodeRep<1, pztopology::TPZPoint>::Read(buf,context);
    }

    void TPZGeoPoint::Write(TPZStream& buf, int withclassid) const {
        TPZNodeRep<1, pztopology::TPZPoint>::Write(buf,withclassid);
    }

};
