/**
 * @file
 * @brief Contains the implementation of the TPZGeoLinear methods.
 */

#include "TPZGeoLinear.h"
#include "pzquad.h"
#include "pzshapelinear.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeolinear");
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {

    // TPZGeoEl * TPZGeoLinear::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
    //     if(side==2) {
    //         TPZManVector<int64_t> nodes(2);
    //         nodes[0] = orig->SideNodeIndex(side,0);
    //         nodes[1] = orig->SideNodeIndex(side,1);
    //         int64_t index;
    //         TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
    //         TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,0)));
    //         TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,1)));
    //         TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
    //         return gel;
    //     }
    //     else if(side==0 || side==1) {
    //         TPZManVector<int64_t> nodeindexes(1);
    //         nodeindexes[0] = orig->SideNodeIndex(side,0);
    //         int64_t index;
    //         TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
    //         TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
    //         return gel;
    //     }
    //     else {
    //         PZError << "TPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
    //     }
        
    //     return 0;
    // }
    
    // TPZGeoEl * TPZGeoLinear::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
    //                                          TPZVec<int64_t>& nodeindexes,
    //                                          int matid,
    //                                          int64_t& index)
    // {
    //     return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
    // }
    
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoLinear::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<int64_t,3> nodeindexes(2);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i=0; i<NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            int j;
            for (j=0; j<co.size(); j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            co.Resize(3);
            for (; j<3; j++) {
                co[j] = shift[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        gmesh.CreateGeoElement(EOned, nodeindexes, matid, index);
    }
    
    int TPZGeoLinear::ClassId() const{
        return Hash("TPZGeoLinear") ^ TPZNodeRep<2, pztopology::TPZLine>::ClassId() << 1;
    }

    void TPZGeoLinear::Read(TPZStream& buf, void* context) {
        TPZNodeRep<2, pztopology::TPZLine>::Read(buf, context);
    }

    void TPZGeoLinear::Write(TPZStream& buf, int withclassid) const {
        TPZNodeRep<2, pztopology::TPZLine>::Write(buf, withclassid);
    }

    
}
