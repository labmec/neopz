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
