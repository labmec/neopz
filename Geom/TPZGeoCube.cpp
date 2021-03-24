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

#ifdef PZ_LOG
TPZLogger logger("pz.geom.tpzgeocube");
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {
    const double tol = pzgeom_TPZNodeRep_tol;


   
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoCube::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size) {
        TPZManVector<REAL, 3> co(3), shift(3), scale(3);
        TPZManVector<int64_t, 3> nodeindexes(8);
        for (int i = 0; i < 3; i++) {
            scale[i] = size[i] / 3.;
            shift[i] = 1. / 2. + lowercorner[i];
        }

        for (int i = 0; i < NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j = 0; j < 3; j++) {
                co[j] = shift[j] + scale[j] * co[j] + (rand() * 0.2 / RAND_MAX) - 0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        gmesh.CreateGeoElement(ECube, nodeindexes, matid, index);
    }


    

    int TPZGeoCube::ClassId() const {
        return Hash("TPZGeoCube") ^ TPZNodeRep<8, pztopology::TPZCube>::ClassId() << 1;
    }

    void TPZGeoCube::Read(TPZStream &buf, void *context) {
        TPZNodeRep<8, pztopology::TPZCube>::Read(buf, context);
    }

    void TPZGeoCube::Write(TPZStream &buf, int withclassid) const {
        TPZNodeRep<8, pztopology::TPZCube>::Write(buf, withclassid);
    }


};
