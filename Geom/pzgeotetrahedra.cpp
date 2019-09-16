/**
 * @file
 * @brief Contains the implementation of the TPZGeoTetrahedra methods. 
 */

#include "pzgeotetrahedra.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "tpzgeoelrefpattern.h"

using namespace pzshape;
using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeotetrahedra"));
#endif

namespace pzgeom {

    const double tol = pzgeom_TPZNodeRep_tol;



    TPZGeoEl *TPZGeoTetrahedra::CreateBCGeoEl(TPZGeoEl *orig, int side, int bc) {
        if (side < 0 || side > 14) {
            cout << "TPZGeoTetrahedra::CreateBCCompEl with bad side = " << side << "not implemented\n";
            return 0;
        }

        if (side == 14) {
            cout << "TPZGeoTetrahedra::CreateBCCompEl with side = 14 not implemented\n";
            return 0;
        }
        if (side < 4) {
            TPZManVector<int64_t> nodeindexes(1);
            //		TPZGeoElPoint *gel;
            nodeindexes[0] = orig->NodeIndex(side);
            int64_t index;
            TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint, nodeindexes, bc, index);
            //		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
            TPZGeoElSide(gel, 0).SetConnectivity(TPZGeoElSide(orig, side));
            return gel;
        } else if (side > 3 && side < 10) {//side =4 a 9 : lados
            TPZManVector<int64_t> nodes(2);
            nodes[0] = orig->SideNodeIndex(side, 0);
            nodes[1] = orig->SideNodeIndex(side, 1);
            int64_t index;
            TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned, nodes, bc, index);
            //		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
            TPZGeoElSide(gel, 0).SetConnectivity(TPZGeoElSide(orig, TPZShapeTetra::ContainedSideLocId(side, 0)));
            TPZGeoElSide(gel, 1).SetConnectivity(TPZGeoElSide(orig, TPZShapeTetra::ContainedSideLocId(side, 1)));
            TPZGeoElSide(gel, 2).SetConnectivity(TPZGeoElSide(orig, side));
            return gel;
        } else if (side > 9) {//side = 10 a 13 : faces
            TPZManVector<int64_t> nodes(3);
            int in;
            for (in = 0; in < 3; in++) {
                nodes[in] = orig->SideNodeIndex(side, in);
            }
            int64_t index;
            TPZGeoEl *gel = orig->CreateGeoElement(ETriangle, nodes, bc, index);
            //		TPZGeoElT2d *gel = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
            for (in = 0; in < 6; in++) {
                TPZGeoElSide(gel, in).SetConnectivity(TPZGeoElSide(orig, TPZShapeTetra::ContainedSideLocId(side, in)));
            }
            TPZGeoElSide(gel, 6).SetConnectivity(TPZGeoElSide(orig, side));
            return gel;
        } else
            PZError << "TPZGeoTetrahedra::CreateBCGeoEl. Side = " << side << endl;
        return 0;
    }

    /** Creates a geometric element according to the type of the father element */
    TPZGeoEl *TPZGeoTetrahedra::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                                 TPZVec<int64_t> &nodeindexes,
                                                 int matid,
                                                 int64_t &index) {
        return CreateGeoElementPattern(mesh, type, nodeindexes, matid, index);
    }

    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoTetrahedra::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner,
                                                TPZVec<REAL> &size) {
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
        CreateGeoElement(gmesh, ETetraedro, nodeindexes, matid, index);
    }

    int TPZGeoTetrahedra::ClassId() const {
        return Hash("TPZGeoTetrahedra") ^ TPZNodeRep<4, pztopology::TPZTetrahedron>::ClassId() << 1;
    }

    void TPZGeoTetrahedra::Read(TPZStream &buf, void *context) {
        TPZNodeRep<4, pztopology::TPZTetrahedron>::Read(buf, context);
    }

    void TPZGeoTetrahedra::Write(TPZStream &buf, int withclassid) const {
        TPZNodeRep<4, pztopology::TPZTetrahedron>::Write(buf, withclassid);
    }


};