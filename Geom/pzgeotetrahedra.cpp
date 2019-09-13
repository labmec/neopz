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


    void TPZGeoTetrahedra::FixSingularity(int side, TPZVec<REAL> &OriginalPoint, TPZVec<REAL> &ChangedPoint) {
        ChangedPoint.Resize(OriginalPoint.NElements(), 0.);
        ChangedPoint = OriginalPoint;

        switch (side) {
            case 4: {
                if (OriginalPoint[1] == 1. || OriginalPoint[2] == 1.) {
                    REAL den = 0.25 + OriginalPoint[1] * OriginalPoint[1] + OriginalPoint[2] * OriginalPoint[2];
                    ChangedPoint[0] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1] * tol) / sqrt(den);
                    ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2] * tol) / sqrt(den);
                }
                break;
            }

            case 5: {
                if (OriginalPoint[0] + OriginalPoint[1] == 0.) {
                    REAL den = 0.5 + OriginalPoint[2] * OriginalPoint[2];
                    ChangedPoint[0] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[1] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2] * tol) / sqrt(den);
                }
                break;
            }

            case 6: {
                if (OriginalPoint[0] == 1. || OriginalPoint[2] == 1.) {
                    REAL den = 0.25 + OriginalPoint[0] * OriginalPoint[0] + OriginalPoint[2] * OriginalPoint[2];
                    ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0] * tol) / sqrt(den);
                    ChangedPoint[1] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[2] = OriginalPoint[2] - (OriginalPoint[2] * tol) / sqrt(den);
                }
                break;
            }

            case 7: {
                if (OriginalPoint[0] == 1. || OriginalPoint[1] == 1.) {
                    REAL den = 0.25 + OriginalPoint[0] * OriginalPoint[0] + OriginalPoint[1] * OriginalPoint[1];
                    ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0] * tol) / sqrt(den);
                    ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1] * tol) / sqrt(den);
                    ChangedPoint[2] = (0.5 * tol) / sqrt(den);
                }
                break;
            }

            case 8: {
                if (OriginalPoint[0] + OriginalPoint[2] == 0.) {
                    REAL den = 0.5 + OriginalPoint[1] * OriginalPoint[1];
                    ChangedPoint[0] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[1] = OriginalPoint[1] - (OriginalPoint[1] * tol) / sqrt(den);
                    ChangedPoint[2] = (0.5 * tol) / sqrt(den);
                }
                break;
            }

            case 9: {
                if (OriginalPoint[1] + OriginalPoint[2] == 0.) {
                    REAL den = 0.5 + OriginalPoint[0] * OriginalPoint[0];
                    ChangedPoint[0] = OriginalPoint[0] - (OriginalPoint[0] * tol) / sqrt(den);
                    ChangedPoint[1] = (0.5 * tol) / sqrt(den);
                    ChangedPoint[2] = (0.5 * tol) / sqrt(den);
                }
                break;
            }

            case 10: {
                if (OriginalPoint[2] == 1.) {
                    ChangedPoint[0] = tol / sqrt(11.);
                    ChangedPoint[1] = tol / sqrt(11.);
                    ChangedPoint[2] = 1. - (3. * tol) / sqrt(11.);
                }
                break;
            }

            case 11: {
                if (OriginalPoint[1] == 1.) {
                    ChangedPoint[0] = tol / sqrt(11.);
                    ChangedPoint[1] = 1. - (3. * tol) / sqrt(11.);
                    ChangedPoint[2] = tol / sqrt(11.);
                }
                break;
            }

            case 12: {
                if (OriginalPoint[0] + OriginalPoint[1] + OriginalPoint[2] == 0.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = tol;
                    ChangedPoint[2] = tol;
                }
                break;
            }

            case 13: {
                if (OriginalPoint[0] == 1.) {
                    ChangedPoint[0] = 1. - (3. * tol) / sqrt(11.);
                    ChangedPoint[1] = tol / sqrt(11.);
                    ChangedPoint[2] = tol / sqrt(11.);
                }
                break;
            }
        }
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
        gmesh.CreateGeoElement(ETetraedro, nodeindexes, matid, index);
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