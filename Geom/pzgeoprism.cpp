/**
 * @file
 * @brief Contains the implementation of the TPZGeoPrism methods. 
 */

#include "pzgeoprism.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapeprism.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoprism"));
#endif
using namespace pzshape;
using namespace std;

namespace pzgeom {

    static const double tol = pzgeom_TPZNodeRep_tol;



    void TPZGeoPrism::FixSingularity(int side, TPZVec<REAL> &OriginalPoint, TPZVec<REAL> &ChangedPoint) {
        ChangedPoint.Resize(OriginalPoint.NElements(), 0.);
        ChangedPoint = OriginalPoint;

        switch (side) {
            case 6: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 1.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = 1. - 2. * tol;
                }
                break;
            }

            case 7: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = tol;
                }
                break;
            }

            case 8: {
                if (OriginalPoint[0] == 1. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = 1. - tol;
                    ChangedPoint[1] = tol / 2.;
                }
                break;
            }

            case 12: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 1.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = 1. - 2. * tol;
                }
                break;
            }

            case 13: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = tol;
                }
                break;
            }

            case 14: {
                if (OriginalPoint[0] == 1. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = 1. - tol;
                    ChangedPoint[1] = tol / 2.;
                }
                break;
            }

            case 16: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 1.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = 1. - 2. * tol;
                }
                break;
            }

            case 17: {
                if (OriginalPoint[0] == 0. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = tol;
                    ChangedPoint[1] = tol;
                }
                break;
            }

            case 18: {
                if (OriginalPoint[0] == 1. && OriginalPoint[1] == 0.) {
                    ChangedPoint[0] = 1. - tol;
                    ChangedPoint[1] = tol / 2.;
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
    void
    TPZGeoPrism::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size) {
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
        gmesh.CreateGeoElement(EPrisma, nodeindexes, matid, index);
    }

    int TPZGeoPrism::ClassId() const {
        return Hash("TPZGeoPrism") ^ TPZNodeRep<6, pztopology::TPZPrism>::ClassId() << 1;
    }

    void TPZGeoPrism::Read(TPZStream &buf, void *context) {
        TPZNodeRep<6, pztopology::TPZPrism>::Read(buf, context);
    }

    void TPZGeoPrism::Write(TPZStream &buf, int withclassid) const {
        TPZNodeRep<6, pztopology::TPZPrism>::Write(buf, withclassid);
    }

};
