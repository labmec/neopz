/**
 * @file
 * @brief Contains the implementation of the TPZGeoPyramid methods. 
 */

#include "pzgeopyramid.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeopyramid");
#endif

#include <cmath>

using namespace pzshape;
using namespace std;

namespace pzgeom {

    const double tol = pzgeom_TPZNodeRep_tol;

//    void TPZGeoPyramid::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
//        if(fabs(pt[0])<1.e-10 && fabs(pt[1])<1.e-10 && pt[2]==1.) {
//            //para testes com transforma�es geometricas-->>Que  o que faz o RefPattern!!
//            //(0,0,1) nunca �um ponto de integra�o
//            phi(0,0)  = 0.;
//            phi(1,0)  = 0.;
//            phi(2,0)  = 0.;
//            phi(3,0)  = 0.;
//            phi(4,0)  = 1.;
//            dphi(0,0)  = -0.25;
//            dphi(1,0)  = -0.25;
//            dphi(2,0)  = -0.25;
//            dphi(0,1)  = 0.25;
//            dphi(1,1)  = -0.25;
//            dphi(2,1)  = -0.25;
//            dphi(0,2)  = 0.25;
//            dphi(1,2)  = 0.25;
//            dphi(2,2)  = -0.25;
//            dphi(0,3)  = -0.25;
//            dphi(1,3)  = 0.25;
//            dphi(2,3)  = -0.25;
//            dphi(0,4)  = 0;
//            dphi(1,4)  = 0;
//            dphi(2,4)  = 1.;
//            
//            
//            
//            return;
//        }
//        
//        REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
//        REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
//        REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
//        REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
//        REAL lmez = (1.-pt[2]);
//        phi(0,0)  = T0xz*T0yz*lmez;
//        phi(1,0)  = T1xz*T0yz*lmez;
//        phi(2,0)  = T1xz*T1yz*lmez;
//        phi(3,0)  = T0xz*T1yz*lmez;
//        phi(4,0)  = pt[2];
//        REAL lmexmez = 1.-pt[0]-pt[2];
//        REAL lmeymez = 1.-pt[1]-pt[2];
//        REAL lmaxmez = 1.+pt[0]-pt[2];
//        REAL lmaymez = 1.+pt[1]-pt[2];
//        dphi(0,0) = -.25*lmeymez / lmez;
//        dphi(1,0) = -.25*lmexmez / lmez;
//        dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
//        
//        dphi(0,1) =  .25*lmeymez / lmez;
//        dphi(1,1) = -.25*lmaxmez / lmez;
//        dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
//        
//        dphi(0,2) =  .25*lmaymez / lmez;
//        dphi(1,2) =  .25*lmaxmez / lmez;
//        dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
//        
//        dphi(0,3) = -.25*lmaymez / lmez;
//        dphi(1,3) =  .25*lmexmez / lmez;
//        dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
//        
//        dphi(0,4) =  0.0;
//        dphi(1,4) =  0.0;
//        dphi(2,4) =  1.0;
//    }




    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void
    TPZGeoPyramid::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size) {
        TPZManVector<REAL, 3> co(3), shift(3), scale(3);
        TPZManVector<int64_t, 5> nodeindexes(5);
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
        gmesh.CreateGeoElement(EPiramide, nodeindexes, matid, index);
    }

    int TPZGeoPyramid::ClassId() const {
        return Hash("TPZGeoPyramid") ^ TPZNodeRep<5, pztopology::TPZPyramid>::ClassId() << 1;
    }

    void TPZGeoPyramid::Read(TPZStream &buf, void *context) {
        TPZNodeRep<5, pztopology::TPZPyramid>::Read(buf, context);
    }

    void TPZGeoPyramid::Write(TPZStream &buf, int withclassid) const {
        TPZNodeRep<5, pztopology::TPZPyramid>::Write(buf, withclassid);
    }

};
