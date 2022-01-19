/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticQuad methods. 
 */
#include "tpzquadraticquad.h"
#include "pzshapequad.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"
#include "tpzgeomid.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadraticquad");
#endif

#include "fad.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

template<class T>
void TPZQuadraticQuad::TShape(const TPZVec<T> &par,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
    
	T qsi = par[0], eta = par[1];
	
	phi(0,0)  = -0.25*(-1. + eta)*(-1. + qsi)*(1. + eta + qsi);
	phi(1,0)  =  0.25*(-1. + eta)*(1. + eta - qsi)*(1. + qsi);
	phi(2,0)  =  0.25*( 1. + qsi)*(1. + eta)*(qsi + eta - 1.);
	phi(3,0)  = -0.25*( 1. + eta)*(-1. + eta - qsi)*(-1. + qsi);
	
	phi(4,0)  =  0.5*(-1. + eta)*(-1. + qsi)*( 1. + qsi);
	phi(5,0)  = -0.5*(-1. + eta)*( 1. + eta)*( 1. + qsi);
	phi(6,0)  = -0.5*( 1. + eta)*(-1. + qsi)*( 1. + qsi);
	phi(7,0)  =  0.5*(-1. + eta)*( 1. + eta)*(-1. + qsi);
	
	dphi(0,0) = -0.25*(-1. + eta)*(eta + 2.*qsi);
	dphi(1,0) = -0.25*(-1. + qsi)*(2.*eta + qsi);
	dphi(0,1) =  0.25*(-1. + eta)*(eta - 2.*qsi);
	dphi(1,1) =  0.25*(2.*eta - qsi)*(1. + qsi);
	dphi(0,2) =  0.25*(1. + eta)*(eta + 2.*qsi);
	dphi(1,2) =  0.25*(1. + qsi)*(2.*eta + qsi);
	dphi(0,3) = -0.25*(1. + eta)*(eta - 2.*qsi);
	dphi(1,3) = -0.25*(2.*eta - qsi)*(-1. + qsi);
	
	dphi(0,4) =  (-1. + eta)*qsi;
	dphi(1,4) =  0.5*(-1. + qsi)*(1. + qsi);
	dphi(0,5) = -0.5*(-1. + eta)*(1. + eta);
	dphi(1,5) = -eta*(1. + qsi);
	dphi(0,6) = -(1. + eta)*qsi;
	dphi(1,6) = -0.5*(-1. + qsi)*(1. + qsi);
	dphi(0,7) =  0.5*(-1. + eta)*(1. + eta);
	dphi(1,7) =  eta*(-1. + qsi);
}


template<class T>
void TPZQuadraticQuad::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<4,T> phi(NNodes,1);
    TPZFNMatrix<8,T> dphi(2,NNodes);
    TShape(loc,phi,dphi);
    int space = nodes.Rows();
    
    for(int i = 0; i < space; i++) {
        x[i] = 0.0;
        for(int j = 0; j < NNodes; j++) {
            x[i] += phi(j,0)*nodes.GetVal(i,j);
        }
    }
    
}

template<class T>
inline void TPZQuadraticQuad::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,2);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 8){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x8." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<3,T> phi(NNodes,1);
    TPZFNMatrix<6,T> dphi(2,NNodes);
    TShape(loc,phi,dphi);
    for(int i = 0; i < NNodes; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
            gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
        }
    }
    
}



/// create an example element based on the topology
/* @param gmesh mesh in which the element should be inserted
 @param matid material id of the element
 @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
 @param size (in) size of space where the element should be created
 */
#include "tpzchangeel.h"

void TPZQuadraticQuad::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZManVector<REAL,3> co(3),shift(3),scale(3);
    TPZManVector<int64_t,4> nodeindexes(NCornerNodes);
    for (int i=0; i<3; i++) {
        scale[i] = size[i]/3.;
        shift[i] = size[i]/2.+lowercorner[i];
    }
    
    for (int i=0; i<NCornerNodes; i++) {
        ParametricDomainNodeCoord(i, co);
        co.Resize(3,0.);
        for (int j=0; j<3; j++) {
            co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
    }
    int64_t index;
    gmesh.CreateGeoElement(EQuadrilateral, nodeindexes, matid, index);
    TPZGeoEl *gel = gmesh.Element(index);
    int nsides = gel->NSides();
    for (int is=0; is<nsides; is++) {
        gel->SetSideDefined(is);
    }
    gel = TPZChangeEl::ChangeToQuadratic(&gmesh, index);
    for (int node = gel->NCornerNodes(); node < gel->NNodes(); node++) {
        TPZManVector<REAL,3> co(3);
        gel->NodePtr(node)->GetCoordinates(co);
        for (int i=0; i<3; i++) {
            co[i] += (0.2*rand())/RAND_MAX - 0.1;
        }
        gel->NodePtr(node)->SetCoord(co);
    }
}

//void TPZQuadraticQuad::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] = -1.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 1.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] = -1.;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  0.;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticQuad

int TPZQuadraticQuad::ClassId() const{
    return Hash("TPZQuadraticQuad") ^ pzgeom::TPZNodeRep<8,pztopology::TPZQuadrilateral>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticQuad>>;

/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<8,TPZQuadraticQuad>;

namespace pzgeom {
    template void TPZQuadraticQuad::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticQuad::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticQuad::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticQuad::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

    template void TPZQuadraticQuad::TShape(const TPZVec<REAL> &param, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

}
