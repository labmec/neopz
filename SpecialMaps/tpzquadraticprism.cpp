/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticPrism methods. 
 */
#include "tpzquadraticprism.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"
#include "pzshapepiram.h"
#include "tpzgeomid.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadraticprism");
#endif

#include "fad.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

template<class T>
void TPZQuadraticPrism::TShape(const TPZVec<T> &par,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
	
	T qsi = par[0], eta = par[1], zeta = par[2];
	
	phi(0,0)   = -0.5*(-1. + eta + qsi)*(-1. + zeta)*(2.*(eta + qsi) + zeta);
	phi(1,0)   =  0.5*qsi*(-1. + zeta)*(2. - 2.*qsi + zeta);
	phi(2,0)   =  0.5*eta*(-1. + zeta)*(2. - 2.*eta + zeta);
	phi(3,0)   = -0.5*(-1. + eta + qsi)*(1. + zeta)*(-2.*(eta + qsi) + zeta);
	phi(4,0)   =  0.5*qsi*(1. + zeta)*(-2. + 2.*qsi + zeta);
	phi(5,0)   =  0.5*eta*(1. + zeta)*(-2. + 2.*eta + zeta);
	phi(6,0)   =  2.*qsi*(-1. + eta + qsi)*(-1. + zeta);
	phi(7,0)   = -2.*eta*qsi*(-1. + zeta);
	phi(8,0)   =  2.*eta*(-1. + eta + qsi)*(-1. + zeta);
	phi(9,0)   =  (-1. + eta + qsi)*(-1. + zeta*zeta);
	phi(10,0)  =  qsi - qsi*zeta*zeta;
	phi(11,0)  =  eta - eta*zeta*zeta;
	phi(12,0)  = -2.*qsi*(-1. + eta + qsi)*(1. + zeta);
	phi(13,0)  =  2.*eta*qsi*(1. + zeta);
    phi(14,0)  = -2.*eta*(-1. + eta + qsi)*(1. + zeta);

	dphi(0,0)  =  (1. - 2.*eta - 2.*qsi - 0.5*zeta)*(-1. + zeta);
	dphi(1,0)  =  (1. - 2.*eta - 2.*qsi - 0.5*zeta)*(-1. + zeta);
	dphi(2,0)  =  (-1. + eta + qsi)*(0.5 - eta - qsi - zeta);
	
	dphi(0,1)  = -1. + qsi*(2. - 2.*zeta) + 0.5*zeta + 0.5*zeta*zeta;
	dphi(1,1)  =  0.;
	dphi(2,1)  =  qsi*(0.5 - 1.*qsi + 1.*zeta);
	
	dphi(0,2)  =  0.;
	dphi(1,2)  = -1. + eta*(2. - 2.*zeta) + 0.5*zeta + 0.5*zeta*zeta;
	dphi(2,2)  =  eta*(0.5 - eta + zeta);
	
	dphi(0,3)  =  (-1. + 2.*eta + 2.*qsi - 0.5*zeta)*(1. + zeta);
	dphi(1,3)  =  (-1. + 2.*eta + 2.*qsi - 0.5*zeta)*(1. + zeta);
	dphi(2,3)  =  (-1. + eta + qsi)*(-0.5 + eta + qsi - zeta);
	
	dphi(0,4)  =  (-1. + 2.*qsi + 0.5*zeta)*(1. + zeta);
	dphi(1,4)  =  0.;
	dphi(2,4)  =  qsi*(-0.5 + qsi + zeta);
	
	dphi(0,5)  =  0.;
	dphi(1,5)  = -1. - 0.5*zeta + 0.5*zeta*zeta + eta*(2. + 2.*zeta);
	dphi(2,5)  =  eta*(-0.5 + eta + zeta);
	
	dphi(0,6)  =  (-2. + 2.*eta + 4.*qsi)*(-1. + zeta);
	dphi(1,6)  =  2.*qsi*(-1. + zeta);
	dphi(2,6)  =  2.*qsi*(-1. + eta + qsi);
	
	dphi(0,7)  = -2.*eta*(-1. + zeta);
	dphi(1,7)  = -2.*qsi*(-1. + zeta);
	dphi(2,7)  = -2.*eta*qsi;
	
	dphi(0,8)  =  2.*eta*(-1. + zeta);
	dphi(1,8)  =  (-2. + 4.*eta + 2.*qsi)*(-1. + zeta);
	dphi(2,8)  =  2.*eta*(-1. + eta + qsi);
	
	dphi(0,9)  = -1. + zeta*zeta;
	dphi(1,9)  = -1. + zeta*zeta;
	dphi(2,9)  =  2.*(-1. + eta + qsi)*zeta;
	
	dphi(0,10) =  1. - zeta*zeta;
	dphi(1,10) =  0.;
	dphi(2,10) = -2.*qsi*zeta;
	
	dphi(0,11) =  0.;
	dphi(1,11) =  1. - zeta*zeta;
	dphi(2,11) = -2.*eta*zeta;
	
	dphi(0,12) =  (2. - 2.*eta - 4.*qsi)*(1. + zeta);
	dphi(1,12) = -2.*qsi*(1. + zeta);
	dphi(2,12) = -2.*qsi*(-1. + eta + qsi);
    
    dphi(0,13) =  2.*eta*(1. + zeta);
	dphi(1,13) =  2.*qsi*(1. + zeta);
	dphi(2,13) =  2.*eta*qsi;
    
    dphi(0,14) = -2.*eta*(1. + zeta);
	dphi(1,14) =  (2. - 4.*eta - 2.*qsi)*(1. + zeta);
	dphi(2,14) = -2.*eta*(-1. + eta + qsi);
}

template<class T>
void TPZQuadraticPrism::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<15,T> phi(NNodes,1);
    TPZFNMatrix<45,T> dphi(3,NNodes);
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
void TPZQuadraticPrism::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,3);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 15){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x15." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<15,T> phi(NNodes,1);
    TPZFNMatrix<45,T> dphi(3,NNodes);
    TShape(loc,phi,dphi);
    for(int i = 0; i < NNodes; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
            gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
            gradx(j,2) += nodes.GetVal(j,i)*dphi(2,i);
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

void TPZQuadraticPrism::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
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
    gmesh.CreateGeoElement(EPrisma, nodeindexes, matid, index);
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

//void TPZQuadraticPrism::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.0;
//            nodeCoord[2] = -1.0;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] = -1.0;
//            break;
//        }
//        case (8):
//        {
//            nodeCoord[0] =  0.0;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] = -1.0;
//            break;
//        }
//        case (9):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (10):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (11):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (12):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.0;
//            nodeCoord[2] =  1.0;
//            break;
//        }
//        case (13):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  1.0;
//            break;
//        }
//        case (14):
//        {
//            nodeCoord[0] =  0.0;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  1.0;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticPrism

int TPZQuadraticPrism::ClassId() const{
    return Hash("TPZQuadraticPrism") ^ TPZNodeRep<15,pztopology::TPZPrism>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticPrism>>;

/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<15,TPZQuadraticPrism>;

namespace pzgeom {
    template void TPZQuadraticPrism::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticPrism::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticPrism::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticPrism::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
