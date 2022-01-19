/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticCube methods. 
 */
#include "tpzquadraticcube.h"
#include "pzshapecube.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"

#include "tpzgeomid.h"

#include "pzlog.h"

#include "tpzchangeel.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadraticcube");
#endif

#include "fad.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

namespace pzgeom
{
    
template<class T>
void TPZQuadraticCube::TShape(const TPZVec<T> &par,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
	T qsi = par[0], eta = par[1], zeta = par[2];
		
	phi(0,0)   =  1./8.*((-1. + eta)*(-1. + qsi)*(-1. + zeta)*(2. + eta + qsi + zeta));
	phi(1,0)   = -1./8.*((-1. + eta)*(1. + qsi)*(-1. + zeta)*(2. + eta - qsi + zeta));
	phi(2,0)   = -1./8.*((1. + eta)*(1. + qsi)*(-2. + eta + qsi - zeta)*(-1. + zeta));
	phi(3,0)   =  1./8.*((1. + eta)*(-1. + qsi)*(-2. + eta - qsi - zeta)*(-1. + zeta));
	
	phi(4,0)   = -1./8.*((-1. + eta)*(-1. + qsi)*(2. + eta + qsi - zeta)*(1. + zeta));
	phi(5,0)   =  1./8.*((-1. + eta)*(1. + qsi)*(2. + eta - qsi - zeta)*(1. + zeta));
	phi(6,0)   =  1./8.*((1. + eta)*(1. + qsi)*(1. + zeta)*(-2. + eta + qsi + zeta));
	phi(7,0)   = -1./8.*((1. + eta)*(-1. + qsi)*(1. + zeta)*(-2. + eta - qsi + zeta));
	
	phi(8,0)   = -1./4.*((-1. + eta)*(-1. + qsi*qsi)*(-1. + zeta));
	phi(9,0)   =  1./4.*((-1. + eta*eta)*(1. + qsi)*(-1. + zeta));
	phi(10,0)  =  1./4.*((1. + eta)*(-1. + qsi*qsi)*(-1. + zeta));
	phi(11,0)  = -1./4.*((-1. + eta*eta)*(-1. + qsi)*(-1. + zeta));
	
	phi(12,0)  = -1./4.*((-1. + eta)*(-1. + qsi)*(-1. + zeta*zeta));
	phi(13,0)  =  1./4.*((-1. + eta)*(1. + qsi)*(-1. + zeta*zeta));
	phi(14,0)  = -1./4.*((1. + eta)*(1. + qsi)*(-1. + zeta*zeta));
	phi(15,0)  =  1./4.*((1. + eta)*(-1. + qsi)*(-1. + zeta*zeta));
	
	phi(16,0)  =  1./4.*((-1. + eta)*(-1. + qsi*qsi)*(1. + zeta));
	phi(17,0)  = -1./4.*((-1. + eta*eta)*(1. + qsi)*(1. + zeta));
	phi(18,0)  = -1./4.*((1. + eta)*(-1. + qsi*qsi)*(1. + zeta));
	phi(19,0)  =  1./4.*((-1. + eta*eta)*(-1. + qsi)*(1. + zeta));
	
	
	dphi(0,0)  =  1./8.*((-1. + eta)*(-1. + zeta)*(1. + eta + 2.*qsi + zeta));
	dphi(1,0)  =  1./8.*((-1. + qsi)*(-1. + zeta)*(1. + 2.*eta + qsi + zeta));
	dphi(2,0)  =  1./8.*((-1. + eta)*(-1. + qsi)*(1. + eta + qsi + 2.*zeta));
	
	dphi(0,1)  = -1./8.*((-1. + eta)*(-1. + zeta)*(1. + eta - 2.*qsi + zeta));
	dphi(1,1)  =  1./8.*((1. + qsi)*(-1. - 2.*eta + qsi - zeta)*(-1. + zeta));
	dphi(2,1)  = -1./8.*((-1. + eta)*(1. + qsi)*(1. + eta - qsi + 2.*zeta));
	
	dphi(0,2)  =  1./8.*((1. + eta)*(-1. + zeta)*(1. - eta - 2.*qsi + zeta));
	dphi(1,2)  =  1./8.*((1. + qsi)*(-1. + zeta)*(1. - 2.*eta - qsi + zeta));
	dphi(2,2)  = -1./8.*((1. + eta)*(1. + qsi)*(-1. + eta + qsi - 2.*zeta));
	
	dphi(0,3)  =  1./8.*((1. + eta)*(-1. + eta - 2.*qsi - zeta)*(-1. + zeta));
	dphi(1,3)  = -1./8.*((-1. + qsi)*(-1. + zeta)*(1. - 2.*eta + qsi + zeta));
	dphi(2,3)  =  1./8.*((1. + eta)*(-1. + qsi)*(-1. + eta - qsi - 2.*zeta));
	
	dphi(0,4)  = -1./8.*((-1. + eta)*(1. + eta + 2.*qsi - zeta)*(1. + zeta));
	dphi(1,4)  = -1./8.*((-1. + qsi)*(1. + 2.*eta + qsi - zeta)*(1. + zeta));
	dphi(2,4)  = -1./8.*((-1. + eta)*(-1. + qsi)*(1. + eta + qsi - 2.*zeta));
	
	dphi(0,5)  =  1./8.*((-1. + eta)*(1. + eta - 2.*qsi - zeta)*(1. + zeta));
	dphi(1,5)  = -1./8.*((1. + qsi)*(1. + zeta)*(-1. - 2.*eta + qsi + zeta));
	dphi(2,5)  =  1./8.*((-1. + eta)*(1. + qsi)*(1. + eta - qsi - 2.*zeta));
	
	dphi(0,6)  =  1./8.*((1. + eta)*(1. + zeta)*(-1. + eta + 2.*qsi + zeta));
	dphi(1,6)  =  1./8.*((1. + qsi)*(1. + zeta)*(-1. + 2.*eta + qsi + zeta));
	dphi(2,6)  =  1./8.*((1. + eta)*(1. + qsi)*(-1. + eta + qsi + 2.*zeta));
	
	dphi(0,7)  = -1./8.*((1. + eta)*(1. + zeta)*(-1. + eta - 2.*qsi + zeta));
	dphi(1,7)  =  1./8.*((-1. + qsi)*(1. - 2.*eta + qsi - zeta)*(1. + zeta));
	dphi(2,7)  =  1./8.*((1. + eta)*(-1. + qsi)*(1. - eta + qsi - 2.*zeta));
	
	dphi(0,8)  = -1./2.*((-1. + eta)*qsi*(-1. + zeta));
	dphi(1,8)  = -1./4.*((-1. + qsi*qsi)*(-1. + zeta));
	dphi(2,8)  = -1./4.*((-1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,9)  =  1./4.*((-1. + eta*eta)*(-1. + zeta));
	dphi(1,9)  =  1./2.*(eta*(1. + qsi)*(-1. + zeta));
	dphi(2,9)  =  1./4.*((-1. + eta*eta)*(1. + qsi));
	
	dphi(0,10) =  1./2.*((1. + eta)*qsi*(-1. + zeta));
	dphi(1,10) =  1./4.*((-1. + qsi*qsi)*(-1. + zeta));
	dphi(2,10) =  1./4.*((1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,11) = -1./4.*((-1. + eta*eta)*(-1. + zeta));
	dphi(1,11) = -1./2.*(eta*(-1. + qsi)*(-1. + zeta));
	dphi(2,11) = -1./4.*((-1. + eta*eta)*(-1. + qsi));
	
	dphi(0,12) = -1./4.*((-1. + eta)*(-1. + zeta*zeta));
	dphi(1,12) = -1./4.*((-1. + qsi)*(-1. + zeta*zeta));
	dphi(2,12) = -1./2.*((-1. + eta)*(-1. + qsi)*zeta);
	
	dphi(0,13) =  1./4.*((-1. + eta)*(-1. + zeta*zeta));
	dphi(1,13) =  1./4.*((1. + qsi)*(-1. + zeta*zeta));
	dphi(2,13) =  1./2.*((-1. + eta)*(1. + qsi)*zeta);
	
	dphi(0,14) = -1./4.*((1. + eta)*(-1. + zeta*zeta));
	dphi(1,14) = -1./4.*((1. + qsi)*(-1. + zeta*zeta));
	dphi(2,14) = -1./2.*((1. + eta)*(1. + qsi)*zeta);
	
	dphi(0,15) =  1./4.*((1. + eta)*(-1. + zeta*zeta));
	dphi(1,15) =  1./4.*((-1. + qsi)*(-1. + zeta*zeta));
	dphi(2,15) =  1./2.*((1. + eta)*(-1. + qsi)*zeta);
	
	dphi(0,16) =  1./2.*((-1. + eta)*qsi*(1. + zeta));
	dphi(1,16) =  1./4.*((-1. + qsi*qsi)*(1. + zeta));
	dphi(2,16) =  1./4.*((-1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,17) = -1./4.*((-1. + eta*eta)*(1. + zeta));
	dphi(1,17) = -1./2.*(eta*(1. + qsi)*(1. + zeta));
	dphi(2,17) = -1./4.*((-1. + eta*eta)*(1. + qsi));
	
	dphi(0,18) = -1./2.*((1. + eta)*qsi*(1. + zeta));
	dphi(1,18) = -1./4.*((-1. + qsi*qsi)*(1. + zeta));
	dphi(2,18) = -1./4.*((1. + eta)*(-1. + qsi*qsi));
	
	dphi(0,19) =  1./4.*((-1. + eta*eta)*(1. + zeta));
	dphi(1,19) =  1./2.*(eta*(-1. + qsi)*(1. + zeta));
	dphi(2,19) =  1./4.*((-1. + eta*eta)*(-1. + qsi));
}


template<class T>
void TPZQuadraticCube::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<20,T> phi(NNodes,1);
    TPZFNMatrix<60,T> dphi(3,NNodes);
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
void TPZQuadraticCube::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,3);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 20){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x20." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<3,T> phi(NNodes,1);
    TPZFNMatrix<6,T> dphi(3,NNodes);
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

void TPZQuadraticCube::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
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
    gmesh.CreateGeoElement(ECube, nodeindexes, matid, index);
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
//void TPZQuadraticCube::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
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
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        case (8):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (9):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (10):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (11):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] = -1.;
//            break;
//        }
//        case (12):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (13):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (14):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (15):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (16):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        case (17):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        case (18):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (19):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  1.;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

    int TPZQuadraticCube::ClassId() const{
        return Hash("TPZQuadraticCube") ^ TPZNodeRep<20,pztopology::TPZCube>::ClassId() << 1;
    }
};

template class TPZRestoreClass< TPZGeoElRefPattern<pzgeom::TPZQuadraticCube>>;

/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<20,pzgeom::TPZQuadraticCube>;

namespace pzgeom {
    template void TPZQuadraticCube::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticCube::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticCube::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticCube::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
