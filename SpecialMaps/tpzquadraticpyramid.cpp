/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticPyramid methods. 
 */
#include "tpzquadraticpyramid.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"
#include "pzshapepiram.h"
#include "tpzgeomid.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadraticpyramid");
#endif

#include "fad.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

template<class T>
void TPZQuadraticPyramid::TShape(const TPZVec<T> &par,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
    
    T qsi = par[0], eta = par[1], zeta = par[2];
    T check = zeta-1.;
    if(fabs(check) < 1.E-6)
    {
        phi(0,0)  = 0.;
        phi(1,0)  = 0.;
        phi(2,0)  = 0.;
        phi(3,0)  = 0.;
        phi(4,0)  = 1.;
        phi(5,0)  = 0.;
        phi(6,0)  = 0.;
        phi(7,0)  = 0.;
        phi(8,0)  = 0.;
        phi(9,0)  = 0.;
        phi(10,0) = 0.;
        phi(11,0) = 0.;
        phi(12,0) = 0.;
        
        dphi(0,0) = 0.;
        dphi(1,0) = 0.;
        dphi(2,0) = -0.25;
        
        dphi(0,1) = 0.;
        dphi(1,1) = 0.;
        dphi(2,1) = -0.25;
        
        dphi(0,2) = 0.;
        dphi(1,2) = 0.;
        dphi(2,2) = -0.25;
        
        dphi(0,3) = 0.;
        dphi(1,3) = 0.;
        dphi(2,3) = -0.25;
        
        dphi(0,4) = 0.;
        dphi(1,4) = 0.;
        dphi(2,4) = 3.;
        
        dphi(0,5) = 0.;
        dphi(1,5) = 0.5;
        dphi(2,5) = 0.5;
        
        dphi(0,6) = -0.5;
        dphi(1,6) = 0.;
        dphi(2,6) = 0.5;
        
        dphi(0,7) = 0.;
        dphi(1,7) = -0.5;
        dphi(2,7) = 0.5;
        
        dphi(0,8) = 0.5;
        dphi(1,8) = 0.;
        dphi(2,8) = 0.5;
        
        dphi(0,9) = -1.;
        dphi(1,9) = -1.;
        dphi(2,9) = -1.;
        
        dphi(0,10) = 1.;
        dphi(1,10) = -1.;
        dphi(2,10) = -1.;
        
        dphi(0,11) = 1.;
        dphi(1,11) = 1.;
        dphi(2,11) = -1.;
        
        dphi(0,12) = -1;
        dphi(1,12) = 1.;
        dphi(2,12) = -1.;
        
        return;
    }

    phi(0,0) = ((1 + eta + qsi)*(-1 + eta + zeta)*(-1 + qsi + zeta))/(4.*(-1 + zeta));
    phi(1,0) = -((1 + eta - qsi)*(1 + qsi - zeta)*(-1 + eta + zeta))/(4.*(-1 + zeta));
    phi(2,0) = -((-1 + eta + qsi)*(1 + eta - zeta)*(1 + qsi - zeta))/(4.*(-1 + zeta));
    phi(3,0) = ((-1 + eta - qsi)*(1 + eta - zeta)*(-1 + qsi + zeta))/(4.*(-1 + zeta));
    phi(4,0) = zeta*(-1 + 2*zeta);
    phi(5,0) = ((-1 + eta + zeta)*(-1 - qsi + zeta)*(-1 + qsi + zeta))/(2.*(-1 + zeta));
    phi(6,0) = ((-1 - eta + zeta)*(-1 + eta + zeta)*(-1 - qsi + zeta))/(2.*(-1 + zeta));
    phi(7,0) = ((qsi*qsi - (-1 + zeta)*(-1 + zeta))*(1 + eta - zeta))/(2.*(-1 + zeta));
    phi(8,0) = ((-(eta*eta) + (-1 + zeta)*(-1 + zeta))*(-1 + qsi + zeta))/(2.*(-1 + zeta));
    phi(9,0) = -((zeta*(-1 + eta + zeta)*(-1 + qsi + zeta))/(-1 + zeta));
    phi(10,0) = -((zeta*(-1 + eta + zeta)*(-1 - qsi + zeta))/(-1 + zeta));
    phi(11,0) = ((1 + eta - zeta)*(1 + qsi - zeta)*zeta)/(1 - zeta);
    phi(12,0) = -((zeta*(-1 - eta + zeta)*(-1 + qsi + zeta))/(-1 + zeta));


    dphi(0,0) = ((-1 + eta + zeta)*(eta + 2*qsi + zeta))/(4.*(-1 + zeta));
    dphi(1,0) = ((-1 + qsi + zeta)*(2*eta + qsi + zeta))/(4.*(-1 + zeta));
    dphi(2,0) = -(eta*eta*qsi + eta*(qsi + qsi*qsi - (-1 + zeta)*(-1 + zeta)) - (1 + qsi)*((-1 + zeta)*(-1 + zeta)))/(4.*((-1 + zeta)*(-1 + zeta)));

    dphi(0,1) = -((-1 + eta + zeta)*(eta - 2*qsi + zeta))/(4.*(-1 + zeta));
    dphi(1,1) = ((1 + qsi - zeta)*(-2*eta + qsi - zeta))/(4.*(-1 + zeta));
    dphi(2,1) = (eta*eta*qsi + eta*(qsi - qsi*qsi + (-1 + zeta)*(-1 + zeta)) - (-1 + qsi)*((-1 + zeta)*(-1 + zeta)))/(4.*((-1 + zeta)*(-1 + zeta)));

    dphi(0,2) = -((1 + eta - zeta)*(eta + 2*qsi - zeta))/(4.*(-1 + zeta));
    dphi(1,2) = -((1 + qsi - zeta)*(2*eta + qsi - zeta))/(4.*(-1 + zeta));
    dphi(2,2) = (eta*eta*qsi + eta*(-qsi + qsi*qsi - (-1 + zeta)*(-1 + zeta)) - (-1 + qsi)*((-1 + zeta)*(-1 + zeta)))/(4.*((-1 + zeta)*(-1 + zeta)));

    dphi(0,3) = ((1 + eta - zeta)*(eta - 2*qsi - zeta))/(4.*(-1 + zeta));
    dphi(1,3) = ((2*eta - qsi - zeta)*(-1 + qsi + zeta))/(4.*(-1 + zeta));
    dphi(2,3) = (-(eta*eta*qsi) + eta*(qsi + qsi*qsi - (-1 + zeta)*(-1 + zeta)) + (1 + qsi)*((-1 + zeta)*(-1 + zeta)))/(4.*((-1 + zeta)*(-1 + zeta)));

    dphi(0,4) = 0;
    dphi(1,4) = 0;
    dphi(2,4) = -1 + 4*zeta;

    dphi(0,5) = -((qsi*(-1 + eta + zeta))/(-1 + zeta));
    dphi(1,5) = (-(qsi*qsi) + (-1 + zeta)*(-1 + zeta))/(2.*(-1 + zeta));
    dphi(2,5) = (-2 + eta + (eta*(qsi*qsi))/((-1 + zeta)*(-1 + zeta)) + 2*zeta)/2.;

    dphi(0,6) = (eta*eta - (-1 + zeta)*(-1 + zeta))/(2.*(-1 + zeta));
    dphi(1,6) = (eta*(1 + qsi - zeta))/(-1 + zeta);
    dphi(2,6) = -1 + qsi*(-0.5 - (eta*eta)/(2.*((-1 + zeta)*(-1 + zeta)))) + zeta;

    dphi(0,7) = (qsi*(1 + eta - zeta))/(-1 + zeta);
    dphi(1,7) = (qsi*qsi - (-1 + zeta)*(-1 + zeta))/(2.*(-1 + zeta));
    dphi(2,7) = -1 + eta*(-0.5 - (qsi*qsi)/(2.*((-1 + zeta)*(-1 + zeta)))) + zeta;

    dphi(0,8) = (-(eta*eta) + (-1 + zeta)*(-1 + zeta))/(2.*(-1 + zeta));
    dphi(1,8) = -((eta*(-1 + qsi + zeta))/(-1 + zeta));
    dphi(2,8) = (-2 + qsi + (eta*eta*qsi)/((-1 + zeta)*(-1 + zeta)) + 2*zeta)/2.;

    dphi(0,9) = -((zeta*(-1 + eta + zeta))/(-1 + zeta));
    dphi(1,9) = -((zeta*(-1 + qsi + zeta))/(-1 + zeta));
    dphi(2,9) = 1 - qsi + eta*(-1 + qsi/((-1 + zeta)*(-1 + zeta))) - 2*zeta;

    dphi(0,10) = (zeta*(-1 + eta + zeta))/(-1 + zeta);
    dphi(1,10) = ((1 + qsi - zeta)*zeta)/(-1 + zeta);
    dphi(2,10) = 1 + qsi + eta*(-1 - qsi/((-1 + zeta)*(-1 + zeta))) - 2*zeta;

    dphi(0,11) = (zeta*(-1 - eta + zeta))/(-1 + zeta);
    dphi(1,11) = (zeta*(-1 - qsi + zeta))/(-1 + zeta);
    dphi(2,11) = 1 + eta + qsi + (eta*qsi)/((-1 + zeta)*(-1 + zeta)) - 2*zeta;

    dphi(0,12) = ((1 + eta - zeta)*zeta)/(-1 + zeta);
    dphi(1,12) = (zeta*(-1 + qsi + zeta))/(-1 + zeta);
    dphi(2,12) = 1 + eta - qsi - (eta*qsi)/((-1 + zeta)*(-1 + zeta)) - 2*zeta;
    
}


template<class T>
void TPZQuadraticPyramid::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<13,T> phi(NNodes,1);
    TPZFNMatrix<39,T> dphi(3,NNodes);
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
void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,3);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 13){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x13." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<13,T> phi(NNodes,1);
    TPZFNMatrix<39,T> dphi(3,NNodes);
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

void TPZQuadraticPyramid::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
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
    gmesh.CreateGeoElement(EPiramide, nodeindexes, matid, index);
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

//void TPZQuadraticPyramid::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
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
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 1.;
//            nodeCoord[2] = 0.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            nodeCoord[2] = 1.;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] = -1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (6):
//        {
//            nodeCoord[0] =  1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (7):
//        {
//            nodeCoord[0] =  0.;
//            nodeCoord[1] =  1.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (8):
//        {
//            nodeCoord[0] = -1.;
//            nodeCoord[1] =  0.;
//            nodeCoord[2] =  0.;
//            break;
//        }
//        case (9):
//        {
//            nodeCoord[0] = -0.5;
//            nodeCoord[1] = -0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (10):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] = -0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (11):
//        {
//            nodeCoord[0] =  0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        case (12):
//        {
//            nodeCoord[0] = -0.5;
//            nodeCoord[1] =  0.5;
//            nodeCoord[2] =  0.5;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticPyramid

int TPZQuadraticPyramid::ClassId() const{
    return Hash("TPZQuadraticPyramid") ^ TPZNodeRep<13,pztopology::TPZPyramid>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticPyramid>>;
/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<13,TPZQuadraticPyramid>;

namespace pzgeom {
    template void TPZQuadraticPyramid::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticPyramid::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
