/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticTrig methods. 
 */
#include "tpzquadratictrig.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzgeotriangle.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "tpztriangle.h"
#include "pznoderep.h"
#include "pzshapetriang.h"
#include "tpzgeoelmapped.h"

#include <iostream>
#include <iostream>
#include <cmath>

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"

#include "tpzgeomid.h"

#include "fad.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;

/**
 / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
 / LabMeC - FEC - UNICAMP
 / 2007
 */

template<class T>
void TPZQuadraticTrig::TShape(const TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi)
{
	T qsi = param[0], eta = param[1];
    
	phi(0,0) = (qsi+eta-1.)*(2.*qsi+2.*eta-1.);
	phi(1,0) = qsi*(2.*qsi-1.);
	phi(2,0) = eta*(2.*eta-1.);
	phi(3,0) = -4.*qsi*(qsi+eta-1.);
	phi(4,0) =  4.*qsi*eta;
	phi(5,0) = -4.*eta*(qsi+eta-1.);
	dphi(0,0) = 1.-4.*(1.-qsi-eta);
    dphi(0,1) = -1.+4.*qsi;
    dphi(0,2) = 0.;
    dphi(0,3) = 4.*(1.-qsi-eta)-4.*qsi;
    dphi(0,4) = 4.*eta;
    dphi(0,5) = -4.*eta;
	dphi(1,0) = 1.-4.*(1.-qsi-eta);
    dphi(1,1) = 0.;
    dphi(1,2) = -1.+4.*eta;
    dphi(1,3) = -4.*qsi;
    dphi(1,4) = 4.*qsi;
    dphi(1,5) = 4.*(1.-qsi-eta)-4.*eta;
}

template<class T>
void TPZQuadraticTrig::X(const TPZFMatrix<REAL> &coord, TPZVec<T>& par, TPZVec< T >& result)
{
	TPZManVector<T,3> parMap(2);
	REAL spacephi[6],spacedphi[12];
	TPZFNMatrix<6,T> phi(NNodes,1);
	TPZFNMatrix<12,T> dphi(2,NNodes);
	TShape(par,phi,dphi);
	for(int i = 0; i < 3; i++)
	{
		result[i] = 0.0;
		for(int j = 0; j < 6; j++) result[i] += phi(j,0)*coord.GetVal(i,j);
	}
}


template<class T>
void TPZQuadraticTrig::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,2);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 6){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x6." << std::endl;
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

void TPZQuadraticTrig::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZManVector<REAL,3> co(3),shift(3),scale(3);
    TPZManVector<int64_t,3> nodeindexes(3);
    for (int i=0; i<3; i++) {
        scale[i] = size[i]/3.;
        shift[i] = 1./2.+lowercorner[i];
    }
    
    for (int i=0; i<NCornerNodes; i++) {
        ParametricDomainNodeCoord(i, co);
        co.Resize(3, 0.);
        for (int j=0; j<co.size(); j++) {
            co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
    }
    int64_t index;
    gmesh.CreateGeoElement(ETriangle, nodeindexes, matid, index);
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


//void TPZQuadraticTrig::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
//{
//    if(node > this->NNodes)
//    {
//        DebugStop();
//    }
//    nodeCoord.Resize(Dimension, 0.);
//    switch (node) {
//        case (0):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 0.;
//            break;
//        }
//        case (1):
//        {
//            nodeCoord[0] = 1.;
//            nodeCoord[1] = 0.;
//            break;
//        }
//        case (2):
//        {
//            nodeCoord[0] = 0.;
//            nodeCoord[1] = 1.;
//            break;
//        }
//        case (3):
//        {
//            nodeCoord[0] = 0.5;
//            nodeCoord[1] = 0.0;
//            break;
//        }
//        case (4):
//        {
//            nodeCoord[0] = 0.5;
//            nodeCoord[1] = 0.5;
//            break;
//        }
//        case (5):
//        {
//            nodeCoord[0] = 0.0;
//            nodeCoord[1] = 0.5;
//            break;
//        }
//        default:
//        {
//            DebugStop();
//            break;
//        }
//    }
//}

///CreateGeoElement -> TPZQuadraticTrig

int TPZQuadraticTrig::ClassId() const{
    return Hash("TPZQuadraticTrig") ^ TPZNodeRep<6,pztopology::TPZTriangle>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTrig>>;

/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<6,TPZQuadraticTrig>;

namespace pzgeom {
    template void TPZQuadraticTrig::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticTrig::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticTrig::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticTrig::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
