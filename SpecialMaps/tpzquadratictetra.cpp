/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticTetra methods. 
 */
#include "tpzquadratictetra.h"
#include "pzshapetetra.h"
#include "tpzgeoelmapped.h"

#include "pzlog.h"

#include "tpzgeomid.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadratictetra");
#endif

#include "fad.h"

using namespace pzgeom;
using namespace pztopology;
using namespace pzshape;
using namespace std;

TPZQuadraticTetra::~TPZQuadraticTetra()
{
}

template<class T>
void TPZQuadraticTetra::TShape(const TPZVec<T> &par, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
{
	T qsi = par[0], eta = par[1], zeta = par[2];
	
    phi(0,0)  =  (zeta + eta + qsi -1.) * (2.*zeta + 2.*eta + 2.*qsi - 1.);
//    phi(0,0)  =  (qsi + eta + zeta -1.) * (2.*qsi + 2.*eta + 2.*zeta - 1.);
//    phi(0,0)  =  2.*qsi*qsi + 2.*eta*eta + 2.*zeta*zeta + 4.*qsi*eta + 4.*qsi*zeta + 4.*eta*zeta - 3.*qsi - 3.*eta - 3.*zeta + 1.;
	phi(1,0)  =  qsi  * (2.*qsi - 1.);
	phi(2,0)  =  eta  * (2.*eta - 1.);
	phi(3,0)  =  zeta * (2.*zeta - 1.);
	phi(4,0)  = -4.*qsi  * (qsi + eta + zeta -1.);
	phi(5,0)  =  4.*qsi  *  eta;
	phi(6,0)  = -4.*eta  * (qsi + eta + zeta -1.);
	phi(7,0)  = -4.*zeta * (qsi + eta + zeta -1.);
	phi(8,0)  =  4.*qsi  *  zeta;
	phi(9,0)  =  4.*eta  *  zeta;
	
	dphi(0,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(1,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(2,0) =  -3. + 4.*qsi + 4.*eta + 4.*zeta;
	dphi(0,1) =  -1. + 4.*qsi;
	dphi(1,1) =   0.;
	dphi(2,1) =   0.;
	dphi(0,2) =   0.;
	dphi(1,2) =  -1. + 4.*eta;
	dphi(2,2) =   0.;
	dphi(0,3) =   0.;
	dphi(1,3) =   0.;
	dphi(2,3) =  -1. + 4.*zeta;
	dphi(0,4) =  -4.*(2.*qsi + eta + zeta - 1.);
	dphi(1,4) =  -4.*qsi;
	dphi(2,4) =  -4.*qsi;
	dphi(0,5) =   4.*eta;
	dphi(1,5) =   4.*qsi;
	dphi(2,5) =   0.;
	dphi(0,6) =  -4.*eta;
	dphi(1,6) =  -4.*(qsi + eta + zeta - 1.)-4.*eta;
	dphi(2,6) =  -4.*eta;
	dphi(0,7) =  -4.*zeta;
	dphi(1,7) =  -4.*zeta;
	dphi(2,7) =  -4.*(qsi + eta + zeta - 1.)-4.*zeta;
	dphi(0,8) =   4.*zeta;
	dphi(1,8) =   0.;
	dphi(2,8) =   4.*qsi;
	dphi(0,9) =   0.;
	dphi(1,9) =   4.*zeta;
	dphi(2,9) =   4.*eta;
}


template<class T>
void TPZQuadraticTetra::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
    
    TPZFNMatrix<10,T> phi(NNodes,1);
    TPZFNMatrix<30,T> dphi(3,NNodes);
    TShape(loc,phi,dphi);
    int space = nodes.Rows();
    for (int i=0; i<3; i++) {
        x[i] = 0.;
    }
    for(int j = 0; j < NNodes; j++) {
        for(int i = 0; i < space; i++) {
            x[i] += phi(j,0)*nodes.GetVal(i,j);
        }
    }
}

template<class T>
void TPZQuadraticTetra::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,3);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 10){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x10." << std::endl;
        DebugStop();
    }
    
#endif
    TPZFNMatrix<10,T> phi(NNodes,1);
    TPZFNMatrix<30,T> dphi(3,NNodes);
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

void TPZQuadraticTetra::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
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
    gmesh.CreateGeoElement(ETetraedro, nodeindexes, matid, index);
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

int TPZQuadraticTetra::ClassId() const{
    return Hash("TPZQuadraticTetra") ^ TPZNodeRep<10,pztopology::TPZTetrahedron>::ClassId() << 1;
}

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticTetra>>;

template class TPZGeoElRefPattern<TPZQuadraticTetra>;
/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<10,TPZQuadraticTetra>;

namespace pzgeom {
    template void TPZQuadraticTetra::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticTetra::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticTetra::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticTetra::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
