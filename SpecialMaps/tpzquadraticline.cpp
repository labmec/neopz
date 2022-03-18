/**
 * @file
 * @brief Contains the implementation of the TPZQuadraticLine methods. 
 */
#include "tpzquadraticline.h"
#include "pzshapequad.h"
#include "tpzgeoblend.h"
#include "tpzgeoelmapped.h"

#include "pzgeoelrefless.h"
#include "tpzgeoelrefpattern.h"
#include "pznoderep.h"

#include "tpzgeomid.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.specialmaps.quadraticline");
#endif

#include "fad.h"

using namespace pzshape;
using namespace pzgeom;
using namespace pztopology;


template<class T>
void TPZQuadraticLine::TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
    
    T qsi = loc[0];
    
    phi(0,0)  = -qsi*(1.-qsi)/2.;
    phi(1,0)  = +qsi*(1.+qsi)/2.;
    phi(2,0)  = (1.-qsi)*(1.+qsi);
    
    dphi(0,0) = qsi-0.5;
    dphi(0,1) = qsi+0.5;
    dphi(0,2) = -2.*qsi;
}

template<class T>
void TPZQuadraticLine::X(const TPZFMatrix<REAL> & coord, TPZVec<T> & loc,TPZVec<T> &result) {
    TPZFNMatrix<9,T> phi(NNodes,1);
    TPZFNMatrix<16,T> dphi(1,NNodes);
    TShape(loc,phi,dphi);
    
    for(int i = 0; i < 3; i++){
        result[i] = 0.0;
        for(int j = 0; j < NNodes; j++) result[i] += phi(j,0)*coord.GetVal(i,j);
    }
}

template<class T>
void TPZQuadraticLine::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
    
    gradx.Resize(3,1);
    gradx.Zero();
    int nrow = nodes.Rows();
    int ncol = nodes.Cols();
#ifdef PZDEBUG
    if(nrow != 3 || ncol  != 3){
        std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
        std::cout << "nodes matrix must be 3x3." << std::endl;
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

        }
        
    }
    
}



#include "tpzgeoelmapped.h"
/// create an example element based on the topology
/* @param gmesh mesh in which the element should be inserted
 @param matid material id of the element
 @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
 @param size (in) size of space where the element should be created
 */
#include "tpzchangeel.h"

void TPZQuadraticLine::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
{
    TPZManVector<REAL,3> co(3),shift(3),scale(3);
    TPZManVector<int64_t,3> nodeindexes(2);
    for (int i=0; i<3; i++) {
        scale[i] = size[i]/3.;
        shift[i] = 1./2.+lowercorner[i];
    }
    
    for (int i=0; i<NCornerNodes; i++) {
        ParametricDomainNodeCoord(i, co);
        int j;
        for (j=0; j<co.size(); j++) {
            co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        co.Resize(3);
        for (; j<3; j++) {
            co[j] = shift[j]+(rand()*0.2/RAND_MAX)-0.1;
        }
        nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
    }
    int64_t index;
    gmesh.CreateGeoElement(EOned, nodeindexes, matid, index);
    TPZGeoEl *gel = gmesh.Element(index);
    int nsides = gel->NSides();
    for (int is = 0; is<nsides; is++) {
        gel->SetSideDefined(is);
    }
    
    gel = TPZChangeEl::ChangeToQuadratic(&gmesh,index);
    TPZGeoNode *nodeptr = gel->NodePtr(2);
    for (int i=0; i<3; i++) {
        co[i] = lowercorner[i]+size[i]/2.;
    }
    gel->NodePtr(2)->SetCoord(co);
}


///CreateGeoElement -> TPZQuadraticLine

int TPZQuadraticLine::ClassId() const{
    return Hash("TPZQuadraticLine") ^ TPZNodeRep<3,pztopology::TPZLine>::ClassId() << 1;
}

template class TPZRestoreClass< TPZGeoElRefPattern<TPZQuadraticLine>>;
/*@orlandini : I REALLY dont know why is this here, so I have commented the following lines.
If it breaks something, I am sorry.*/
//template class pzgeom::TPZNodeRep<3,TPZQuadraticLine>;

namespace pzgeom {
    template void TPZQuadraticLine::X(const TPZFMatrix<REAL>&, TPZVec<REAL>&, TPZVec<REAL>&);
    template void TPZQuadraticLine::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx);

    template void TPZQuadraticLine::X(const TPZFMatrix<REAL>&, TPZVec<Fad<REAL> >&, TPZVec<Fad<REAL> >&);
    template void TPZQuadraticLine::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<Fad<REAL> > &loc, TPZFMatrix<Fad<REAL> > &gradx);

}
