/**
 * @file
 * @brief Contains the implementation of the TPZGeoQuad methods. 
 */

#include "pzgeoquad.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeoquad"));
#endif

using namespace std;

namespace pzgeom {
    const REAL tol = pzgeom_TPZNodeRep_tol;

    template<class T>
    void TPZGeoQuad::CalcSideInfluence(const int &side, const TPZVec<T> &xi, T &correctionFactor,
                                       TPZVec<T> &corrFactorDxi){
#ifdef PZDEBUG
        std::ostringstream sout;
        if(side < NNodes || side >= NSides){
            sout<<"The side\t"<<side<<"is invalid. Aborting..."<<std::endl;

            PZError<<std::endl<<sout.str()<<std::endl;
            DebugStop();
        }

        if(!IsInParametricDomain(xi,tol)){
            sout<<"The method CalcSideInfluence expects the point xi to correspond to coordinates of a point";
            sout<<" inside the parametric domain. Aborting...";
            PZError<<std::endl<<sout.str()<<std::endl;
            #ifdef LOG4CXX
            LOGPZ_FATAL(logger,sout.str().c_str());
            #endif
            DebugStop();
        }
#endif
        TPZFNMatrix<4,T> phi(NNodes,1);
        TPZFNMatrix<8,T> dphi(Dimension,NNodes);
        TPZGeoQuad::TShape(xi,phi,dphi);
        corrFactorDxi.Resize(TPZGeoQuad::Dimension, (T) 0);
        int i = -1;
        switch(side){
            case 0:
            case 1:
            case 2:
            case 3:
                correctionFactor = 0;
                return;
            case 4:
                i = 0;
                break;
            case 5:
                i = 1;
                break;
            case 6:
                i = 2;
                break;
            case 7:
                i = 3;
                break;
            case 8:
                correctionFactor = 1;
                return;
        }
        correctionFactor = phi(i,0) + phi((i+1)%NNodes,0);
        corrFactorDxi[0] = dphi(0,i) + dphi(0,(i+1)%NNodes);
        corrFactorDxi[1] = dphi(1,i) + dphi(1,(i+1)%NNodes);
    }
	
	void TPZGeoQuad::VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result){
		if(v1.NElements()!=3||v2.NElements()!=3)
		{
			cout << " o tamanho do vetores eh diferente de 3"<< endl;
		}
		REAL x1=v1[0], y1=v1[1],z1=v1[2];
		REAL x2=v2[0], y2=v2[1],z2=v2[2];
		result.Resize(v1.NElements());
		result[0]=y1*z2-z1*y2;
		result[1]=z1*x2-x1*z2;	
		result[2]=x1*y2-y1*x2;	
	}
	
	void TPZGeoQuad::ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result){
		TPZVec<REAL> v1(3);
		TPZVec<REAL> v2(3);
		TPZVec<REAL> normal(3);
		v1[0]=p1[0]-p2[0];
		v1[1]=p1[1]-p2[1];
		v1[2]=p1[2]-p2[2];
		v2[0]=p2[0]-p3[0];
		v2[1]=p2[1]-p3[1];
		v2[2]=p2[2]-p3[2];
		VectorialProduct(v1,v2,normal);
		VectorialProduct(v1,normal,result);	
	}
	
	
	TPZGeoEl *TPZGeoQuad::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		if(side==8) {//8
			TPZManVector<int64_t> nodes(4);
			int i;
			for (i=0;i<4;i++) {
				nodes[i] = orig->SideNodeIndex(side,i);
			}
			int64_t index;
			TPZGeoEl *gel = orig->CreateGeoElement(EQuadrilateral,nodes,bc,index);
			int iside;
			for (iside = 0; iside <8; iside++){
				TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,iside)));
			}
			TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side>-1 && side<4) {//side = 0,1,2,3
			TPZManVector<int64_t> nodeindexes(1);
			nodeindexes[0] = orig->SideNodeIndex(side,0);
			int64_t index;
			TPZGeoEl *gel = orig->CreateGeoElement(EPoint,nodeindexes,bc,index);
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side>3 && side<8) {
			TPZManVector<int64_t> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			int64_t index;
			TPZGeoEl *gel = orig->CreateGeoElement(EOned,nodes,bc,index);
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,pztopology::TPZQuadrilateral::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else PZError << "TPZGeoQuad::CreateBCCompEl has no bc.\n";
		return 0;
	}
	
    void TPZGeoQuad::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<int64_t,4> nodeindexes(4);
        TPZManVector<REAL,3> co(lowercorner);
        for (int i=0; i<3; i++) {
            co[i] += 0.2*size[i];   
        }
        
        nodeindexes[0] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[0]].Initialize(co, gmesh);
        co[0] += 0.6*size[0];
        nodeindexes[1] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[1]].Initialize(co, gmesh);
        co[1] += 0.6*size[0];
        co[0] += 0.1*size[0];
        co[2] += 0.3*size[0];
        nodeindexes[2] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[2]].Initialize(co, gmesh);
        for (int i=0; i<3; i++) co[i] = lowercorner[i]+0.2*size[i];
        co[1] += 0.4*size[1];
        co[2] -= 0.2*size[2];
        nodeindexes[3] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[3]].Initialize(co, gmesh);
        int64_t index;
        CreateGeoElement(gmesh, EQuadrilateral, nodeindexes, matid, index);
    }

    
	TPZGeoEl *TPZGeoQuad::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										   TPZVec<int64_t>& nodeindexes,
										   int matid,
										   int64_t& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
        
        int TPZGeoQuad::ClassId() const{
            return Hash("TPZGeoQuad") ^ TPZNodeRep<4, pztopology::TPZQuadrilateral>::ClassId() << 1;
    }

    void TPZGeoQuad::Read(TPZStream& buf, void* context) {
        TPZNodeRep<4, pztopology::TPZQuadrilateral>::Read(buf, context);
    }

    void TPZGeoQuad::Write(TPZStream& buf, int withclassid) const {
        TPZNodeRep<4, pztopology::TPZQuadrilateral>::Write(buf, withclassid);
    }

    template void TPZGeoQuad::CalcSideInfluence<REAL>(const int &, const TPZVec<REAL> &, REAL &, TPZVec<REAL> &);

};

#ifdef _AUTODIFF
template<class T=REAL>
class Fad;

template void pzgeom::TPZGeoQuad::CalcSideInfluence<Fad<REAL>>(const int &, const TPZVec<Fad<REAL>> &, Fad<REAL> &,
        TPZVec<Fad<REAL>> &);
#endif
