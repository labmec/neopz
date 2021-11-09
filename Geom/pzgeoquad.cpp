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

#ifdef PZ_LOG
static TPZLogger logger("pz.geom.pzgeoquad");
#endif

using namespace std;

namespace pzgeom {
    const REAL tol = pzgeom_TPZNodeRep_tol;
	
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
	
	
	
    void TPZGeoQuad::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<int64_t,4> nodeindexes(4);
        TPZManVector<REAL,3> co(lowercorner);
        for (int i=0; i<3; i++) {
            co[i] += 0.2*size[i];
        }

        nodeindexes[0] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[0]].Initialize(co, gmesh);
        co[0] += 0.6 * size[0];
        nodeindexes[1] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[1]].Initialize(co, gmesh);
        co[1] += 0.6 * size[0];
        co[0] += 0.1 * size[0];
        co[2] += 0.3 * size[0];
        nodeindexes[2] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[2]].Initialize(co, gmesh);
        for (int i = 0; i < 3; i++) co[i] = lowercorner[i] + 0.2 * size[i];
        co[1] += 0.4 * size[1];
        co[2] -= 0.2 * size[2];
        nodeindexes[3] = gmesh.NodeVec().AllocateNewElement();
        gmesh.NodeVec()[nodeindexes[3]].Initialize(co, gmesh);
        int64_t index;
        gmesh.CreateGeoElement(EQuadrilateral, nodeindexes, matid, index);
    }

        
        int TPZGeoQuad::ClassId() const{
            return Hash("TPZGeoQuad") ^ TPZNodeRep<4, pztopology::TPZQuadrilateral>::ClassId() << 1;
    }

    void TPZGeoQuad::Read(TPZStream &buf, void *context) {
        TPZNodeRep<4, pztopology::TPZQuadrilateral>::Read(buf, context);
    }

    void TPZGeoQuad::Write(TPZStream &buf, int withclassid) const {
        TPZNodeRep<4, pztopology::TPZQuadrilateral>::Write(buf, withclassid);
    }

};
