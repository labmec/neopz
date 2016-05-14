/**
 * @file
 * @brief Contains the implementation of the TPZGeoPoint methods. 
 */

#include "pzgeopoint.h"
#include "pzquad.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"

using namespace std;

namespace pzgeom {
	
	
	TPZGeoEl *TPZGeoPoint::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
		if(side==0) {
			TPZManVector<long> nodeindexes(1);
			nodeindexes[0] = orig->NodeIndex(0);
			long index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,0));
			return gel;
		}
		else PZError << "TPZGeoPoint::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	/** Creates a geometric element according to the type of the father element */
	TPZGeoEl *TPZGeoPoint::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											TPZVec<long>& nodeindexes,
											int matid,
											long& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoPoint::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<long,1> nodeindexes(1);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i = 0; i < NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j=0; j< co.size(); j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        long index;
        CreateGeoElement(gmesh, EPoint, nodeindexes, matid, index);
    }
    

};
