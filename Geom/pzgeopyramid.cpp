/**
 * @file
 * @brief Contains the implementation of the TPZGeoPyramid methods. 
 */

#include "pzgeopyramid.h"
#include "pzfmatrix.h"
#include "pzgeoel.h"
#include "pzquad.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.geom.pzgeopyramid"));
#endif

#include <cmath>

using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	const double tol = pzgeom_TPZNodeRep_tol;
	
//    void TPZGeoPyramid::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
//        if(fabs(pt[0])<1.e-10 && fabs(pt[1])<1.e-10 && pt[2]==1.) {
//            //para testes com transforma�es geometricas-->>Que  o que faz o RefPattern!!
//            //(0,0,1) nunca �um ponto de integra�o
//            phi(0,0)  = 0.;
//            phi(1,0)  = 0.;
//            phi(2,0)  = 0.;
//            phi(3,0)  = 0.;
//            phi(4,0)  = 1.;
//            dphi(0,0)  = -0.25;
//            dphi(1,0)  = -0.25;
//            dphi(2,0)  = -0.25;
//            dphi(0,1)  = 0.25;
//            dphi(1,1)  = -0.25;
//            dphi(2,1)  = -0.25;
//            dphi(0,2)  = 0.25;
//            dphi(1,2)  = 0.25;
//            dphi(2,2)  = -0.25;
//            dphi(0,3)  = -0.25;
//            dphi(1,3)  = 0.25;
//            dphi(2,3)  = -0.25;
//            dphi(0,4)  = 0;
//            dphi(1,4)  = 0;
//            dphi(2,4)  = 1.;
//            
//            
//            
//            return;
//        }
//        
//        REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
//        REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
//        REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
//        REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
//        REAL lmez = (1.-pt[2]);
//        phi(0,0)  = T0xz*T0yz*lmez;
//        phi(1,0)  = T1xz*T0yz*lmez;
//        phi(2,0)  = T1xz*T1yz*lmez;
//        phi(3,0)  = T0xz*T1yz*lmez;
//        phi(4,0)  = pt[2];
//        REAL lmexmez = 1.-pt[0]-pt[2];
//        REAL lmeymez = 1.-pt[1]-pt[2];
//        REAL lmaxmez = 1.+pt[0]-pt[2];
//        REAL lmaymez = 1.+pt[1]-pt[2];
//        dphi(0,0) = -.25*lmeymez / lmez;
//        dphi(1,0) = -.25*lmexmez / lmez;
//        dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
//        
//        dphi(0,1) =  .25*lmeymez / lmez;
//        dphi(1,1) = -.25*lmaxmez / lmez;
//        dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
//        
//        dphi(0,2) =  .25*lmaymez / lmez;
//        dphi(1,2) =  .25*lmaxmez / lmez;
//        dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
//        
//        dphi(0,3) = -.25*lmaymez / lmez;
//        dphi(1,3) =  .25*lmexmez / lmez;
//        dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
//        
//        dphi(0,4) =  0.0;
//        dphi(1,4) =  0.0;
//        dphi(2,4) =  1.0;
//    }
    
	
	TPZGeoEl *TPZGeoPyramid::CreateBCGeoEl(TPZGeoEl *orig,int side,int bc) {
		if(side<0 || side>18) {
			cout << "TPZGeoPyramid::CreateBCGeoEl Bad parameter side = " 
			<< side << "not implemented\n";
			return 0;
		}
		
		if(side==18) {
			cout << "TPZGeoPyramid::CreateBCCompEl with side = 18 not implemented\n";
			return 0;
		}
		
		if(side<5) {
			TPZManVector<int64_t> nodeindexes(1);
			//		TPZGeoElPoint *gel;
			nodeindexes[0] = orig->NodeIndex(side);
			int64_t index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//		gel = new TPZGeoElPoint(nodeindexes,bc,*orig->Mesh());
			TPZGeoElSide origside(orig,side);
			TPZGeoElSide(gel,0).SetConnectivity(origside);
			return gel;
		} 
		else if (side > 4 && side < 13) {//side =5 a 12 : lados
			TPZManVector<int64_t> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			int64_t index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			//		TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if (side > 12) {//side = 13 a 17 : faces
			TPZManVector<int64_t> nodes(4);//4o = -1 para face triangular
			int iside;
			for (iside=0;iside<4;iside++){
				nodes[iside] = orig->SideNodeIndex(side,iside);
			}
			if(side==13) {
				int64_t index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EQuadrilateral,nodes,bc,index);
				//      		gelq = new TPZGeoElQ2d(nodes,bc,*orig->Mesh());
				for (iside=0; iside<8; iside++){
					TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gel,8).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			} 
			else {
				nodes.Resize(3);
				int64_t index;
				TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(ETriangle,nodes,bc,index);
				//			gelt = new TPZGeoElT2d(nodes,bc,*orig->Mesh());
				for (iside=0; iside<6; iside++){
					TPZGeoElSide(gel,iside).SetConnectivity(TPZGeoElSide(orig,TPZShapePiram::ContainedSideLocId(side,iside)));
				}
				TPZGeoElSide(gel,6).SetConnectivity(TPZGeoElSide(orig,side));
				return gel;
			}
		} 
		else 
			PZError << "TPZGeoPyramid::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	void TPZGeoPyramid::FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
	{
		ChangedPoint.Resize(OriginalPoint.NElements(),0.);
		ChangedPoint = OriginalPoint;
		
		switch(side)
		{
			case 5:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 6:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 7:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 8:
			{
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 9:
			{
				if( ChangedPoint[0] == 1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == 1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 10:
			{
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == 1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 11:
			{
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] ==  1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] =  1. - tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 12:
			{
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] ==  1. )
				{
					ChangedPoint[0] = 1. - tol;
					ChangedPoint[1] = 1. - tol;
				}
				if( ChangedPoint[0] ==  1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] =  1. - tol;
					ChangedPoint[1] = -1. + tol;
				}
				if( ChangedPoint[0] == -1. && ChangedPoint[1] == -1. )
				{
					ChangedPoint[0] = -1. + tol;
					ChangedPoint[1] = -1. + tol;
				}
				if(OriginalPoint[2] ==  1.)
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 13:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				break;
			}
				
			case 14:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  1.) )
				{
					ChangedPoint[1] = 1. - tol;
				}
				break;
			}
				
			case 15:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] ==  -1.) )
				{
					ChangedPoint[0] = -1. + tol;
				}
				break;
			}
				
			case 16:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] == -1.) || (ChangedPoint[0] ==  -1. && ChangedPoint[0] == -1.) )
				{
					ChangedPoint[1] = -1. + tol;
				}
				break;
			}
				
			case 17:
			{
				if( ChangedPoint[2] ==  1. )
				{
					ChangedPoint[2] = 1. - tol;
				}
				if( (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  1.) || (ChangedPoint[0] ==  1. && ChangedPoint[0] ==  -1.) )
				{
					ChangedPoint[0] = 1. - tol;
				}
				break;
			}
		}
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoPyramid::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											  TPZVec<int64_t>& nodeindexes,
											  int matid,
											  int64_t& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}

    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoPyramid::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<int64_t,3> nodeindexes(8);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i=0; i<NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j=0; j<3; j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        int64_t index;
        CreateGeoElement(gmesh, EPiramide, nodeindexes, matid, index);
    }
    
    int TPZGeoPyramid::ClassId() const{
        return Hash("TPZGeoPyramid") ^ TPZNodeRep<5, pztopology::TPZPyramid>::ClassId() << 1;
    }
    
    void TPZGeoPyramid::Read(TPZStream& buf, void* context) {
        TPZNodeRep<5, pztopology::TPZPyramid>::Read(buf,context);
    }

    void TPZGeoPyramid::Write(TPZStream& buf, int withclassid) const {
        TPZNodeRep<5, pztopology::TPZPyramid>::Write(buf,withclassid);
    }

    
};
