/**
 * @file
 * @brief Contains the implementation of the TPZGeoLinear methods. 
 */

#include "TPZGeoLinear.h"
#include "pzquad.h"
#include "pzshapelinear.h"
#include "pzgeoel.h"
#include "tpzgeoelrefpattern.h"

#include "pzlog.h"

#ifdef LOG4CXX
static log4cxx::LoggerPtr logger(Logger::getLogger("pz.geom.pzgeolinear"));
#endif

using namespace pzshape;
using namespace std;

namespace pzgeom {
	
	void TPZGeoLinear::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
		REAL x = pt[0];
		phi(0,0) = (1-x)/2.;
		phi(1,0) = (1+x)/2.;
		dphi(0,0) = -0.5;
		dphi(0,1) = 0.5;
	}
	
	
	TPZGeoEl *TPZGeoLinear::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
		if(side==2) {
			TPZManVector<int> nodes(2);
			nodes[0] = orig->SideNodeIndex(side,0);
			nodes[1] = orig->SideNodeIndex(side,1);
			int index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
			//      TPZGeoEl1d *gel = new TPZGeoEl1d(nodes,bc,*orig->Mesh());
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,0)));
			TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,1)));
			TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else if(side==0 || side==1) {
			TPZManVector<int> nodeindexes(1);
			//      TPZGeoElPoint *gel;
			nodeindexes[0] = orig->SideNodeIndex(side,0);
			int index;
			TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
			//      gel = new TPZGeoElPoint(nodeindexes,bc,*(orig->Mesh()));
			TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
			return gel;
		}
		else PZError << "TPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
		return 0;
	}
	
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	TPZGeoEl *TPZGeoLinear::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
											 TPZVec<int>& nodeindexes,
											 int matid,
											 int& index)
	{
		return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
	}
	
    void TPZGeoLinear::Jacobian(TPZFMatrix<REAL> &coord,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
								TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) {
        
        //VERSAO FUNCIONAL
        jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
        int ic;
        REAL v1[3] = {0.};
        int nrow = coord.Rows();
        REAL mod1 = 0.;
        for(ic=0; ic<nrow; ic++) {
            v1[ic] = (coord(ic,1)-coord(ic,0))*0.5;
            mod1 += v1[ic]*v1[ic];
        }
        mod1 = sqrt(mod1);
        jacobian(0,0) = mod1;
        detjac = mod1;
        if(IsZero(detjac))
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << "Singular Jacobian " << detjac;
            LOGPZ_ERROR(logger, sout.str())
            detjac = ZeroTolerance();
        }
        jacinv(0,0) = 1./mod1;
        
        for(ic=0; ic<3; ic++) {
            axes(0,ic) = v1[ic]/mod1;
        }
    }
	
};
