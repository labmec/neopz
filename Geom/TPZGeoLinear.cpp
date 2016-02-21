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

    TPZGeoEl * TPZGeoLinear::CreateBCGeoEl(TPZGeoEl *orig, int side,int bc){
        if(side==2) {
            TPZManVector<long> nodes(2);
            nodes[0] = orig->SideNodeIndex(side,0);
            nodes[1] = orig->SideNodeIndex(side,1);
            long index;
            TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EOned,nodes,bc,index);
            TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,0)));
            TPZGeoElSide(gel,1).SetConnectivity(TPZGeoElSide(orig,TPZShapeLinear::ContainedSideLocId(side,1)));
            TPZGeoElSide(gel,2).SetConnectivity(TPZGeoElSide(orig,side));
            return gel;
        }
        else if(side==0 || side==1) {
            TPZManVector<long> nodeindexes(1);
            nodeindexes[0] = orig->SideNodeIndex(side,0);
            long index;
            TPZGeoEl *gel = orig->Mesh()->CreateGeoElement(EPoint,nodeindexes,bc,index);
            TPZGeoElSide(gel,0).SetConnectivity(TPZGeoElSide(orig,side));
            return gel;
        }
        else {
            PZError << "TPZGeoLinear::CreateBCGeoEl. Side = " << side << endl;
        }
        
        return 0;
    }
    
    TPZGeoEl * TPZGeoLinear::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                             TPZVec<long>& nodeindexes,
                                             int matid,
                                             long& index)
    {
        return CreateGeoElementPattern(mesh,type,nodeindexes,matid,index);
    }
    
    
    void TPZGeoLinear::Jacobian(const TPZFMatrix<REAL> &coord,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
                                TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) {
        
        jacobian.Resize(1,1); axes.Resize(1,3); jacinv.Resize(1,1);
        
        // Computing the Gradient of X
        TPZFNMatrix<3,REAL> gradx(3,1);
        gradx.Zero();
        
        int nrow = coord.Rows();
        REAL mod1 = 0.;
        for(int i = 0; i < nrow; i++) {
            gradx(i,0) = (coord.GetVal(i,1)-coord.GetVal(i,0))*0.5;
            mod1 += gradx(i,0)*gradx(i,0);
        }
        
        mod1 = sqrt(mod1);
        jacobian(0,0) = mod1;
        detjac = mod1;
        
        if(IsZero(detjac))
        {
            
#ifdef PZDEBUG
            std::stringstream sout;
            sout << "Singular Jacobian " << detjac;
            LOGPZ_ERROR(logger, sout.str())
#endif
            detjac = ZeroTolerance();
        }
        
        jacinv(0,0) = 1./detjac;
        for(int i=0; i < 3; i++) {
            axes(0,i) = gradx(i,0)/mod1;
        }
    }
    
}