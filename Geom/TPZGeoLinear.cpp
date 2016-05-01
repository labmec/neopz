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
        GradX(coord, param, gradx);
        
        int nrow = gradx.Rows();
        REAL module = 0.;
        for(int i = 0; i < nrow; i++) {
            module += gradx(i,0)*gradx(i,0);
        }
        
        module = sqrt(module);
        jacobian(0,0) = module;
        detjac = module;
        
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
            axes(0,i) = gradx(i,0)/module;
        }
    }
    
    /// create an example element based on the topology
    /* @param gmesh mesh in which the element should be inserted
     @param matid material id of the element
     @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
     @param size (in) size of space where the element should be created
     */
    void TPZGeoLinear::InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size)
    {
        TPZManVector<REAL,3> co(3),shift(3),scale(3);
        TPZManVector<long,3> nodeindexes(2);
        for (int i=0; i<3; i++) {
            scale[i] = size[i]/3.;
            shift[i] = 1./2.+lowercorner[i];
        }
        
        for (int i=0; i<NCornerNodes; i++) {
            ParametricDomainNodeCoord(i, co);
            for (int j=0; j<co.size(); j++) {
                co[j] = shift[j]+scale[j]*co[j]+(rand()*0.2/RAND_MAX)-0.1;
            }
            nodeindexes[i] = gmesh.NodeVec().AllocateNewElement();
            gmesh.NodeVec()[nodeindexes[i]].Initialize(co, gmesh);
        }
        long index;
        CreateGeoElement(gmesh, EOned, nodeindexes, matid, index);
        lowercorner = co;
    }
    

}