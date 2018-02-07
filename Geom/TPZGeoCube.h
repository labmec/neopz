/**
 * @file
 * @brief Contains the TPZGeoCube class which implements the geometry of hexahedra element.
 */
// $Id: TPZGeoCube.h,v 1.12 2011-05-11 01:38:40 phil Exp $

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpzcube.h"
#include "pzfmatrix.h"

#include <string>

class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
    
    /**
     * @ingroup geometry
     * @brief Implements the geometry of hexahedra element. \ref geometry "Geometry"
     */
    class TPZGeoCube : public TPZNodeRep<8, pztopology::TPZCube> {
        
    public:
        /** @brief Number of corner nodes */
        enum {NNodes = 8};
        
        /** @brief Constructor with list of nodes */
        TPZGeoCube(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoCube::ClassId), TPZNodeRep<NNodes, pztopology::TPZCube>(nodeindexes)
        {
        }
        
        /** @brief Empty constructor */
        TPZGeoCube() : TPZRegisterClassId(&TPZGeoCube::ClassId),TPZNodeRep<NNodes, pztopology::TPZCube>()
        {
        }
        
        /** @brief Constructor with node map */
        TPZGeoCube(const TPZGeoCube &cp,
                   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoCube::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZCube>(cp,gl2lcNdMap)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoCube(const TPZGeoCube &cp) : TPZRegisterClassId(&TPZGeoCube::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoCube(const TPZGeoCube &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoCube::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
        {
        }
        
        static bool IsLinearMapping(int side)
        {
            return true;
        }
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "Hexahedron";}
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        
        /* @brief Compute x mapping from local parametric coordinates */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &x) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,x);
        }
        
        /** @brief Compute gradient of x mapping from local parametric coordinates */
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &loc, TPZFMatrix<T> &gradx) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
//            int nrow = coord.Rows();
//            int ncol = coord.Cols();
//            TPZFNMatrix<3*NNodes,T> nodes(nrow,ncol);
//            for(int i = 0; i < nrow; i++)
//            {
//                for(int j = 0; j < ncol; j++)
//                {
//                    nodes(i,j) = coord(i,j);
//                }
//            }
            
            GradX(coord,loc,gradx);
        }
        
        /** @brief Compute x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodecoordinates,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodecoordinates,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
        
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        
        public:
virtual int ClassId() const;
        void Read(TPZStream& buf, void* context);
        void Write(TPZStream& buf, int withclassid) const;

    public:
        /** @brief Creates a geometric element according to the type of the father element */
        static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                          TPZVec<int64_t>& nodeindexes,
                                          int matid,
                                          int64_t& index);
        
    };
    
    template<class T>
    inline void TPZGeoCube::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1] , zeta  = loc[2];
        
        T x[2],dx[2],y[2],dy[2],z[2],dz[2];
        x[0]    = (1.-qsi)/2.;
        x[1]    = (1.+qsi)/2.;
        dx[0]   = -0.5;
        dx[1]   = +0.5;
        y[0]    = (1.-eta)/2.;
        y[1]    = (1.+eta)/2.;
        dy[0]   = -0.5;
        dy[1]   = +0.5;
        z[0]    = (1.-zeta)/2.;
        z[1]    = (1.+zeta)/2.;
        dz[0]   = -0.5;
        dz[1]   = +0.5;
        
        phi(0,0) = x[0]*y[0]*z[0];
        phi(1,0) = x[1]*y[0]*z[0];
        phi(2,0) = x[1]*y[1]*z[0];
        phi(3,0) = x[0]*y[1]*z[0];
        phi(4,0) = x[0]*y[0]*z[1];
        phi(5,0) = x[1]*y[0]*z[1];
        phi(6,0) = x[1]*y[1]*z[1];
        phi(7,0) = x[0]*y[1]*z[1];
        dphi(0,0) = dx[0]*y[0]*z[0];
        dphi(1,0) = x[0]*dy[0]*z[0];
        dphi(2,0) = x[0]*y[0]*dz[0];
        dphi(0,1) = dx[1]*y[0]*z[0];
        dphi(1,1) = x[1]*dy[0]*z[0];
        dphi(2,1) = x[1]*y[0]*dz[0];
        dphi(0,2) = dx[1]*y[1]*z[0];
        dphi(1,2) = x[1]*dy[1]*z[0];
        dphi(2,2) = x[1]*y[1]*dz[0];
        dphi(0,3) = dx[0]*y[1]*z[0];
        dphi(1,3) = x[0]*dy[1]*z[0];
        dphi(2,3) = x[0]*y[1]*dz[0];
        dphi(0,4) = dx[0]*y[0]*z[1];
        dphi(1,4) = x[0]*dy[0]*z[1];
        dphi(2,4) = x[0]*y[0]*dz[1];
        dphi(0,5) = dx[1]*y[0]*z[1];
        dphi(1,5) = x[1]*dy[0]*z[1];
        dphi(2,5) = x[1]*y[0]*dz[1];
        dphi(0,6) = dx[1]*y[1]*z[1];
        dphi(1,6) = x[1]*dy[1]*z[1];
        dphi(2,6) = x[1]*y[1]*dz[1];
        dphi(0,7) = dx[0]*y[1]*z[1];
        dphi(1,7) = x[0]*dy[1]*z[1];
        dphi(2,7) = x[0]*y[1]*dz[1];
        
    }
    
    template<class T>
    inline void TPZGeoCube::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        TPZFNMatrix<8,T> phi(NNodes,1);
        TPZFNMatrix<24,T> dphi(3,NNodes);
        TShape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < NNodes; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
    }
    
    
    template<class T>
    inline void TPZGeoCube::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 || ncol  != 8){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x8." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<4,T> phi(NNodes,1);
        TPZFNMatrix<12,T> dphi(3,NNodes);
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
    
};

#endif
