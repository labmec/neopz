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

#include <string>

template<class TVar>
class TPZFMatrix;
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
		TPZGeoCube(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZCube>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoCube() : TPZNodeRep<NNodes, pztopology::TPZCube>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoCube(const TPZGeoCube &cp,
				   std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoCube(const TPZGeoCube &cp) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoCube(const TPZGeoCube &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
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
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &x) const
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
            int nrow = coord.Rows();
            int ncol = coord.Cols();
            TPZFMatrix<T> nodes(nrow,ncol);
            for(int i = 0; i < nrow; i++)
            {
                for(int j = 0; j < ncol; j++)
                {
                    nodes(i,j) = coord(i,j);
                }
            }
            
            GradX(nodes,loc,gradx);
        }
        
        /** @brief Compute x mapping from element nodes and local parametric coordinates */
        static void X(const TPZFMatrix<REAL> &nodecoordinates,TPZVec<REAL> &loc,TPZVec<REAL> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<T> &nodecoordinates,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
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
    
    inline void TPZGeoCube::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &x){
        
        TPZFNMatrix<8,REAL> phi(8,1);
        TPZFNMatrix<24,REAL> dphi(3,8);
        Shape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 8; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
    }
    
    
    template<class T>
    inline void TPZGeoCube::GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 && ncol  != 8){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x8." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<4,T> phi(8,1);
        TPZFNMatrix<12,T> dphi(3,8);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 8; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
                gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
                gradx(j,2) += nodes.GetVal(j,i)*dphi(2,i);
            }
        }
        
    }
    
//    void TPZGeoCube::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> & loc,TPZVec<REAL> &result){
//        
//        int nrow = nodes.Rows();
//        int ncol = nodes.Cols();
//        if(nrow != 3 || ncol != 8){//8x3 n� por linhas
//            cout << "TPZGeoCube::X nodes matrix error size\n";
//        }
//        int i,j;
//        TPZFNMatrix<12> phi(8,1);
//        TPZFNMatrix<30> dphi(3,8);
//        Shape(loc,phi,dphi);
//        for(j=0;j<3;j++) {
//            result[j] = 0.0;
//            for(i=0;i<8;i++) result[j] += nodes.GetVal(j,i)*phi(i,0);
//        }
//    }
    
//    void TPZGeoCube::Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
//        
//        REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
//        x[0] = (1.-pt[0])/2.;
//        x[1] = (1.+pt[0])/2.;
//        dx[0] = -0.5;
//        dx[1] = 0.5;
//        y[0] = (1.-pt[1])/2.;
//        y[1] = (1.+pt[1])/2.;
//        dy[0] = -0.5;
//        dy[1] = 0.5;
//        z[0] = (1.-pt[2])/2.;
//        z[1] = (1.+pt[2])/2.;
//        dz[0] = -0.5;
//        dz[1] = 0.5;
//        
//        phi(0,0) = x[0]*y[0]*z[0];
//        phi(1,0) = x[1]*y[0]*z[0];
//        phi(2,0) = x[1]*y[1]*z[0];
//        phi(3,0) = x[0]*y[1]*z[0];
//        phi(4,0) = x[0]*y[0]*z[1];
//        phi(5,0) = x[1]*y[0]*z[1];
//        phi(6,0) = x[1]*y[1]*z[1];
//        phi(7,0) = x[0]*y[1]*z[1];
//        dphi(0,0) = dx[0]*y[0]*z[0];
//        dphi(1,0) = x[0]*dy[0]*z[0];
//        dphi(2,0) = x[0]*y[0]*dz[0];
//        dphi(0,1) = dx[1]*y[0]*z[0];
//        dphi(1,1) = x[1]*dy[0]*z[0];
//        dphi(2,1) = x[1]*y[0]*dz[0];
//        dphi(0,2) = dx[1]*y[1]*z[0];
//        dphi(1,2) = x[1]*dy[1]*z[0];
//        dphi(2,2) = x[1]*y[1]*dz[0];
//        dphi(0,3) = dx[0]*y[1]*z[0];
//        dphi(1,3) = x[0]*dy[1]*z[0];
//        dphi(2,3) = x[0]*y[1]*dz[0];
//        dphi(0,4) = dx[0]*y[0]*z[1];
//        dphi(1,4) = x[0]*dy[0]*z[1];
//        dphi(2,4) = x[0]*y[0]*dz[1];
//        dphi(0,5) = dx[1]*y[0]*z[1];
//        dphi(1,5) = x[1]*dy[0]*z[1];
//        dphi(2,5) = x[1]*y[0]*dz[1];
//        dphi(0,6) = dx[1]*y[1]*z[1];
//        dphi(1,6) = x[1]*dy[1]*z[1];
//        dphi(2,6) = x[1]*y[1]*dz[1];
//        dphi(0,7) = dx[0]*y[1]*z[1];
//        dphi(1,7) = x[0]*dy[1]*z[1];
//        dphi(2,7) = x[0]*y[1]*dz[1];
//    }
	
};

#endif
