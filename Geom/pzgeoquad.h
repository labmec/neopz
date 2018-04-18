/**
 * @file TPZGeoQuad
 * @brief Implements the geometry of a quadrilateral element.
 */

#ifndef TPZGEOQUADH
#define TPZGEOQUADH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"
#include "pzfmatrix.h"

#include <string>

class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a quadrilateral element. \ref geometry "Geometry"
	 */
	class TPZGeoQuad  : public TPZNodeRep<4, pztopology::TPZQuadrilateral>
	{
	public:
		/** @brief Number of corner nodes */
		enum {NNodes = 4};

		/** @brief Constructor with list of nodes */
		TPZGeoQuad(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoQuad::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoQuad() : TPZRegisterClassId(&TPZGeoQuad::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoQuad(const TPZGeoQuad &cp,
				   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoQuad::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoQuad(const TPZGeoQuad &cp) : TPZRegisterClassId(&TPZGeoQuad::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoQuad(const TPZGeoQuad &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoQuad::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
        static bool IsLinearMapping(int side)
        {
            return true;
        }
        
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Quad";}
        
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
//            TPZFMatrix<T> nodes(nrow,ncol);
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
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
		
        /**
         * @brief Method which creates a geometric boundary condition
         * element based on the current geometric element, \n
         * a side and a boundary condition number
         */
        static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
        
        
		/** @brief Implementation of normal vector to Hdiv space*/
		/** 
		 construct the normal vector for element Hdiv
		 */
		static void VecHdiv(TPZFMatrix<REAL> & coord,TPZFMatrix<REAL> &NormalVec,TPZVec<int> & VectorSide);
		/** @brief Computes the vecorial product of the two vectors*/ 
		static void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);
		/** @brief Computes normal vector to plane determinated by three points */
		static void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);
        
		/* brief compute the coordinate of a point given in parameter space */
        void VecHdiv(const TPZGeoEl &gel,TPZFMatrix<REAL> &NormalVec,TPZVec<int> & VectorSide) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            VecHdiv(coord,NormalVec,VectorSide);
        }

        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
           @param matid material id of the element
           @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
           @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int64_t>& nodeindexes,
										  int matid,
										  int64_t& index);
                virtual int ClassId() const;
                void Read(TPZStream& buf, void* context);
                void Write(TPZStream& buf, int withclassid) const;
	};
    
    template<class T>
    inline void TPZGeoQuad::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1];
        
        phi(0,0) = 0.25*(1.-qsi)*(1.-eta);
        phi(1,0) = 0.25*(1.+qsi)*(1.-eta);
        phi(2,0) = 0.25*(1.+qsi)*(1.+eta);
        phi(3,0) = 0.25*(1.-qsi)*(1.+eta);
        
        dphi(0,0) = 0.25*(eta-1.);
        dphi(1,0) = 0.25*(qsi-1.);
        
        dphi(0,1) = 0.25*(1.-eta);
        dphi(1,1) =-0.25*(1.+qsi);
        
        dphi(0,2) = 0.25*(1.+eta);
        dphi(1,2) = 0.25*(1.+qsi);
        
        dphi(0,3) =-0.25*(1.+eta);
        dphi(1,3) = 0.25*(1.-qsi);
        
        
    }
    
    template<class T>
    inline void TPZGeoQuad::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        TPZFNMatrix<4,T> phi(NNodes,1);
        TPZFNMatrix<8,T> dphi(2,NNodes);
        TShape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 4; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
        
    }
    
    template<class T>
    inline void TPZGeoQuad::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        int space = nodes.Rows();
        int ncol = nodes.Cols();
        
        gradx.Resize(space,2);
        gradx.Zero();
        
#ifdef PZDEBUG
        if(/* nrow != 3 || */ ncol  != 4){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be spacex4." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<3,T> phi(4,1);
        TPZFNMatrix<6,T> dphi(2,4);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < space; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
                gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
            }
        }
        
    }
    
    
};

#endif 
