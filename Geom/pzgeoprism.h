/**
 * @file
 * @brief Contains the TPZGeoPrism class which implements the geometry of a prism element.
 */

#ifndef TPZGEOPRISMH
#define TPZGEOPRISMH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzprism.h"
#include "pzfmatrix.h"

class TPZGeoEl;
class TPZGeoMesh;

#include <string>


namespace pzgeom {
    
    /**
     * @ingroup geometry
     * @brief Implements the geometry of a prism element. \ref geometry "Geometry"
     */
    class TPZGeoPrism : public TPZNodeRep<6, pztopology::TPZPrism>
    {
    public:
        /** @brief Number of corner nodes */
        enum {NNodes = 6};
        /** @brief Constructor with list of nodes */
        TPZGeoPrism(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoPrism::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPrism>(nodeindexes)
        {
        }
        
        /** @brief Empty constructor */
        TPZGeoPrism() : TPZRegisterClassId(&TPZGeoPrism::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPrism>()
        {
        }
        
        /** @brief Constructor with node map */
        TPZGeoPrism(const TPZGeoPrism &cp,
                    std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoPrism::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPrism>(cp,gl2lcNdMap)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoPrism(const TPZGeoPrism &cp) : TPZRegisterClassId(&TPZGeoPrism::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoPrism(const TPZGeoPrism &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoPrism::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
        {
        }
        
        static bool IsLinearMapping(int side)
        {
            return true;
        }
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "Prism";}
        
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
         * element based on the current geometric element,
         * a side and a boundary condition number
         */
        static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
        
        virtual int ClassId() const;
        
        void Read(TPZStream& buf, void* context);
        
        void Write(TPZStream& buf, int withclassid) const;
        
    protected:
        /**
         * @brief This method apply an infinitesimal displacement in some points
         * to fix singularity problems when using MapToSide() method!
         */
        /**
         * This points are CornerNodes, when projected in the opposing side
         */
        static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);
        
        
    public:
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        
        /** @brief Creates a geometric element according to the type of the father element */
        static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
                                          TPZVec<int64_t>& nodeindexes,
                                          int matid, int64_t& index);
    };
    
    template<class T>
    inline void TPZGeoPrism::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1] , zeta  = loc[2];
        
        phi(0,0)  = .5*(1.-qsi-eta)*(1.-zeta);
        phi(1,0)  = .5*qsi*(1.-zeta);
        phi(2,0)  = .5*eta*(1.-zeta);
        phi(3,0)  = .5*(1.-qsi-eta)*(1.+zeta);
        phi(4,0)  = .5*qsi*(1.+zeta);
        phi(5,0)  = .5*eta*(1.+zeta);
        
        dphi(0,0) = -.5*(1.-zeta);
        dphi(1,0) = -.5*(1.-zeta);
        dphi(2,0) = -.5*(1.-qsi-eta);
        dphi(0,1) =  .5*(1.-zeta);
        dphi(1,1) =  .0;
        dphi(2,1) = -.5*qsi;
        dphi(0,2) =  .0;
        dphi(1,2) =  .5*(1.-zeta);
        dphi(2,2) = -.5*eta;
        dphi(0,3) = -.5*(1.+zeta);
        dphi(1,3) = -.5*(1.+zeta);
        dphi(2,3) =  .5*(1.-qsi-eta);
        dphi(0,4) =  .5*(1.+zeta);
        dphi(1,4) =  .0;
        dphi(2,4) =  .5*qsi;
        dphi(0,5) =  .0;
        dphi(1,5) =  .5*(1.+zeta);
        dphi(2,5) =  .5*eta;
        
    }
    
    template<class T>
    inline void TPZGeoPrism::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        TPZFNMatrix<6,T> phi(6,1);
        TPZFNMatrix<18,T> dphi(3,6);
        TShape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 6; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
    }
    
    
    template<class T>
    inline void TPZGeoPrism::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 || ncol  != 6){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x6." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<6,T> phi(6,1);
        TPZFNMatrix<18,T> dphi(3,6);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 6; i++)
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
