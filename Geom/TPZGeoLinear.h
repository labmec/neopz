
#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH

/**
 * @file TPZGeoLinear
 * @brief Implements the geometry of a one dimensional linear element.
 */

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"
#include "tpzline.h"

#include <string>

class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom
{
    
    
    /**
     * @ingroup geometry
     * @brief Implements the geometry of a one dimensional linear element. \ref geometry "Geometry"
     */
    class TPZGeoLinear : public TPZNodeRep<2, pztopology::TPZLine> {
        
    public:
        
        /** @brief Number of corner nodes */
        enum {NNodes = 2};
        
        /** @brief Constructor with list of nodes */
        TPZGeoLinear(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoLinear::ClassId), TPZNodeRep<NNodes, pztopology::TPZLine>(nodeindexes)
        {
        }
        
        /** @brief Empty constructor */
        TPZGeoLinear() : TPZRegisterClassId(&TPZGeoLinear::ClassId),TPZNodeRep<NNodes, pztopology::TPZLine>()
        {
        }
        
        /** @brief Constructor with node map */
        TPZGeoLinear(const TPZGeoLinear &cp, std::map<int64_t,int64_t> & gl2lcNdMap) : 
        TPZRegisterClassId(&TPZGeoLinear::ClassId), TPZNodeRep<NNodes, pztopology::TPZLine>(cp,gl2lcNdMap)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoLinear(const TPZGeoLinear &cp) : TPZRegisterClassId(&TPZGeoLinear::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
        {
        }
        
        /** @brief Copy constructor */
        TPZGeoLinear(const TPZGeoLinear &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoLinear::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
        {
        }
        
        /** @brief answer if the element side is a linear map */
        static bool IsLinearMapping(int side)
        {
            return true;
        }
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "Linear";}
        
        /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        
        /* @brief Compute X mapping from local parametric coordinates */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &x) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,x);
        }
        
        /** @brief Compute gradient of X mapping from local parametric coordinates */
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
        
        /** @brief Compute X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        /** @brief Compute the jacoabina associated to the x mapping from local parametric coordinates  */
        static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
                             TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
        
        /**
         * @brief Method which creates a geometric boundary condition
         * element based on the current geometric element, \n
         * a side and a boundary condition number
         */
        static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
        
        public:
        virtual int ClassId() const;
        void Read(TPZStream& buf, void* context);
        void Write(TPZStream& buf, int withclassid) const;

        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        
        /** @brief Creates a geometric element according to the type of the father element */
        static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh,
                                          MElementType type,
                                          TPZVec<int64_t>& nodeindexes,
                                          int matid,
                                          int64_t& index);
        
    };
    
    template<class T>
    inline void TPZGeoLinear::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T x = loc[0];
        phi(0,0) = (1.0-x)/2.;
        phi(1,0) = (1.0+x)/2.;
        dphi(0,0) = -0.5;
        dphi(0,1) = 0.5;
    }
    
    template<class T>
    inline void TPZGeoLinear::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        T xi = loc[0];
        int nrow = nodes.Rows();
        for(int i = 0; i < nrow; i++)
        {
            x[i] = nodes.GetVal(i,0)*(1.-xi)*0.5+nodes.GetVal(i,1)*(1.+xi)*0.5;
        }
        
    }
    
    
    template<class T>
    inline void TPZGeoLinear::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
        
        gradx.Resize(nrow,1);
        gradx.Zero();
        
        TPZFNMatrix<3,T> phi(2,1);
        TPZFNMatrix<6,T> dphi(2,2);
        TShape(loc,phi,dphi);
        for(int i = 0; i < ncol; i++)
        {
            for(int j = 0; j < nrow; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);

            }
        }
        
    }
    
};

#endif
