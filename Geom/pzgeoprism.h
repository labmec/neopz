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
        typedef pztopology::TPZPrism Top;
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

        
        
        /** @brief Compute x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        
        // /**
        //  * @brief Method which creates a geometric boundary condition
        //  * element based on the current geometric element,
        //  * a side and a boundary condition number
        //  */
        // static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
        
        int ClassId() const override;
        
        void Read(TPZStream &buf, void *context) override;
        
        void Write(TPZStream &buf, int withclassid) const override;
        
        
    public:
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        
        /** @brief Creates a geometric element according to the type of the father element */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t>& nodeindexes,
        //                                   int matid, int64_t& index);
    };

    
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
