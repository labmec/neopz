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
        typedef pztopology::TPZCube Top;
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

        /** @brief Compute x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodecoordinates,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodecoordinates,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
       
        
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        
        public:
int ClassId() const override;
        void Read(TPZStream &buf, void *context) override;
        void Write(TPZStream &buf, int withclassid) const override;

    public:
        /** @brief Creates a geometric element according to the type of the father element */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t>& nodeindexes,
        //                                   int matid,
        //                                   int64_t& index);
        
    };

    
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
