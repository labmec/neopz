/**
 * @file
 * @brief Contains the TPZNodeRep class which implements ...
 * Clase intermediaria que guarda 
 */

#ifndef PZNODEREPH
#define PZNODEREPH

#include "pzvec.h"

#include <map>
#include <memory.h>

#include "pzmanvector.h"
#include "pztrnsform.h"

class TPZGeoMesh;
class TPZGeoElSide;
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr lognoderep(Logger::getLogger("pz.geom.tpznoderep"));
#endif

namespace pzgeom {
    
    /**
     * @ingroup topology
     * @brief Initializing tolerance to TPZNodeRep.
     */
    const double pzgeom_TPZNodeRep_tol = 1.E-12;
    
    /**
     * @ingroup geometry
     * @ingroup topology
     * @brief Implements ... \ref geometry "Geometry" \ref topology "Topology"
     */
    template<int N, class Topology>
    class TPZNodeRep : public Topology
    {
        
    private:
        
        /** @brief Verifies if pt (in parametric domain of the side) is within boundaries */
        bool IsInSideParametricDomain(int side, TPZVec<REAL> &pt, REAL tol);
        
    public:
        
        virtual void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform<> &trans) {
            std::cout << "Element that is NOT TPZGeoBlend trying to Set Neighbour Information on Geometric Mesh!\n";
            std::cout << "See TPZGeoElRefLess::SetNeighbourInfo() Method!\n";
            DebugStop();
        }
        
        bool IsLinearMapping() const;
        
        bool IsGeoBlendEl() const
        {
            return false;
        }
        
        static const int NNodes=N;
        /** @brief Node indexes of the element */
        int64_t fNodeIndexes[N];
        /** @brief Constructor with list of nodes */
        TPZNodeRep(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZNodeRep::ClassId)
        {
            int nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
#ifdef PZDEBUG
            if(nn<N)
            {
                DebugStop();
            }
#endif
            memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int64_t));
            int i;
            for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
        }
        
        /** @brief Empty constructor */
        TPZNodeRep() : TPZRegisterClassId(&TPZNodeRep::ClassId)
        {
            int64_t i;
            for(i=0; i<N; i++) fNodeIndexes[i]=-1;
        }
        
        /** @brief Constructor with node map */
        TPZNodeRep(const TPZNodeRep &cp,
                   std::map<int64_t,int64_t> & gl2lcNdMap);
        
        /** @brief Copy constructor */
        TPZNodeRep(const TPZNodeRep<N,Topology> &cp) : TPZRegisterClassId(&TPZNodeRep::ClassId)
        {
            memcpy(fNodeIndexes,cp.fNodeIndexes,N*sizeof(int64_t));
        }
        
        void Read(TPZStream &buf, void *context) {
            Topology::Read(buf, context);
            buf.Read(fNodeIndexes, NNodes);
        }
        
        virtual void Write(TPZStream &buf, int withclassid) const { 
            Topology::Write(buf, withclassid);
            buf.Write(fNodeIndexes, NNodes);
        }
        
        void Initialize(TPZVec<int64_t> &nodeindexes)
        {
            int64_t nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
#ifndef NODEBUG
            if(nodeindexes.NElements() != N)
            {
                std::stringstream sout;
                sout << __PRETTY_FUNCTION__ << " Nodeindexes have wrong size " << nodeindexes.NElements() << " but should be " << N;
#ifdef LOG4CXX
                LOGPZ_ERROR(lognoderep,sout.str().c_str());
#else
                std::cout << sout.str().c_str() << std::endl;
#endif
            }
#endif
            memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int64_t));
            int64_t i;
            for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
            
        }
        void Initialize(TPZGeoEl *refel)
        {
        }
        
        /** @brief Gets the corner node coordinates in coord */
        void CornerCoordinates(const TPZGeoEl &gel, TPZFMatrix<REAL> &coord) const
        {
            TPZGeoNode *np;
            int i,j;
            for(i=0;i<NNodes;i++) {
                np = gel.NodePtr(i);
                for(j=0;j<3;j++) {
                    coord(j,i) = np->Coord(j);
                }
            }
        }
        
        /**
         * @brief Projects point pt (in parametric coordinate system) in the element parametric domain.
         * @return Returns the side where the point was projected.
         * @note Observe that if the point is already in the parametric domain, the method will return
         * \f$ NSides() - 1 \f$
         */
        int ProjectInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain){
            const int nsides = Topology::NSides;
            if(this->IsInParametricDomain(qsi,0.)){///it is already in the domain
                qsiInDomain = qsi;
                return nsides-1;
            }//if
            
            int winnerSide = -1;
            REAL winnerDistance = 1e12;
            TPZManVector<REAL,3> pt1(qsi.NElements()), pt2(qsi.NElements());
            for(int is = 0; is < nsides-1; is++){
                
                ///Go from NSides-1 to side is
                TPZTransform<> T1 = Topology::SideToSideTransform(nsides-1, is);
                T1.Apply(qsi,pt1);
                
                ///Check if the point is within side boundaries
                bool IsInSideDomain = this->IsInSideParametricDomain(is,pt1,0.);
                if(!IsInSideDomain) continue;
                
                ///Come back from side is to \f$ NSides-1 \f$
                TPZTransform<> T2 = Topology::SideToSideTransform(is,nsides-1);
                T2.Apply(pt1,pt2);
                
                ///Compare original to mapped point
                REAL distance = 0.;
                for(int i = 0; i < qsi.NElements(); i++){
                    REAL val = qsi[i]-pt2[i];
                    distance += val*val;
                }//i
                distance = sqrt(distance);
                
                ///The closest side point to the original is the projected point
                if(distance < winnerDistance){
                    winnerDistance = distance;
                    winnerSide = is;
                    qsiInDomain = pt2;
                }
            }//for is
            
            return winnerSide;
            
        }//method
        
        // project the point towards the center of the element
        // find the intersecting side
        int ProjectBissectionInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain){
            const int nsides = Topology::NSides;
            const int dim = Topology::Dimension;
            
            qsiInDomain.Resize(dim);
            
            if(this->IsInParametricDomain(qsi,0.))///it is already in the domain
            {
                qsiInDomain = qsi;
                return nsides-1;
            }
            
            REAL tol = 1.e-10;
            
            ///first, will be made a project to center direction
            TPZVec<REAL> OutPt(qsi), InnPt(dim);
            Topology::CenterPoint(nsides-1,InnPt);
            
            REAL dist = 0.;
            for(int c = 0; c < dim; c++)
            {
                dist += (InnPt[c] - OutPt[c])*(InnPt[c] - OutPt[c]);
            }
            dist = sqrt(dist);
            
            while(dist > tol)
            {
                for(int c = 0; c < dim; c++)
                {
                    qsiInDomain[c] = (InnPt[c] + OutPt[c])/2.;
                }
                if(this->IsInParametricDomain(qsiInDomain,0.))
                {
                    InnPt = qsiInDomain;
                }
                else
                {
                    OutPt = qsiInDomain;
                }
                dist = 0.;
                for(int c = 0; c < dim; c++)
                {
                    dist += (InnPt[c] - OutPt[c])*(InnPt[c] - OutPt[c]);
                }
                dist = sqrt(dist);
            }
            
            ///found in witch side the projection belongs
            int winnerSide = -1;
            TPZManVector<REAL,3> pt1(dim), pt2(dim);
            for(int is = 0; is < nsides-1; is++)
            {
                ///Go orthogonally from \f$ NSides-1 \f$ to side is
                TPZTransform<> T1 = Topology::SideToSideTransform(nsides-1, is);
                T1.Apply(qsiInDomain,pt1);
                
                ///Come back from side is to \f$ NSides-1 \f$
                TPZTransform<> T2 = Topology::SideToSideTransform(is,nsides-1);
                T2.Apply(pt1,pt2);
                
                ///Compare ptInDomain to transformed point
                dist = 0.;
                for(int c = 0; c < dim; c++)
                {
                    dist += (qsiInDomain[c]-pt2[c]) * (qsiInDomain[c]-pt2[c]);
                }//i
                dist = sqrt(dist);
                
                ///Closest side
                if(dist < tol)
                {
                    winnerSide = is;
                    qsiInDomain = pt2;
                    break;
                }
            }//for is
            
            return winnerSide;
            
        }//method
        
        void Print(std::ostream &out) const
        {
            int nn;
            out << "Nodeindices: ";
            for(nn=0; nn<N; nn++)
            {
                out << fNodeIndexes[nn] << ' ';
            }
            out << std::endl;
        }
        
    public:
        virtual int ClassId() const;
        
    protected:
        /**
         * @brief This method is redefined in TPZGeoTriangle, TPZGeoPrism, TPZGeoTetrahedra, TPZGeoPyramid \n
         * to fix singularity problems when using MapToSide() method!
         */
        static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
        {
            ChangedPoint.Resize(OriginalPoint.NElements(),0.);
            ChangedPoint = OriginalPoint;
        }
        template<class T>
        static void TFixSingularity(int side, TPZVec<T>& OriginalPoint, TPZVec<T>& ChangedPoint)
        {
            ChangedPoint.Resize(OriginalPoint.NElements(),0.);
            ChangedPoint = OriginalPoint;
        }
    };
    
    template<int N, class Topology>
    int TPZNodeRep<N, Topology>::ClassId() const{
        return Hash("TPZNodeRep") ^ Topology::ClassId() << 1 ^ (N << 2);
    }
};

#include "pznoderep.h.h"

#endif
