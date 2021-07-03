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

class TPZGeoElSide;
class TPZGeoEl;

#include "pzlog.h"

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

        bool ResetBlendConnectivity(const int64_t &side, const int64_t &index){
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
        
        void Read(TPZStream &buf, void *context) override{
            Topology::Read(buf, context);
            buf.Read(fNodeIndexes, NNodes);
        }
        
        void Write(TPZStream &buf, int withclassid) const override{
            Topology::Write(buf, withclassid);
            buf.Write(fNodeIndexes, NNodes);
        }
        
        void Initialize(TPZVec<int64_t> &nodeindexes)
        {
            int64_t nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
#ifndef PZNODEBUG
            if(nodeindexes.NElements() != N)
            {
                std::stringstream sout;
                sout << __PRETTY_FUNCTION__ << " Nodeindexes have wrong size " << nodeindexes.NElements() << " but should be " << N;
                std::cout << sout.str().c_str() << std::endl;
                DebugStop();
            }
#endif
            memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int64_t));
            int64_t i;
            for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
            
        }
        
        void Initialize(TPZGeoEl *)
        {
            
        }
        
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
        int ClassId() const override;
        
    protected:
    };
    
    template<int N, class Topology>
    int TPZNodeRep<N, Topology>::ClassId() const{
        return Hash("TPZNodeRep") ^ Topology::ClassId() << 1 ^ (N << 2);
    }
};

#include "pznoderep.h.h"

#endif
