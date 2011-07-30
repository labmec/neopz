//
// C++ Interface: pznoderep
//
// Description:
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef PZNODEREPH
#define PZNODEREPH

#include "pzvec.h"

#include <map>
#include <memory.h>

#include "pzmanvector.h"
#include "pztrnsform.h"

class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr lognoderep(Logger::getLogger("pz.geom.tpznoderep"));
#endif

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @ingroup topology
	 */
	const double pzgeom_TPZNodeRep_tol = 1.E-6;
	
	template<int N, class Topology>
	class TPZNodeRep : public Topology
	{
		
	private:
		
		/** @brief Verifies if pt (in parametric domain of the side) is within boundaries
		 */
		bool IsInSideParametricDomain(int side, TPZVec<REAL> &pt, REAL tol);
		
	public:
		
		virtual void SetNeighbourInfo(int side, TPZGeoElSide &neigh, TPZTransform &trans) {
			std::cout << "Element that is NOT TPZGeoBlend trying to Set Neighbour Information on Geometric Mesh!\n";
			std::cout << "See TPZGeoElRefLess::SetNeighbourInfo() Method!\n";
			DebugStop();
		}
		
		bool IsLinearMapping() const { return true; }
		bool IsGeoBlendEl() const 
		{ 
			return false; 
		}
		
		static const int NNodes=N;
		/** @brief Index of the nodes of the element */
		int fNodeIndexes[N];
		/**
		 * @brief Constructor with list of nodes
		 */
		TPZNodeRep(TPZVec<int> &nodeindexes)
		{
			int nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
			memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int));
			int i;
			for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
		}
		
		/** @brief Empty constructor */
		TPZNodeRep()
		{
			int i;
			for(i=0; i<N; i++) fNodeIndexes[i]=-1;
		}
		
		/** @brief Constructor with node map */
		TPZNodeRep(const TPZNodeRep &cp,
				   std::map<int,int> & gl2lcNdMap);
		
		/** @brief Copy constructor */
		TPZNodeRep(const TPZNodeRep<N,Topology> &cp)
		{
			memcpy(fNodeIndexes,cp.fNodeIndexes,N*sizeof(int));
		}
		
		void Initialize(TPZVec<int> &nodeindexes)
		{
			int nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
#ifndef NDEBUG
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
			memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int));
			int i;
			for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
			
		}
		void Initialize(TPZGeoEl *refel)
		{
		}
		/**
		 * Create an element along the side with material id bc
		 */
		/*
		 TPZGeoEl *CreateBCGeoEl(TPZGeoEl *orig, int side,int bc)
		 {
		 int ns = orig->NSideNodes(side);
		 TPZManVector<int> nodeindices(ns);
		 int in;
		 for(in=0; in<ns; in++)
		 {
		 nodeindices[in] = orig->SideNodeIndex(side,in);
		 }
		 int index;
		 TPZGeoEl *newel = CreateGeoElementPattern(*(orig->Mesh()),Topology::Type(side),nodeindices,bc,index);
		 TPZGeoElSide me(orig,side);
		 TPZGeoElSide newelside(newel,newel->NSides()-1);
		 newelside.InsertConnectivity(me);
		 return newel;
		 }
		 */
		/** Verifies if the parametric point pt is in the element parametric domain
		 */
		/*
		 bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6){
		 const bool result = Topology::IsInParametricDomain(pt, tol);
		 return result;
		 }//method
		 */
		
		/** @brief Projects point pt (in parametric coordinate system) in the element parametric domain.
		 * @return Returns the side where the point was projected.
		 * @note Observe that if the point is already in the parametric domain, the method will return
		 * NSides() - 1
		 */
		int ProjectInParametricDomain(TPZVec<REAL> &pt, TPZVec<REAL> &ptInDomain){
			const int nsides = Topology::NSides;
			if(this->IsInParametricDomain(pt,0.)){///it is already in the domain
				ptInDomain = pt;
				return nsides-1;
			}//if
			
			int winnerSide = -1;
			REAL winnerDistance = 1e12;
			TPZManVector<REAL,3> pt1(pt.NElements()), pt2(pt.NElements());
			for(int is = 0; is < nsides-1; is++){
				
				///Go from NSides-1 to side is
				TPZTransform T1 = Topology::SideToSideTransform(nsides-1, is);
				T1.Apply(pt,pt1);
				
				///Check if the point is within side boundaries
				bool IsInSideDomain = this->IsInSideParametricDomain(is,pt1,0.);
				if(!IsInSideDomain) continue;
				
				///Come back from side is to NSides-1
				TPZTransform T2 = Topology::SideToSideTransform(is,nsides-1);
				T2.Apply(pt1,pt2);
				
				///Compare original to mapped point
				REAL distance = 0.;
				for(int i = 0; i < pt.NElements(); i++){
					REAL val = pt[i]-pt2[i];
					distance += val*val;
				}//i
				distance = sqrt(distance);
				
				///The closest side point to the original is the projected point
				if(distance < winnerDistance){
					winnerDistance = distance;
					winnerSide = is;
					ptInDomain = pt2;
				}
			}//for is
			
			return winnerSide;
			
		}//method
		
		void Print(std::ostream &out) const
		{
			int nn;
			out << "Nodeindices :";
			for(nn=0; nn<N; nn++)
			{
				out << fNodeIndexes[nn] << ' ';
			}
			out << std::endl;
		}
		
    protected:
		/**
		 * @brief This method is redefined in TPZGeoTriangle, TPZGeoPrism, TPZGeoTetrahedra, TPZGeoPyramid
		 * to fix singularity problems when using MapToSide() method!
		 */
		static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint)
		{
			ChangedPoint.Resize(OriginalPoint.NElements(),0.);
			ChangedPoint = OriginalPoint;
		}
		
	};
	
};

#include "pznoderep.h.h"

#endif
