/**
 * @file
 * @brief Contains the implementation of the TPZPoint methods. 
 */

#include "tpzpoint.h"
#include "pzquad.h"
#include "pzeltype.h"

#ifdef _AUTODIFF
#include "fad.h"
#endif

#include "pzcreateapproxspace.h"

namespace pztopology {
	
	TPZIntPoints *TPZPoint::CreateSideIntegrationRule(int side, int order)
	{
		return new IntruleType(order);
	}
	
	
	MElementType TPZPoint::Type()
	{
		return EPoint;
	}
	
	MElementType TPZPoint::Type(int side)
	{
		switch(side) {
			case 0:
				return EPoint;
			default:
				return ENoType;
		}
	}
    
    int TPZPoint::NBilinearSides()
    {return 0;}
    
    template<class T>
    bool TPZPoint::MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide) {
		SidePar.Resize(0); JacToSide.Resize(0,0);
		return true;
	}
    
    void TPZPoint::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {

        if(node > NCornerNodes)
        {
            DebugStop();
        }
        nodeCoord.Resize(Dimension, 0.);
        switch (node) {
            case (0):
            {
                return;
                break;
            }
            default:
            {
                DebugStop();
                break;
            }
        }
    }
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZPoint::GetTransformId(TPZVec<int64_t> &id)
	{
		return 0;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZPoint::GetTransformId(int side, TPZVec<int64_t> &id)
	{
		return 0;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
//	void TPZPoint::GetSideHDivPermutation(int side, TPZVec<int64_t> &id, TPZVec<int> &permgather)
//	{
//		permgather[0] = 0;
//		return;
//	}
    void TPZPoint::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        directions.Redim(3, 0);
        sidevectors.Resize(0);
    }
	
    int TPZPoint::ClassId() const{
        return Hash("TPZPoint");
    }

    void TPZPoint::Read(TPZStream& buf, void* context) {

    }

    void TPZPoint::Write(TPZStream& buf, int withclassid) const {

    }
    
}

template
bool pztopology::TPZPoint::MapToSide<REAL>(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);

#ifdef _AUTODIFF
template
bool pztopology::TPZPoint::MapToSide<Fad<REAL> >(int side, TPZVec<Fad<REAL> > &InternalPar, TPZVec<Fad<REAL> > &SidePar, TPZFMatrix<Fad<REAL> > &JacToSide);
#endif
