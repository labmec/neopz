#include "TPZTopologyUtils.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzpyramid.h"
#include "tpzprism.h"

namespace pztopology{
    REAL GetTolerance(){return gTolerance;}

    void SetTolerance(const REAL &tol){
        if(tol > 0) gTolerance = tol;
        else {
            typedef std::numeric_limits< REAL > dbl;

            std::cout.precision(dbl::max_digits10);
            std::cout<<"Invalid tolerance parameter for topologies. Trying to set: "<<tol<<std::endl;
            std::cout<<"This value will be ignored. Tolerance is set at: "<<gTolerance<<std::endl;
        }
    }

    template<class Topology>
    void GetPermutation(const int permutationIndex, TPZVec<int> &permutation){
        #ifdef PZDEBUG
        if(permutationIndex < 0 || permutationIndex >= Topology::NPermutations){
            PZError<<"GetPermutation: invalid parameter: permutationIndex = "<<permutationIndex<<std::endl;
            PZError<<"Aborting..."<<std::endl;
            DebugStop();
        }
        #endif
        if(permutation.size() != Topology::NSides)  permutation.Resize(Topology::NSides,-1);
        for(int i = 0; i < Topology::NSides; i++) permutation[i] = Topology::fPermutations[permutationIndex][i];
    }

}

template void pztopology::GetPermutation<pztopology::TPZPoint>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZLine>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZTriangle>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZQuadrilateral>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZTetrahedron>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZCube>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZPrism>(const int permutationIndex, TPZVec<int> &permutation);
template void pztopology::GetPermutation<pztopology::TPZPyramid>(const int permutationIndex, TPZVec<int> &permutation);



