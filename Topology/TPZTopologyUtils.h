#ifndef TOPOLOGYUTILS_H
#define TOPOLOGYUTILS_H
#include "pzreal.h"
#include <limits>
#include <pzerror.h>
#include "pzvec.h"

template<class T>
class TPZVec;

namespace pztopology{
//    class Settings{
//    private:
//        Settings() = default;
//        REAL gTolerance;
//    public:
//        static Settings &GetSettings(){static Settings fOnlyInstance;return fOnlyInstance;}
//        Settings(Settings const&)   = delete;
//        void operator=(Settings const&) = delete;
//    };// if, in the future, there are more topology settings to be adjusted, this model of singleton can be used.

    typedef std::numeric_limits< REAL > dbl;
    static REAL gTolerance = pow(10,(-1 * (dbl::max_digits10- 5)));

    REAL GetTolerance();

    void SetTolerance(const REAL &tol);

    template<class Topology>
    extern void GetPermutation(const int permutationIndex, TPZVec<int> &permutation);
}

#endif