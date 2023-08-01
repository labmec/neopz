#include "TPZParallelUtils.h"
#ifdef USING_MKL
#include <mkl.h>
#endif

namespace pzutils{
    int SetNumThreadsLocalMKL(const int nt)
    {
#ifdef USING_MKL
        return mkl_set_num_threads_local(nt);
#endif
        return 0;
    }
};

template int pzutils::AtomicAdd<int>(int&, int);
template int64_t pzutils::AtomicAdd<int64_t>(int64_t&, int64_t);
template float pzutils::AtomicAdd<float>(float&, float);
template double pzutils::AtomicAdd<double>(double&, double);
template long double pzutils::AtomicAdd<long double>(long double&, long double);
