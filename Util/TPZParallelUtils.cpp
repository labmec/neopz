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

template void pzutils::AtomicAdd<int>(int&, int);
template void pzutils::AtomicAdd<int64_t>(int64_t&, int64_t);
template void pzutils::AtomicAdd<float>(float&, float);
template void pzutils::AtomicAdd<double>(double&, double);
template void pzutils::AtomicAdd<long double>(long double&, long double);
template void pzutils::AtomicAdd<std::complex<float>>(std::complex<float>&, std::complex<float>);
template void pzutils::AtomicAdd<std::complex<double>>(std::complex<double>&, std::complex<double>);
template void pzutils::AtomicAdd<std::complex<long double>>(std::complex<long double>&, std::complex<long double>);
