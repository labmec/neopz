#include "TPZParallelUtils.h"

template int AtomicAdd<int>(int&, int);
template int64_t AtomicAdd<int64_t>(int64_t&, int64_t);
template float AtomicAdd<float>(float&, float);
template double AtomicAdd<double>(double&, double);
template long double AtomicAdd<long double>(long double&, long double);
