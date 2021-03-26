#include "TPZParallelUtils.h"

template int pzutils::AtomicAdd<int>(int&, int);
template int64_t pzutils::AtomicAdd<int64_t>(int64_t&, int64_t);
template float pzutils::AtomicAdd<float>(float&, float);
template double pzutils::AtomicAdd<double>(double&, double);
template long double pzutils::AtomicAdd<long double>(long double&, long double);
