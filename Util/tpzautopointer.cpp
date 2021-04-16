/** 
 * @file 
 * @brief Initialize the std::recursive_mutex g_ap_mut and a few template instantiations. 
 */

#include "tpzautopointer.h"
//namespace pzinternal{
//    std::recursive_mutex g_ap_mut;
//    std::mutex g_diag_mut;
//}


#include "pzmatrix.h"

template class TPZAutoPointer<TPZMatrix<float>>;
template class TPZAutoPointer<TPZMatrix<double>>;
template class TPZAutoPointer<TPZMatrix<long double>>;
template class TPZAutoPointer<TPZMatrix<std::complex<float> >>;
template class TPZAutoPointer<TPZMatrix<std::complex<double> >>;
template class TPZAutoPointer<TPZMatrix<std::complex<long double> >>;

#include "pzfunction.h"
template class TPZAutoPointer<TPZFunction<float>>;
template class TPZAutoPointer<TPZFunction<double>>;
template class TPZAutoPointer<TPZFunction<long double>>;
template class TPZAutoPointer<TPZFunction<std::complex<float> >>;
template class TPZAutoPointer<TPZFunction<std::complex<double> >>;
template class TPZAutoPointer<TPZFunction<std::complex<long double> >>;

#include "TPZRefPattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
class TPZGeoMesh;
template class TPZAutoPointer<TPZGeoMesh>;
class TPZCompMesh;
template class TPZAutoPointer<TPZCompMesh>;
class TPZRefPattern;
template class TPZAutoPointer<TPZRefPattern>;
