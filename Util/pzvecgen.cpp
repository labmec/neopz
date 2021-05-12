/** 
 * @file 
 * @brief Creates a vectors to integer, REAL, char and TPZString. 
 */

#include "pzvec.h"

template class TPZVec<float>;
template class TPZVec<float * >;
template class TPZVec<double>;
template class TPZVec<double * >;
template class TPZVec<long double>;
template class TPZVec<long double * >;
template class TPZVec<int>;
template class TPZVec<int64_t>;
template class TPZVec<int64_t *>;
template class TPZVec<int *>;
template class TPZVec<char *>;
template class TPZVec<void *>;
template class TPZVec<char>;


class TPZGeoEl;
class TPZGeoNode;
struct TPZGeoNodeBC;
struct TPZGeoElBC;
template class TPZVec<TPZGeoEl *>;
template class TPZVec<TPZGeoNode *>;
template class TPZVec<TPZGeoNodeBC *>;
template class TPZVec<TPZGeoElBC *>;
#include "pzgeoelside.h"
template class TPZVec<TPZGeoElSide>;


#include "pzcompel.h"
class TPZConnect;
struct TPZConnectBC;
template class TPZVec<TPZCompEl *>;
template class TPZVec<TPZConnect *>;
template class TPZVec<TPZConnectBC *>;
template class TPZVec<TPZCompElSide>;

class TPZMaterial;
class TPZBndCond;
template class TPZVec<TPZMaterial *>;
template class TPZVec<TPZBndCond *>;
#include "TPZMaterialDataT.h"
template class TPZVec<TPZMaterialDataT<STATE>>;
template class TPZVec<TPZMaterialDataT<CSTATE>>;

struct TPZElementMatrix;
template class TPZVec<TPZElementMatrix *>;

class TPZGraphEl;
class TPZGraphNode;
class TPZCompMesh;
template class TPZVec<TPZGraphEl *>;
template class TPZVec<TPZGraphNode *>;
template class TPZVec<TPZCompMesh *>;