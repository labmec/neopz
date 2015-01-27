//---------------------------------------------------------------------------

#ifndef TSWXGridDataToPZH
#define TSWXGridDataToPZH
#include "map"
#include "pzvec.h"
class TPZGeoMesh;
class TSWXGridData;

namespace swx{

  TPZGeoMesh * ConvertGridDataToPZ( TSWXGridData & grid, std::map<int,std::string> &matId2GroupName,
                                    TPZVec<int> &PZGeoElIndex2ItascaID   );

  void ConvertPZToGridData(TPZGeoMesh * gmesh, TSWXGridData & grid);

};
//---------------------------------------------------------------------------
#endif
