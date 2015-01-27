//---------------------------------------------------------------------------

#include "TSWXGridDataToPZ.h"
#include "TSWXIMPGRIDData.h"
#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzmanvector.h"

using namespace pzgeom;
using namespace pzrefine;
using namespace pzshape;

TPZGeoMesh * swx::ConvertGridDataToPZ( TSWXGridData & grid, std::map<int,std::string> &matId2GroupName,
                                       TPZVec<int> &PZGeoElIndex2ItascaID ){
  TPZGeoMesh * gmesh = new TPZGeoMesh();

  ///nodes
  const unsigned int nnodes = grid.GetSizeGridPoint();
  gmesh->NodeVec().Resize( nnodes );
  TPZManVector<REAL,3> coord(3);
  for(unsigned int i = 0; i < nnodes; i++){
    const TSWXGridPoint & point = grid.GridPoints()[i];
    unsigned int id = point.fId;
#ifdef DEBUG
  unsigned int index = point.fIndex;
  if(index != i) DebugStop();
#endif
    for(int j = 0; j < 3; j++) coord[j] = point.fX[j];
    gmesh->NodeVec()[i].Initialize(id,coord,*gmesh);
  }///for inodes

  ///elements
  ///mapping PZindex -> id
  std::map<int,int> mappingId2PZIndex;
  PZGeoElIndex2ItascaID.Resize(grid.GetSizeB8Element() + grid.GetSizeW6Element() );
  PZGeoElIndex2ItascaID.Fill(-1);

  ///elements : B8
  {///escopo B8
    const unsigned int nB8els = grid.GetSizeB8Element();
    const int matid = 8511965;///depois vou trocar com base no zgroup
    long index;
    TPZManVector<long,8> cornerindexes(8);
    const int mapping[8] = {1,2,5,3,4,7,8,6};
    for(unsigned int i = 0; i < nB8els; i++){
      const TSWXGridElement<8> & b8el = grid.B8Elements()[i];
      for(unsigned int k = 0; k < 8; k++){
        const int nodeid = b8el.fIndices[ mapping[k]-1 ];
        const int nodeindex = grid.GridPointMapping()[ nodeid ];
        cornerindexes[k] = nodeindex;
      }

      TPZGeoElRefPattern < TPZGeoCube > *gel = new TPZGeoElRefPattern < TPZGeoCube >(cornerindexes, matid, *gmesh, index);
      gel->SetId( b8el.fId );
      mappingId2PZIndex[ b8el.fId ] = index;
      PZGeoElIndex2ItascaID[ index ] = b8el.fId;
    }///for ielem
  }///escopo B8

  ///elements : W6
  {///escopo W6
    const unsigned int nW6els = grid.GetSizeW6Element();
    const int matid = 8511965;///depois vou trocar com base no zgroup
    long index;
    TPZManVector<long,6> cornerindexes(6);
    const int mapping[6] = {1,2,4,3,5,6};
    for(unsigned int i = 0; i < nW6els; i++){
      const TSWXGridElement<6> & w6el = grid.W6Elements()[i];
      for(unsigned int k = 0; k < 6; k++){
        const int nodeid = w6el.fIndices[ mapping[k]-1 ];
        const int nodeindex = grid.GridPointMapping()[ nodeid ];
        cornerindexes[k] = nodeindex;
      }

      TPZGeoElRefPattern < TPZGeoPrism > *gel = new TPZGeoElRefPattern < TPZGeoPrism >(cornerindexes, matid, *gmesh, index);
      gel->SetId( w6el.fId );
      mappingId2PZIndex[ w6el.fId ] = index;
      PZGeoElIndex2ItascaID[ index ] = w6el.fId;
    }///for ielem
  }///escopo W6

  ///zgroups
  matId2GroupName.clear();
  const unsigned int ngroups = grid.GetSizeZoneGroup();
  for(unsigned int i = 0; i < ngroups; i++){
    const TSWXZoneGroup & group = grid.ZoneGroups()[i];
    const int matid = i+1;
    matId2GroupName[ matid ] = group.fGroupName;
    const unsigned int nels = group.fEls.size();
    for(int iel = 0; iel < nels; iel++){
      const int elId = group.fEls[iel];
      const int elIndex = mappingId2PZIndex[elId];
#ifdef DEBUG
      if(gmesh->ElementVec()[elIndex] == NULL) DebugStop();
#endif
      gmesh->ElementVec()[elIndex]->SetMaterialId(matid);
    }///iel
  }///igrupos

  gmesh->BuildConnectivity();

  return gmesh;
}///method


void swx::ConvertPZToGridData(TPZGeoMesh * gmesh, TSWXGridData & grid)
{
  int nnodes = gmesh->NNodes();

  for(unsigned int n = 0; n < nnodes; n++)
  {
    const double x = gmesh->NodeVec()[n].Coord(0);
    const double y = gmesh->NodeVec()[n].Coord(1);
    const double z = gmesh->NodeVec()[n].Coord(2);
    grid.AddGridPoint(n,x,y,z);
  }

  int nelem = gmesh->NElements();
  for(int el = 0; el < nelem; el++)
  {
    TPZGeoEl * gel = gmesh->ElementVec()[el];

    if(gel && gel->Dimension() == 3 && gel->HasSubElement() == false)
    {
      if(gel->NNodes() == 6)//Prism
      {
        const unsigned int id = gel->Id();
        std::vector<unsigned int> indices(6);
        indices[0] = gel->NodeIndex(0);
        indices[1] = gel->NodeIndex(1);
        indices[2] = gel->NodeIndex(3);
        indices[3] = gel->NodeIndex(2);
        indices[4] = gel->NodeIndex(4);
        indices[5] = gel->NodeIndex(5);

        grid.AddW6Element(id,indices);
      }
      else if(gel->NNodes() == 8)//Cube
      {
        const unsigned int id = gel->Id();
        std::vector<unsigned int> indices(8);
        indices[0] = gel->NodeIndex(0);
        indices[1] = gel->NodeIndex(1);
        indices[2] = gel->NodeIndex(3);
        indices[3] = gel->NodeIndex(4);
        indices[4] = gel->NodeIndex(2);
        indices[5] = gel->NodeIndex(7);
        indices[6] = gel->NodeIndex(5);
        indices[7] = gel->NodeIndex(6);

        grid.AddB8Element(id,indices);
      }
      else
      {
        continue;
      }
    }
  }

  //Faltam os GROUPS!!!
}

