//---------------------------------------------------------------------------

#pragma hdrstop

#include "SmallerGMeshGenerator.h"

#include "pzgmesh.h"
#include "pzreal.h"
#include "pzmanvector.h"
#include "pzgnode.h"

#include "pzerror.h"

#include <iostream>
#include <algorithm>

//---------------------------------------------------------------------------
#pragma package(smart_init)

SmallerGMeshGenerator::SmallerGMeshGenerator()
{
}

SmallerGMeshGenerator::~SmallerGMeshGenerator()
{
}

TPZGeoMesh * SmallerGMeshGenerator::GetGMesh(TPZVec<REAL> & XYZanchor,
                                             double Lx, double Ly, double Lz, double Lcharact)
{
  TPZGeoMesh * gmesh = new TPZGeoMesh;

  int nX = std::max(1.,Lx/Lcharact+0.5);
  int nY = std::max(1.,Ly/Lcharact+0.5);
  int nZ = std::max(1.,Lz/Lcharact+0.5);

  double deltaX = Lx / nX;
  double deltaY = Ly / nY;
  double deltaZ = Lz / nZ;

  int nnodesX = nX + 1;
  int nnodesY = nY + 1;
  int nnodesZ = nZ + 1;

  int nnodes = nnodesX * nnodesY * nnodesZ;

  gmesh->NodeVec().Resize(nnodes);

  int n = 0;
  for(int yy = 0; yy < nnodesY; yy++)
  {
    for(int zz = 0; zz < nnodesZ; zz++)
    {
      for(int xx = 0; xx < nnodesX; xx++)
      {
        double coordX = XYZanchor[0] + xx*deltaX;
        double coordY = XYZanchor[1] + yy*deltaY;
        double coordZ = XYZanchor[2] + zz*deltaZ;

        TPZManVector<REAL> nodeCoord(3);
        nodeCoord[0] = coordX;
        nodeCoord[1] = coordY;
        nodeCoord[2] = coordZ;

        TPZGeoNode node;
        node.Initialize(nodeCoord, *gmesh);

        if(n >= gmesh->NNodes())
        {
          DebugStop();
        }
        gmesh->NodeVec()[n] = node;
        n++;
      }
    }
  }

  long el = 0;
  int nnodesPerXZplane = nnodesX*nnodesZ;
  TPZManVector<long> topol(8);
  for(int yy = 0; yy < nY; yy++)
  {
    for(int zz = 0; zz < nZ; zz++)
    {
      for(int xx = 0; xx < nX; xx++)
      {
        int _A = yy*nnodesPerXZplane + zz*nnodesX + xx;
        int _B = _A + 1;
        int _C = _B + nnodesPerXZplane;
        int _D = _A + nnodesPerXZplane;

        int _E = _A + nnodesX;
        int _F = _B + nnodesX;
        int _G = _C + nnodesX;
        int _H = _D + nnodesX;

        topol[0] = _A;
        topol[1] = _B;
        topol[2] = _C;
        topol[3] = _D;
        topol[4] = _E;
        topol[5] = _F;
        topol[6] = _G;
        topol[7] = _H;

        int matid = 1;
        gmesh->CreateGeoElement(ECube,topol,matid,el);
        el++;
      }
    }
  }

  gmesh->BuildConnectivity();
  return gmesh;
}
