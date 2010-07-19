//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{Edimar Cesar Rylo} <%{cesar@labmec.fec.unicamp.br}>, (C) %{2005}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <fstream>
#include <iostream>

#include <pzgeoel.h>
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <TPZRefPattern.h>


void QuadTwoAdjacentRibRefine(TPZGeoMesh *gmesh)
{
  //malha 3 triangulos
  const int nelem = 5;
  //nmero de n�
  const int ntotal_coord = 6;
  TPZVec<REAL> new_node_coord(ntotal_coord,0.);
  REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} ,
                                              { 1.,0.,0.} ,
                                              { 1.,1.,0.} ,
                                              { 0.,1.,0.} ,
                                              {0.5,0.,0.} ,
                                              {1.,0.5,0.}
                                             };

  int Connect[nelem][4] = {
                            {0,1,2,3},
                            {0,4,3,-1},
                            {4,1,5,-1},
                            {5,2,3,-1},
                            {3,4,5,-1}
                          };
  int nConnect[nelem] = {4,3,3,3,3};

  // criar um objeto tipo malha geometrica
  TPZGeoMesh *Wmesh = new TPZGeoMesh();
  // criar nos
  int i,j;
  for(i=0; i<(ntotal_coord); i++)
  {
    int nodind = Wmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(3);
    for (j=0; j<3; j++) {
      new_node_coord[j] = Total_Nodes_Coord[i][j];
    }
    Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
  }

  //int index=0;
  TPZGeoEl *gel[nelem];
  for(i=0;i<nelem;i++)
  {
    TPZVec<int> indices(nConnect[i]);
    for(j=0;j<nConnect[i];j++)
    {
      indices[j] = Connect[i][j];
    }
    int index;
    if(nConnect[i] == 3)
    {
      //std::cout << "Creating tria with cornernodes = " << indices << std::endl;
      gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
    }
    else
    {
      //cout << "Creating tria with cornernodes = " << indices << endl;
      gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
    }
  }

  //Wmesh->Print(cout);
  TPZRefPattern  *patt = new TPZRefPattern(*Wmesh) ;
  patt->SetId(0);
  std::ofstream teste("qua2adjribref.txt");
//  patt->CreateFile(teste);
//  patt->WritePattern(teste);
  delete Wmesh;
  Wmesh = 0;
  patt->InsertPermuted(/**gmesh*/);  
}
