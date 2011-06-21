#include "tpzchangeel.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>

#include "tpzmathtools.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZGeoElement.h.h"
#include "pzgeoelside.h"
#include "pzstack.h"
#include "tpzgeoelrefpattern.h"

#include <sstream>
using namespace std;
using namespace pzgeom;

TPZChangeEl::TPZChangeEl()
{
}

TPZChangeEl::~TPZChangeEl()
{
}

TPZGeoEl * TPZChangeEl::ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex)
{
     TPZGeoEl *OldElem = NULL;
     TPZGeoEl *NewElem = NULL;
     OldElem = Mesh->ElementVec()[ElemIndex];
     if(OldElem)
     {
          if(OldElem->TypeName() == "Triangle" && OldElem->NNodes() == 3)
          {
               /** Creating Midnodes */
               TPZGeoNode Node3, Node4, Node5;
               TPZVec<REAL> Coord3(3), Coord4(3), Coord5(3);

               /** Setting Midnodes Coordinates */
               TPZGeoNode Node0 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(0)];
               TPZGeoNode Node1 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(1)];
               TPZGeoNode Node2 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(2)];
               for(int i = 0; i < 3; i++)
               {
                    Coord3[i] = (Node0.Coord(i) + Node1.Coord(i)) / 2.;
                    Coord4[i] = (Node1.Coord(i) + Node2.Coord(i)) / 2.;
                    Coord5[i] = (Node2.Coord(i) + Node0.Coord(i)) / 2.;
               }
               Node3.SetCoord(Coord3); Node4.SetCoord(Coord4); Node5.SetCoord(Coord5);

               /** Setting Midnodes Id's */
               int NewNodeId = Mesh->CreateUniqueNodeId();
               Mesh->SetNodeIdUsed(NewNodeId + 2);
               Node3.SetNodeId(NewNodeId); Node4.SetNodeId(NewNodeId+1); Node5.SetNodeId(NewNodeId+2);

               /** Allocating Memory for MidNodes and Pushing Them */
               TPZVec <int> NodesSequence(6) ;
               NodesSequence[0]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(0);
               NodesSequence[1]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(1);
               NodesSequence[2]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(2);
               NodesSequence[3]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[3]] = Node3;
               NodesSequence[4]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[4]] = Node4;
               NodesSequence[5]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[5]] = Node5;

               /** Inserting New Element in Mesh and Deleting Old Element */
               
               NewElem = new TPZGeoElRefPattern<TPZQuadraticTrig> (OldElem->Id(),NodesSequence,OldElem->MaterialId(),*Mesh);

               TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);

               delete OldElem;
               OldElem = NULL;
          }

          else if(OldElem->TypeName() == "Quad" && OldElem->NNodes() == 4)
          {
               /** Creating Midnodes */
               TPZGeoNode Node4, Node5, Node6, Node7;
               TPZVec<REAL> Coord4(3), Coord5(3), Coord6(3), Coord7(3);

               /** Setting Midnodes Coordinates */
               TPZGeoNode Node0 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(0)];
               TPZGeoNode Node1 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(1)];
               TPZGeoNode Node2 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(2)];
               TPZGeoNode Node3 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(3)];
               for(int i = 0; i < 3; i++)
               {
                    Coord4[i] = (Node0.Coord(i) + Node1.Coord(i)) / 2.;
                    Coord5[i] = (Node1.Coord(i) + Node2.Coord(i)) / 2.;
                    Coord6[i] = (Node2.Coord(i) + Node3.Coord(i)) / 2.;
                    Coord7[i] = (Node3.Coord(i) + Node0.Coord(i)) / 2.;
               }
               Node4.SetCoord(Coord4); Node5.SetCoord(Coord5); Node6.SetCoord(Coord6); Node7.SetCoord(Coord7);

               /** Setting Midnodes Id's */
               int NewNodeId = Mesh->CreateUniqueNodeId();
               Mesh->SetNodeIdUsed(NewNodeId + 3);
               Node4.SetNodeId(NewNodeId); Node5.SetNodeId(NewNodeId+1); Node6.SetNodeId(NewNodeId+2); Node7.SetNodeId(NewNodeId+3);

               /** Allocating Memory for MidNodes and Pushing Them */
               TPZVec <int> NodesSequence(8) ;
               NodesSequence[0]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(0);
               NodesSequence[1]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(1);
               NodesSequence[2]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(2);
               NodesSequence[3]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(3);
               NodesSequence[4]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[4]] = Node4;
               NodesSequence[5]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[5]] = Node5;
               NodesSequence[6]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[6]] = Node6;
               NodesSequence[7]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[7]] = Node7;

               /** Inserting New Element in Mesh and Deleting Old Element */
               NewElem = new TPZGeoElRefPattern<TPZQuadraticQuad> (OldElem->Id(),NodesSequence,OldElem->MaterialId(),*Mesh);
               
               TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);              

               delete OldElem;
               OldElem = NULL;
          }
          else if(OldElem->NNodes() == 2)
          {
               /** Creating Midnodes */
               TPZGeoNode Node2;
               TPZVec<REAL> Coord2(3);

               /** Setting Midnodes Coordinates */
               TPZGeoNode Node0 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(0)];
               TPZGeoNode Node1 = Mesh->NodeVec()[Mesh->ElementVec()[ElemIndex]->NodeIndex(1)];
               for(int i = 0; i < 3; i++)
               {
                    Coord2[i] = (Node0.Coord(i) + Node1.Coord(i)) / 2.;
               }
               Node2.SetCoord(Coord2);

               /** Setting Midnodes Id's */
               int NewNodeId = Mesh->CreateUniqueNodeId();
               Mesh->SetNodeIdUsed(NewNodeId);
               Node2.SetNodeId(NewNodeId);

               /** Allocating Memory for MidNodes and Pushing Them */
               TPZVec <int> NodesSequence(3) ;
               NodesSequence[0]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(0);
               NodesSequence[1]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(1);
               NodesSequence[2]  = Mesh->NodeVec().AllocateNewElement(); Mesh->NodeVec()[NodesSequence[2]] = Node2;               

               /** Inserting New Element in Mesh and Deleting Old Element */
               NewElem = new TPZGeoElRefPattern<TPZQuadraticLine> (OldElem->Id(),NodesSequence,OldElem->MaterialId(),*Mesh);
               
               TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);              

               delete OldElem;
               OldElem = NULL;
          }
          else { cout << "Element type don't recognized!\nSee ChangeToQuadratic Method!\n"; exit(-1);}
     }
     OldElem = NULL;
     return NewElem;
}

TPZGeoEl * TPZChangeEl::ChangeToLinear(TPZGeoMesh *Mesh, int ElemIndex)
{
     TPZGeoEl *OldElem = NULL;
     TPZGeoEl *NewElem = NULL;
     OldElem = Mesh->ElementVec()[ElemIndex];
     if(OldElem)
     {
          if(OldElem->TypeName() == "Triangle" && OldElem->NNodes() == 6)
          {
               TPZVec <int> NodesSequence(3) ;
               NodesSequence[0]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(0);
               NodesSequence[1]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(1);
               NodesSequence[2]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(2);

               /** Inserting New Element in Mesh and Deleting Old Element */
               
               NewElem = new TPZGeoElRefPattern<TPZGeoTriangle> (OldElem->Id(),NodesSequence,OldElem->MaterialId(),*Mesh);
               
               TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);

               delete OldElem;
               OldElem = NULL;
          }
          else if(OldElem->TypeName() == "Quad" && OldElem->NNodes() == 8)
          {
               TPZVec <int> NodesSequence(4) ;
               NodesSequence[0]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(0);
               NodesSequence[1]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(1);
               NodesSequence[2]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(2);
               NodesSequence[3]  = Mesh->ElementVec()[ElemIndex]->NodeIndex(3);

               /** Inserting New Element in Mesh and Deleting Old Element */
               
               NewElem = new TPZGeoElRefPattern<TPZGeoTriangle> (OldElem->Id(),NodesSequence,OldElem->MaterialId(),*Mesh);

               TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);

               delete OldElem;
               OldElem = NULL;
          }
          else { cout << "Element type don't recognized!\nSee ChangeToLinear Method!\n"; exit(-1);}
     }
     OldElem = NULL;
     return NewElem;
}

TPZGeoEl* TPZChangeEl::QuarterPoints(TPZGeoMesh *Mesh, int ElemIndex, int side)
{
     TPZGeoEl *OldElem = Mesh->ElementVec()[ElemIndex];
    if(!OldElem)
    {
        DebugStop();
    }
     if(OldElem->TypeName() == "Triangle" && OldElem->NNodes() == 3)
     {
          if(side < 0 || side > 6) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
          Mesh->SetNodeIdUsed(2);
          TPZGeoEl * newEl = ChangeToQuadratic(Mesh,ElemIndex);
          switch(side)
          {
               case 0:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                      SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                      SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 1:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                      SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 2:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 3:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 4:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 5:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(3)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               case 6:
               // Do nothing!
               break;
          }
          return newEl;
     }
     else if(OldElem->TypeName() == "Quad" && OldElem->NNodes() == 4)
     {
          if(side < 0 || side > 8) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
          Mesh->SetNodeIdUsed(3);
          TPZGeoEl* newEl = ChangeToQuadratic(Mesh,ElemIndex);
          switch(side)
          {
               case 0:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 1:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));
               }
               break;
               case 2:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 3:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 4:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 5:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 6:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(5)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(7)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               break;
               case 7:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(4)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));

                    Mesh->NodeVec()[newEl->NodeIndex(6)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(2)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(3)].Coord(i));
               }
               case 8:
               // Do nothing!
               break;
          }
        return newEl;
     }
     else if(OldElem->NNodes() == 2)
     {
          if(side < 0 || side > 2) { cout << "Invalid Side to Compute Quarter Points!\nSee QuarterPoints Method!\n"; exit(-1);}
          TPZGeoEl* newEl = ChangeToQuadratic(Mesh,ElemIndex);
          switch(side)
          {
               case 0:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(2)].
                    SetCoord(i, 0.75*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.25*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
               }
               break;
               case 1:
               for(int i = 0; i < 3; i++)
               {
                    Mesh->NodeVec()[newEl->NodeIndex(2)].
                    SetCoord(i, 0.25*Mesh->NodeVec()[newEl->NodeIndex(0)].Coord(i) + 0.75*Mesh->NodeVec()[newEl->NodeIndex(1)].Coord(i));
               }
               break;
          }
        return newEl;
     }
     else { cout << "Element type don't recognized!\nSee QuarterPoints Method!\n"; exit(-1); }
}

TPZGeoEl * TPZChangeEl::ChangeToGeoBlend(TPZGeoMesh *Mesh, int ElemIndex)
{
  TPZGeoEl * OldElem = Mesh->ElementVec()[ElemIndex];
  if(!OldElem){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " - NULL geometric element.\n";
    return NULL;
  }
  
  const int nnodes = OldElem->NNodes();
  TPZManVector<int> nodeindexes(nnodes);
  for(int i = 0; i < nnodes; i++) nodeindexes[i] = OldElem->NodeIndex(i);
  int newindex;
  TPZGeoEl * NewElem = Mesh->CreateGeoBlendElement(OldElem->Type(),nodeindexes,OldElem->MaterialId(),newindex);
  TPZChangeEl::AdjustNeighbourhood(OldElem,NewElem);
  NewElem->BuildBlendConnectivity();

  delete OldElem;

  return NewElem;
}///method

void TPZChangeEl::AdjustNeighbourhood(TPZGeoEl* OldElem, TPZGeoEl*NewElem){
  for(int j = 0; j < OldElem->NSides(); j++)
  { 
    TPZGeoElSide currSide (OldElem,j);
    TPZGeoElSide newElSide (NewElem,j);
    TPZGeoElSide neighbour (currSide.Neighbour());
    
    if (neighbour == currSide)
    {
        newElSide.SetNeighbour(newElSide);
        continue;
    }
    newElSide.SetNeighbour(neighbour);
  
    while (neighbour != currSide)
    {
      TPZGeoElSide next = neighbour.Neighbour();
      if (next == currSide)
      {
          neighbour.SetNeighbour(newElSide);
          break;
      }
      neighbour = neighbour.Neighbour();
    }
    
  }
}///method






// #include "pznoderep.h.h"


//          template class pzgeom::TPZNodeRep<8,TPZGeoBlend<TPZGeoCube> >;
//          template class pzgeom::TPZNodeRep<6,TPZGeoBlend<TPZGeoPrism> >;

         


