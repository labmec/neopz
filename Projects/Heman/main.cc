#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "pzelcq2d.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include <TPZRefPattern.h>

#include "pzmaterial.h"
#include "pzelasmat.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
//#include "c0-simplequad.cpp" 
#include "c0-simplequad.cpp"

/*#include <dl_entities.h>
#include <dl_dxf.h>
#include <dl_attributes.h>
*/
#include <vector>
 

using namespace std ;

TPZGeoMesh *ReadGeoMesh(ifstream &arq);


void neighbourcheck (int ribsideindex , TPZGeoEl *elemento , TPZVec <int> &neighbourindex ){


  int ribsidenode1 = elemento->SideNodeIndex(ribsideindex,0);
  int ribsidenode2 = elemento->SideNodeIndex(ribsideindex,1);

  int nsides = elemento->NSides() ;
  int dimensao , node1index , node2index , k=0 ;
 
  for (int control=0;control<nsides;control++){

    dimensao = elemento->SideDimension(control);
       if ((dimensao==1)&&(control!=ribsideindex)){
	 node1index = elemento->SideNodeIndex(control,0);
         node2index = elemento->SideNodeIndex(control,1);
	     if ((node1index==ribsidenode1)||(node1index==ribsidenode2)||(node2index==ribsidenode1)||(node2index==ribsidenode2)){

	       neighbourindex[k]=control;
               k++;

             }

       }

  }

  neighbourindex.Resize(k);

}






void siderecog (TPZGeoMesh *Gmesh, int indice_elemento){

  // Para reconhecer o lado dividido ao meio, basta trabalharmos com o primeiro elemento da malha (pai).
  TPZGeoEl *elemento = Gmesh->ElementVec()[0];
  // Acima é gerado o GelEl do elemento pai .
  
  int nnodes = elemento->NNodes();
  int nsides = elemento->NSides();
  int lado = 0;
  int dimensao , no1index , no2index , ladodividido ;
  double detcheck ; 
  double smallnumber  = 1.e-10;
  int ribsideindex , ribsidenode1 , ribsidenode2 ; 
  TPZVec <double> no1coord(3);
  TPZVec <double> no2coord(3);
  TPZVec <double> nomediocoord(3);
  TPZManVector <int> neighbourindex(nsides);

  int nomedioindex = nnodes ;
  TPZGeoNode *nomedio = &(Gmesh->NodeVec()[nomedioindex]);

  // No loop abaixo testamos todos os sides, e naqueles que são veificados como sendo arestas, fazemos o teste de alinhamento dos pontos .

  for (int controle=0;controle<nsides;controle++){
    dimensao = elemento->SideDimension(controle);
    if (dimensao==1){
      lado++;
      no1index = elemento->SideNodeIndex(controle,0);
      no2index = elemento->SideNodeIndex(controle,1);
      TPZGeoNode *no1 = &(elemento->Mesh()->NodeVec()[no1index]);//elemento->NodePtr(no1index);
      TPZGeoNode *no2 = &(elemento->Mesh()->NodeVec()[no2index]);//elemento->NodePtr(no2index);
      for (int i=0;i<3;i++){
        no1coord[i] = no1->Coord(i);
        no2coord[i] = no2->Coord(i);
        nomediocoord[i] = nomedio->Coord(i);
      }
      // Calcularemos o determinante da matriz contendo as coordenadas de 
      // três pontos para checar alinhamento (os dois extremos do lado e o nó médio)
      detcheck = ((nomediocoord[1]*no1coord[2])+(nomediocoord[2]*no2coord[1])+(no1coord[1]*no2coord[2])-(no1coord[2]*no2coord[1])-(nomediocoord[2]*no1coord[1])-(nomediocoord[1]*no2coord[2]))*((nomediocoord[1]*no1coord[2])+(nomediocoord[2]*no2coord[1])+(no1coord[1]*no2coord[2])-(no1coord[2]*no2coord[1])-(nomediocoord[2]*no1coord[1])-(nomediocoord[1]*no2coord[2]))+((nomediocoord[2]*no1coord[0])+(nomediocoord[0]*no2coord[2])+(no1coord[2]*no2coord[0])-(no1coord[0]*no2coord[2])-(nomediocoord[0]*no1coord[2])-(nomediocoord[2]*no2coord[0]))*((nomediocoord[2]*no1coord[0])+(nomediocoord[0]*no2coord[2])+(no1coord[2]*no2coord[0])-(no1coord[0]*no2coord[2])-(nomediocoord[0]*no1coord[2])-(nomediocoord[2]*no2coord[0]))+((nomediocoord[0]*no1coord[1])+(nomediocoord[1]*no2coord[0])+(no1coord[0]*no2coord[1])-(no1coord[1]*no2coord[0])-(nomediocoord[1]*no1coord[0])-(nomediocoord[0]*no2coord[1]))*((nomediocoord[0]*no1coord[1])+(nomediocoord[1]*no2coord[0])+(no1coord[0]*no2coord[1])-(no1coord[1]*no2coord[0])-(nomediocoord[1]*no1coord[0])-(nomediocoord[0]*no2coord[1]));
  
      // Caso o alinhamento seja verificado, guardamos algumas referências abaixo.
      if (fabs(detcheck)<smallnumber){
        ladodividido = lado;
        ribsideindex = controle;
        ribsidenode1 = no1index;
        ribsidenode2 = no2index;
        neighbourcheck (ribsideindex,elemento, neighbourindex );
      }
    }
  }

  int considesnumber = neighbourindex.NElements();

  cout << "Index geral do lado dividido entre os sides :" << ribsideindex << endl;
  cout << "Numero do lado, contando a partir de 1 :" << ladodividido << endl;
  cout << "Indice de um dos nos do lado :"  << ribsidenode1 << endl ;
  cout << "Indice de outro no do lado :" << ribsidenode2 << endl ;
  cout << endl;
  cout << endl;
  cout << "Existem " << considesnumber << " lados conectados ao RibSide." << endl;
  cout << "Estes lados são os sides de índice :" << endl;
  for (int p=0;p<considesnumber;p++){
    cout << neighbourindex[p] << endl;
  }
}


void linemarker (TPZGeoMesh *Gmesh, map < set<int> , TPZRefPattern* > &MyMap, TPZGeoEl *elemento, TPZStack<int> &refinesides){
  set<int> lados;
  int lados_nelements=0;
  
  int nsides = elemento->NSides();
  int ncornernodes = elemento->NCornerNodes();
  
  
  TPZVec<int> cornerindic(elemento->NCornerNodes(),0);
  
  int no;
  for(no=0; no<elemento->NCornerNodes(); no++) 
  {
    TPZGeoElSide gels(elemento,no);
    TPZGeoElSide neigh = gels.Neighbour();
    while(neigh != gels)
    {
      if(neigh.Element()->MaterialId() < 0) 
      {
        cornerindic[no] = 1;
        break;
      }
      neigh = neigh.Neighbour();
    }
  }
  int is;
  for(is= elemento->NCornerNodes(); is<elemento->NSides(); is++)
  {
    if(elemento->SideDimension(is) != 1) continue;
    TPZStack<int> smallsides;
    elemento->LowerDimensionSides(is,smallsides);
    int nsmall = smallsides.NElements();
    int iss;
    for(iss=0; iss<nsmall; iss++)
    {
      if(cornerindic[smallsides[iss]] != 0) break;
    }
    if(iss < nsmall)
    {
      TPZGeoElSide gels(elemento,is);
      TPZGeoElSide neigh = gels.Neighbour();
      while(neigh != gels)
      {
        if(neigh.Element()->MaterialId() < 0) 
        {
          break;
        }
        neigh = neigh.Neighbour();
      }
      if(neigh == gels) 
      {
        refinesides.Push(is);    
        lados.insert(is);
        lados_nelements++;
        cout << "Sides to insert\n";
        cout << is << endl;
        cout << endl;
      } 
    }
  }   
 
  /* Observação: No STACK são retornados os lados a serem refinados mas com a numeração local de cada lado em relação ao elemento (não a numeração global) */
  // Abaixo serão criados os RefPatterns e associados ao set no mapa.
 
  if (ncornernodes==3)
  {
    cout << "TRIANGLE SELECTED TO CREATE A NEW REFINEMENT PATTERN" << endl;
    if(lados_nelements==1)
    {
      int permut = refinesides[refinesides.NElements()-1]-3;
      
      //malha 2 triangulos
      const int nelem = 2;
      //número de nós 
      const int ntotal_coord = 4;
      TPZVec<REAL> new_node_coord(ntotal_coord,0.);
      
      REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} , 
                                                  { 1.,0.,0.} ,
                                                  { 0.,1.,0.} , 
                                                  { 0.,0.,0.} };
      
      switch (refinesides[refinesides.NElements()-1])
      {
        case(3):
          Total_Nodes_Coord[3][0] = 0.5 ;
          Total_Nodes_Coord[3][1] = 0. ;
          Total_Nodes_Coord[3][2] = 0. ;
          break;
        
        case(4):
          Total_Nodes_Coord[3][0] = 0.5 ;
          Total_Nodes_Coord[3][1] = 0.5 ;
          Total_Nodes_Coord[3][2] = 0. ;
          break;
          
        case(5):
          Total_Nodes_Coord[3][0] = 0. ;
          Total_Nodes_Coord[3][1] = 0.5 ;
          Total_Nodes_Coord[3][2] = 0. ;
          break;
        
        default:
          cout << "Erro " << endl;
          cout << endl ;
      }

      int Connect[nelem][4] = { {(permut+1)%3,(permut+2)%3,3,-1},
                                {(permut)%3,(permut+2)%3,3,-1} };
      int nConnect[nelem] = {3,3};
  
      // criar um objeto tipo malha geometrica
      TPZGeoMesh *Wmesh = new TPZGeoMesh();
  
      // criar nos
      int i,j;
      for(i=0; i<(ntotal_coord); i++) {
        int nodind = Wmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        for (j=0; j<3; j++) {
          new_node_coord[j] = Total_Nodes_Coord[i][j];
        }
        Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
      }
      
      int index=0;
              
      TPZVec <int> indicesfather (3);
      indicesfather [0]= 0;
      indicesfather [1]= 1;
      indicesfather [2]= 2;
      TPZGeoEl *father = Wmesh->CreateGeoElement(ETriangle,indicesfather,1,index);
      // criação dos elementos
      TPZGeoEl *gel[nelem];
    
      for(i=0;i<nelem;i++) {  
        TPZVec<int> indices(nConnect[i]);
        for(j=0;j<nConnect[i];j++) {
          indices[j] = Connect[i][j];
        }
        int index;
        switch (nConnect[i]){
        case (4): 
          cout << "Creating quad with cornernodes = " << indices << endl;
          gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        case(3):
          cout << "Creating tria with cornernodes = " << indices << endl;
          gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        default:
          cout << "Erro : elemento não implementado" << endl;
        }
      }
      
      Wmesh->Print(cout);
      TPZRefPattern  *patt = new TPZRefPattern(Wmesh) ;
/*      cout << "Refinement pattern data:\n";
      cout << "NNodes = " << patt->NNodes() << endl;
      cout << "NSubel = " << patt->NSubElements() << endl;
      cout << endl;*/
      Gmesh->InsertRefPattern(patt);
      //MyMap[lados] = patt ;
      elemento->SetRefPattern(patt);
      TPZVec <TPZGeoEl*> ElVec ;
      elemento->Divide(ElVec);
    }

    if(lados_nelements==2)
    {
      cout << "Number of sides selected to divide = 2" << endl;
      int maior, menor;
       
      if ((refinesides[refinesides.NElements()-1])>(refinesides[refinesides.NElements()-2]))
      {
        maior = refinesides[refinesides.NElements()-1];
        menor = refinesides[refinesides.NElements()-2];
      }
      else
      {
        maior = refinesides[refinesides.NElements()-2];
        menor = refinesides[refinesides.NElements()-1];
      }
      int snlocid_comum;
      int snlocid_dif_menor, snlocid_dif_maior;
      
      int k,v ;
      for (k=0; k<2; k++)
      {
        for (v=0; v<2; v++)
	{
          if (elemento->SideNodeLocIndex(menor,k)==elemento->SideNodeLocIndex(maior,v))
          {
            snlocid_comum = elemento->SideNodeLocIndex(menor,k);
          }
        }
      }
      switch (snlocid_comum)
      {
        case(0):
	  snlocid_dif_menor=1;
	  snlocid_dif_maior=2;
        break;

        case(1):
          snlocid_dif_menor=0;
          snlocid_dif_maior=2;
        break;
  
        case(2):
          snlocid_dif_menor=1;
          snlocid_dif_maior=0;
        break;
        
        default:
          cout << "Error..." << endl;
      }

      //malha 2 elementos
      const int nelem = 2;
      //número de nós 
      const int ntotal_coord = 5;
      TPZVec<REAL> new_node_coord(ntotal_coord,0.);
     
      REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} , 
                                                  { 1.,0.,0.} ,
	  	                                  { 0.,1.,0.} , 
                                                  { 0.,0.,0.} , 
                                                  { 0.,0.,0.} };
      switch (menor)
      {
        case(3):
          Total_Nodes_Coord[3][0] = 0.5 ;
          Total_Nodes_Coord[3][1] = 0. ;
          Total_Nodes_Coord[3][2] = 0. ;
        break;
       
        case(4):
          Total_Nodes_Coord[3][0] = 0.5 ;
          Total_Nodes_Coord[3][1] = 0.5 ;
          Total_Nodes_Coord[3][2] = 0. ;
        break;
         
        case(5):
          Total_Nodes_Coord[3][0] = 0. ;
          Total_Nodes_Coord[3][1] = 0.5 ;
          Total_Nodes_Coord[3][2] = 0. ;
        break;
       
        default:
            cout << "Erro " << endl;
            cout << endl ;
      }
     
      switch (maior)
      {
        case(3):
          Total_Nodes_Coord[4][0] = 0.5 ;
          Total_Nodes_Coord[4][1] = 0. ;
          Total_Nodes_Coord[4][2] = 0. ;
        break;
        
        case(4):
          Total_Nodes_Coord[4][0] = 0.5 ;
          Total_Nodes_Coord[4][1] = 0.5 ;
          Total_Nodes_Coord[4][2] = 0. ;
        break;
          
        case(5):
          Total_Nodes_Coord[4][0] = 0. ;
          Total_Nodes_Coord[4][1] = 0.5 ;
          Total_Nodes_Coord[4][2] = 0. ;
        break;
        
        default:
          cout << "Erro " << endl;
          cout << endl ;
      }
     
      int Connect[nelem][4] = { {3,4,snlocid_comum, -1},
	                        {3,4,snlocid_dif_menor,snlocid_dif_maior} };
      int nConnect[nelem] = {3,4};
  
      // criar um objeto tipo malha geometrica
      TPZGeoMesh *Wmesh = new TPZGeoMesh();
  
      // criar nos
      int i,j;
      for(i=0; i<(ntotal_coord); i++) {
        int nodind = Wmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        for (j=0; j<3; j++) {
          new_node_coord[j] = Total_Nodes_Coord[i][j];
        }
        Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
      }
      int index=0;
      TPZVec <int> indicesfather (3);
      indicesfather [0]= 0;
      indicesfather [1]= 1;
      indicesfather [2]= 2;
      TPZGeoEl *father = Wmesh->CreateGeoElement(ETriangle,indicesfather,1,index);
      // criação dos elementos
      TPZGeoEl *gel[nelem];

      for(i=0;i<nelem;i++) {  
        TPZVec<int> indices(nConnect[i]);
        for(j=0;j<nConnect[i];j++) {
          indices[j] = Connect[i][j];
        }
        int index;
        switch (nConnect[i]){
        case (4): 
          gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        case(3):
          gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        default:
          cout << "Erro : elemento não implementado" << endl;
        }
      }
 
      cout << "Triangle to debug..." << endl;
      Wmesh->Print(cout);
       
     
      TPZRefPattern  *patt = new TPZRefPattern(Wmesh) ;
  /*     cout << patt->NNodes() << endl;
      cout << patt->NSubElements() << endl;
      cout << endl;*/
      // TPZGeoMesh, string Nome, int código
      //patt->Print1(Wmesh, cout);
      Gmesh->InsertRefPattern(patt);
      //MyMap[lados] = patt ;             
      elemento->SetRefPattern(patt);
      TPZVec <TPZGeoEl*> ElVec ;
      elemento->Divide(ElVec);
    }
  }
  //Neste else considera-se o caso de ser um quadrado
  else
  {
    if (lados_nelements==1)
    {
      int permut = refinesides[refinesides.NElements()-1]-3;
      
      //malha 3 triangulos
      const int nelem = 3;
      //número de nós 
      const int ntotal_coord = 5;
      TPZVec<REAL> new_node_coord(ntotal_coord,0.);
      
      REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} , 
                                                  { 1.,0.,0.} , 
                                                  { 1.,1.,0.} ,
                                                  { 0.,1.,0.} , 
                                                  { 0.,0.,0.} };
      
      switch (refinesides[refinesides.NElements()-1])
      {
        case(4):
          Total_Nodes_Coord[4][0] = 0.5 ;
          Total_Nodes_Coord[4][1] = 0. ;
          Total_Nodes_Coord[4][2] = 0. ;
          break;
        
        case(5):
          Total_Nodes_Coord[4][0] = 1. ;
          Total_Nodes_Coord[4][1] = 0.5 ;
          Total_Nodes_Coord[4][2] = 0. ;
          break;
          
        case(6):
          Total_Nodes_Coord[4][0] = 0.5 ;
          Total_Nodes_Coord[4][1] = 1. ;
          Total_Nodes_Coord[4][2] = 0. ;
          break;
          
        case(7):
          Total_Nodes_Coord[4][0] = 0. ;
          Total_Nodes_Coord[4][1] = 0.5 ;
          Total_Nodes_Coord[4][2] = 0. ;
          break;
        
        default:
            cout << "Erro " << endl;
            cout << endl ;
      
      }
      int Connect[nelem][4] = { {(permut+1)%4,(permut+2)%4,4,-1},
                                {(permut+2)%4,(permut+3)%4,4,-1},
                                {(permut+4)%4,(permut+3)%4,4,-1} };
      int nConnect[nelem] = {3,3,3};
    
      // criar um objeto tipo malha geometrica
      TPZGeoMesh *Wmesh = new TPZGeoMesh();
    
      // criar nos
      int i,j;
      for(i=0; i<(ntotal_coord); i++) {
        int nodind = Wmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        for (j=0; j<3; j++) {
          new_node_coord[j] = Total_Nodes_Coord[i][j];
        }
        Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
      }
  
      int index=0;
              
      TPZVec <int> indicesfather (4);
      indicesfather [0]= 0;
      indicesfather [1]= 1;
      indicesfather [2]= 2;
      indicesfather [3]= 3;
      TPZGeoEl *father = Wmesh->CreateGeoElement(EQuadrilateral,indicesfather,1,index);
      
      // criação dos elementos
      TPZGeoEl *gel[nelem];
  
      for(i=0;i<nelem;i++) {  
        TPZVec<int> indices(nConnect[i]);
        for(j=0;j<nConnect[i];j++) {
          indices[j] = Connect[i][j];
        }
        int index;
        switch (nConnect[i]){
        case (4): 
          gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        case(3):
          gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
          father->SetSubElement( i , gel[i]);      
          gel[i]->SetFather(0);
          gel[i]->SetFather(father);
          break;
        default:
          cout << "Erro : elemento não implementado" << endl;
        }
      }
                    
      Wmesh->Print(cout);
        
      TPZRefPattern  *patt = new TPZRefPattern(Wmesh) ;
/*      cout << patt->NNodes() << endl;
      cout << patt->NSubElements() << endl;
      cout << endl;*/
      // TPZGeoMesh, string Nome, int código
      //patt->Print1(Wmesh, cout);
      Gmesh->InsertRefPattern(patt);
      //MyMap[lados] = patt ;
      elemento->SetRefPattern(patt);
      TPZVec <TPZGeoEl*> ElVec ;
      elemento->Divide(ElVec);
    }
  
    if (lados_nelements==2)
    {
      if((refinesides[refinesides.NElements()-1])-(refinesides[refinesides.NElements()-2]))
      {
        //Neste caso temos que os lados marcados são opostos.
        int permut ;
        //As condições abaixo garantem que o permut será calculado com o lado marcado de maior índice.
        if (refinesides[refinesides.NElements()-1]>refinesides[refinesides.NElements()-2])
        {
          permut = (refinesides[refinesides.NElements()-1]) - 4 ;
        }
        else 
        {
          permut = (refinesides[refinesides.NElements()-2]) - 4 ;
        } 
        const int nelem = 2;
        //número de nós 
        const int ntotal_coord = 6;
        TPZVec<REAL> new_node_coord(ntotal_coord,0.);
      
        REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} , 
                                                    { 1.,0.,0.} , 
                                                    { 1.,1.,0.} ,
                                                    { 0.,1.,0.} , 
                                                    { 0.,0.,0.} , 
                                                    { 0.,0.,0.} };
        
        switch (refinesides[refinesides.NElements()-1])
        {
          case(4):
          case(6):
            Total_Nodes_Coord[4][0] = 0.5 ;
            Total_Nodes_Coord[4][1] = 0. ;
            Total_Nodes_Coord[4][2] = 0. ;
          
            Total_Nodes_Coord[5][0] = 0.5 ;
            Total_Nodes_Coord[5][1] = 1. ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
        
          case(5):
          case(7):
            Total_Nodes_Coord[4][0] = 1. ;
            Total_Nodes_Coord[4][1] = 0.5 ;
            Total_Nodes_Coord[4][2] = 0. ;
          
            Total_Nodes_Coord[5][0] = 0. ;
            Total_Nodes_Coord[5][1] = 0.5 ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
          
          default:
              cout << "Erro " << endl;
              cout << endl ;
        }
        int Connect[nelem][4] = { {(permut+1)%4,(permut+2)%4,4,5 },
                                  {(permut)%4,(permut+3)%4,4,5 } };
        int nConnect[nelem] = {4,4};
    
        // criar um objeto tipo malha geometrica
        TPZGeoMesh *Wmesh = new TPZGeoMesh();
        // criar nos
        int i,j;
        for(i=0; i<(ntotal_coord); i++) {
          int nodind = Wmesh->NodeVec().AllocateNewElement();
          TPZVec<REAL> coord(3);
          for (j=0; j<3; j++) {
            new_node_coord[j] = Total_Nodes_Coord[i][j];
          }
          Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
        }
  
        int index=0;
              
        TPZVec <int> indicesfather (4);
        indicesfather [0]= 0;
        indicesfather [1]= 1;
        indicesfather [2]= 2;
        indicesfather [3]= 3;
        TPZGeoEl *father = Wmesh->CreateGeoElement(EQuadrilateral,indicesfather,1,index);
        //cout << "Created Father Element: " << endl;
        //father->Print(cout);
        //TPZRefPattern *father_patt (this);
        // criação dos elementos
        TPZGeoEl *gel[nelem];
        for(i=0;i<nelem;i++) {  
          TPZVec<int> indices(nConnect[i]);
          for(j=0;j<nConnect[i];j++) {
            indices[j] = Connect[i][j];
          }
          switch (nConnect[i]){
          case (4): 
            gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
            father->SetSubElement( i , gel[i]);      
            gel[i]->SetFather(0);
            gel[i]->SetFather(father);
            break;
          case(3):
            gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
            gel[i]->SetFather(0);
            gel[i]->SetFather(father);
            break;
          default:
            cout << "Erro : elemento não implementado" << endl;
          }
  //          cout << "Created son element: " << i << endl;
  //          gel[i]->Print(cout);
        }
        
        //father->SetSubElementConnectivities();
        //Wmesh->BuildConnectivity();
        Wmesh->Print(cout);
        TPZRefPattern  *patt = new TPZRefPattern(Wmesh) ;
/*        cout << patt->NNodes() << endl;
        cout << patt->NSubElements() << endl;
        cout << endl;*/
        // TPZGeoMesh, string Nome, int código
        //patt->Print1(Wmesh, cout);
        Gmesh->InsertRefPattern(patt);
        //MyMap[lados] = patt ;
        elemento->SetRefPattern(patt);
        TPZVec <TPZGeoEl*> ElVec ;
        elemento->Divide(ElVec);
      }
      else
      {
        //Já neste caso iremos trabalhar com lados marcados consecutivos.
        int maior, menor;
        
        if ((refinesides[refinesides.NElements()-1])>(refinesides[refinesides.NElements()-2]))
        {
          maior = refinesides[refinesides.NElements()-1];
          menor = refinesides[refinesides.NElements()-2];
        }
        else
        {
          maior = refinesides[refinesides.NElements()-2];
          menor = refinesides[refinesides.NElements()-1];
        }
        int snlocid_comum, snlocid_dif_menor, snlocid_dif_maior, snlocid_nao_usado;
        int k, v ;
        for (k=0; k<2; k++)
        {
          for (v=0; v<2; v++)
          {
            if (elemento->SideNodeLocIndex(menor,k)==elemento->SideNodeLocIndex(maior,v))
            {
              snlocid_comum = elemento->SideNodeLocIndex(menor,k);
            }
          }
        }
        switch (snlocid_comum)
        {
          case(0):
            snlocid_dif_menor=1;
            snlocid_dif_maior=3;
            break;

          case(1):
            snlocid_dif_menor=0;
            snlocid_dif_maior=2;
            break;
            
          case(2):
            snlocid_dif_menor=1;
            snlocid_dif_maior=3;
            break;
  
          case(3):
            snlocid_dif_menor=2;
            snlocid_dif_maior=0;
            break;
            
          default:
              cout << "Error..." << endl;
        }
        for (int d=0; d<4 ; d++)
        {
          if ( (d!=snlocid_comum) && (d!=snlocid_dif_menor) && (d!=snlocid_dif_maior) )
          {
            snlocid_nao_usado = d ;
          }
        }
        
        //malha 4 triângulos
        const int nelem = 4;
        //número de nós 
        const int ntotal_coord = 6;
        TPZVec<REAL> new_node_coord(ntotal_coord,0.);
      
        REAL Total_Nodes_Coord[ntotal_coord][3] = { { 0.,0.,0.} , 
                                                    { 1.,0.,0.} , 
                                                    { 1.,1.,0.} ,
                                                    { 0.,1.,0.} , 
                                                    { 0.,0.,0.} , 
                                                    { 0.,0.,0.} };
        switch (menor)
        {
          case(4):
            Total_Nodes_Coord[4][0] = 0.5 ;
            Total_Nodes_Coord[4][1] = 0. ;
            Total_Nodes_Coord[4][2] = 0. ;
            break;
          
          case(5):
            Total_Nodes_Coord[4][0] = 1. ;
            Total_Nodes_Coord[4][1] = 0.5 ;
            Total_Nodes_Coord[4][2] = 0. ;
            break;
            
          case(6):
            Total_Nodes_Coord[4][0] = 0.5 ;
            Total_Nodes_Coord[4][1] = 1. ;
            Total_Nodes_Coord[4][2] = 0. ;
            break;
                  
          case(7):
            Total_Nodes_Coord[4][0] = 0. ;
            Total_Nodes_Coord[4][1] = 0.5 ;
            Total_Nodes_Coord[4][2] = 0. ;
            break;
              
          default:
              cout << "Erro " << endl;
              cout << endl ;
              break;
        }
      
        switch (maior)
        {
        case(4):
            Total_Nodes_Coord[5][0] = 0.5 ;
            Total_Nodes_Coord[5][1] = 0. ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
          
          case(5):
            Total_Nodes_Coord[5][0] = 1. ;
            Total_Nodes_Coord[5][1] = 0.5 ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
            
          case(6):
            Total_Nodes_Coord[5][0] = 0.5 ;
            Total_Nodes_Coord[5][1] = 1. ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
                  
          case(7):
            Total_Nodes_Coord[5][0] = 0. ;
            Total_Nodes_Coord[5][1] = 0.5 ;
            Total_Nodes_Coord[5][2] = 0. ;
            break;
              
          default:
              cout << "Erro " << endl;
              cout << endl ;
              break;
        }
        int Connect[nelem][4] = { { menor,snlocid_comum,maior, -1 },
                                  { menor,maior,snlocid_nao_usado, -1 },
                                  { menor,snlocid_nao_usado,snlocid_dif_menor, -1 },
                                  { maior,snlocid_nao_usado,snlocid_dif_maior, -1 } };
        int nConnect[nelem] = {3,3,3,3};
    
        // criar um objeto tipo malha geometrica
        TPZGeoMesh *Wmesh = new TPZGeoMesh();
    
        // criar nos
        int i,j;
        for(i=0; i<(ntotal_coord); i++) {
          int nodind = Wmesh->NodeVec().AllocateNewElement();
          TPZVec<REAL> coord(3);
          for (j=0; j<3; j++) {
            new_node_coord[j] = Total_Nodes_Coord[i][j];
          }
          Wmesh->NodeVec()[nodind] = TPZGeoNode(i,new_node_coord,*Wmesh);
        }
  
        int index=0;
              
        TPZVec <int> indicesfather (4);
        indicesfather [0]= 0;
        indicesfather [1]= 1;
        indicesfather [2]= 2;
        indicesfather [3]= 3;
        TPZGeoEl *father = Wmesh->CreateGeoElement(EQuadrilateral,indicesfather,1,index,1);
        // criação dos elementos
        TPZGeoEl *gel[nelem];
  
        for(i=0;i<nelem;i++) {  
          TPZVec<int> indices(nConnect[i]);
          for(j=0;j<nConnect[i];j++) {
            indices[j] = Connect[i][j];
          }
          int index;
          switch (nConnect[i]){
          case (4): 
            gel[i] = Wmesh->CreateGeoElement(EQuadrilateral,indices,1,index,1);
            father->SetSubElement( i , gel[i]);      
            gel[i]->SetFather(0);
            gel[i]->SetFather(father);
            break;
          case(3):
            gel[i] = Wmesh->CreateGeoElement(ETriangle,indices,1,index,1);
            father->SetSubElement( i , gel[i]);      
            gel[i]->SetFather(0);
            gel[i]->SetFather(father);
            break;
          default:
            cout << "Erro : elemento não implementado" << endl;
          }
        }
        
        Wmesh->Print(cout);
        
        TPZRefPattern  *patt = new TPZRefPattern(Wmesh) ;
//         cout << patt->NNodes() << endl;
//         cout << patt->NSubElements() << endl;
//         cout << endl;
        // TPZGeoMesh, string Nome, int código
        //patt->Print1(Wmesh, cout);
        Gmesh->InsertRefPattern(patt);
        //MyMap[lados] = patt ;
        elemento->SetRefPattern(patt);
        TPZVec <TPZGeoEl*> ElVec ;
        elemento->Divide(ElVec);
      }
    }
  }
}


int main ()
{
  TPZCompMesh *cmesh = CreateSillyMesh();
  TPZGeoMesh *Gmesh = cmesh->Reference();
  
  //Gmesh->Print(cout);
  //cout << endl;

  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Hexa_Unif.rpt");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/ic/luis/Documents/ic/catalogacao/meuRefPattern/Hexaedro_5");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Piram_Rib_Side_7.rpt");
  //TPZGeoMesh *Gmesh = reftetra->Mesh();
  
  //dxfdraw(Gmesh);  
  //dxfdrawsep(Gmesh);

  //Teste da função de reconhecimento de lados
  //siderecog(Gmesh);
  //Fim do teste
 
  
  int nelements = Gmesh->NElements();
  int el;  
  map<set<int>, TPZRefPattern*> MyMap;
  TPZStack<int> refinesides;
  for (el=0; el<nelements; el++)
  {
    TPZGeoEl *elemento = Gmesh->ElementVec()[el];
    if (elemento->MaterialId()>0)
    {
      linemarker (Gmesh, MyMap, elemento, refinesides); 
    }
  }
  //int nelement_stack = refinesides.NElements();
  //int elem;
  /*for (elem=0; elem<nelement_stack;elem++)
  {
    cout << refinesides[elem] << endl ;
  
  }
  */
  //Gmesh->Print(cout);
  return 0 ;
}



