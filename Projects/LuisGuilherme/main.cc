#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzintel.h"
#include "pzcompel.h"

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

/*#include <dl_entities.h>
#include <dl_dxf.h>
#include <dl_attributes.h>
*/
#include <vector>
 

using namespace std ;

TPZGeoMesh *ReadGeoMesh(ifstream &arq);



void polyfacemesh (ostream &arquivo,int meshindex){

   
    arquivo << "  0" << endl;
    arquivo << "POLYLINE" << endl;
    arquivo << "  100" << endl;
    arquivo << "AcDbEntity" << endl;
    arquivo << "  8" << endl;
    arquivo << "Mesh" << meshindex << endl;
    arquivo << "  100" << endl;
    arquivo << "AcDbPolyFaceMesh" << endl;
    arquivo << "  66" << endl;
    arquivo << "1" << endl;
    arquivo << "  10" << endl;
    arquivo << "0.0" << endl;
    arquivo << "  20" << endl;
    arquivo << "0.0" << endl;
    arquivo << "  30" << endl;
    arquivo << "0.0" << endl;
    arquivo << "  70" << endl;
    arquivo << "64" << endl;



}


void facemaker (ostream &arquivo, int faces ,TPZGeoEl *gel , int el , TPZGeoMesh *Gmesh){

  int dxfnodeindex, sidenodes , nodeindex , finalnodeindex  ;


  int nsides = gel->NSides();
  int sidedim ;
  
 
  for (int controlsides=0 ; controlsides<nsides ; controlsides++){
      sidedim = gel->SideDimension(controlsides);
      if (sidedim==2){
	   

   
     dxfnodeindex = 70 ;

     arquivo << "  0" << endl;
     arquivo << "VERTEX" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbEntity" << endl;
     arquivo << "  8" << endl;
     arquivo << "Faces"  << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbFaceRecord" << endl;
     arquivo << "  10" << endl;
     arquivo << "0.0" << endl;
     arquivo << "  20" << endl;
     arquivo << "0.0" << endl;
     arquivo << "  30" << endl;
     arquivo << "0.0" << endl;
     arquivo << "  70" << endl;
     arquivo << "128" << endl;

     //int numerodenos = gel->NNodes() ;
     sidenodes = gel->NSideNodes(controlsides) ; 
     

     for (int c = 0 ; c<sidenodes ; c++){
     //for (int c = 0 ; c<numerodenos ; c++){  
       dxfnodeindex++ ; 
       nodeindex = gel->SideNodeLocIndex(controlsides,c);
       //nodeindex = gel->NodeIndex(c);
       finalnodeindex = nodeindex + 1 ;
     
       arquivo << "  " << dxfnodeindex << endl;
       arquivo << finalnodeindex << endl ;
       

     }

   
  }

  }


}



void dxfdraw (TPZGeoMesh *Gmesh)
{
 int nnodes = Gmesh->NNodes();
 int node;
 TPZVec <double> node_coord(3);

 string nome ;

 int meshindex = 0 ;
  
  
 int el,n;
 int nelem = Gmesh->ElementVec().NElements();
 
 ofstream arquivo("prismunif1.dxf");

 arquivo << "  0" << endl ;
 arquivo << "SECTION" << endl ;
 arquivo << "  2" << endl;
 arquivo << "ENTITIES" << endl;


 for (el=0;el<nelem;el++){
    
   


   //ofstream arquivo(nome.c_str());
 
    meshindex = el ; 
  
    polyfacemesh (arquivo,meshindex);
   
   
   TPZGeoEl *gel = Gmesh->ElementVec()[el];
   if (!gel) continue;
   int ncorner = gel->NCornerNodes();

    arquivo << "  71" << endl;
    arquivo << ncorner << endl;





    int faces = 0 , nsides , controlsides = 0 , sidedim ;
 
   // Identificando o n�mero de faces do elemento.
  
    nsides = gel->NSides();
    

    for (controlsides=0 ; controlsides<nsides ; controlsides++){
      sidedim = gel->SideDimension(controlsides);
      if (sidedim==2){
	faces++;
      }
    }

   // Fim da identifica��o do n�mero de faces.


    arquivo << "  72" << endl;
    arquivo << faces << endl;
   
 





   for (n=0;n<ncorner;n++){
     arquivo << "  0" << endl;
     arquivo << "VERTEX" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbEntity" << endl;
     arquivo << "  8" << endl;
     arquivo << "Vertices" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbVertex" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbPolyFaceMeshVertex" << endl;
     
     
     TPZGeoNode *node = gel->NodePtr(n);
     
     TPZGeoNode *geonode = &(Gmesh->NodeVec()[node->Id()]);
     int i,k;

     for (int i=0; i<3 ; i++){
       node_coord[i] = geonode->Coord(i);
       k = (i*10)+10 ;       
       arquivo << k << endl;
       arquivo << node_coord[i] << endl;        
     }
     arquivo << "  70" << endl;
     arquivo << "192" << endl;



     
   }


   // Defini��o do tipo de elemento e chamada de fun��o construtora de faces

   
   facemaker (arquivo , faces , gel , el , Gmesh);


   // Fim da constru��o de faces


   arquivo << "  0" << endl;
   arquivo << "SEQEND" << endl;

 }



 // Inicio da insercao de indice dos pontos

 int f=0;

 //int numerel = Gmesh->ElementVec().NElements();

 //for (int control=0;control<numerel;control++){

 // TPZGeoEl *gel = Gmesh->ElementVec()[control];
 //if (!gel) continue;
 //int ncorner = gel->NCornerNodes();

   int ncorners = Gmesh->NNodes();

   for (int control2=0;control2<ncorners;control2++){

    arquivo << "  0" << endl;
    arquivo << "TEXT" << endl;
    arquivo << "  100"  << endl;
    arquivo << "AcDbText" << endl;
    arquivo << "  8" << endl;
    arquivo << "Indices" << endl;
    arquivo << "  62" << endl;
    arquivo << "5" << endl;
   
   
    TPZGeoNode *no = &(Gmesh->NodeVec()[control2]);

    //TPZGeoNode *node = gel->NodePtr(control2);
     
    //TPZGeoNode *geonode = &(Gmesh->NodeVec()[node->Id()]);
    int z,k;
     for (int z=0; z<3 ; z++){
       k=z;
       node_coord[z] = no->Coord(z);
       k=(z*10)+10;
       arquivo << k << endl;
       arquivo << node_coord[z] << endl;
    }

   
   arquivo << "  40" << endl;
   arquivo << "0.1" << endl;
   arquivo << "  1" << endl;
   arquivo << f << endl;
   f++ ;
  }


   //ncorners = Gmesh->NNodes();
   int numerel = Gmesh->ElementVec().NElements();
   int lay = 0 , colo = 0 ;

   for (int y=0;y<numerel;y++){

     int w=0;
     lay++ ;
     colo++ ;

     TPZGeoEl *gel = Gmesh->ElementVec()[y];
   
     int numNodes = gel->NNodes();
   
     for (int control=0;control<numNodes;control++){

     

      arquivo << "  0" << endl;
      arquivo << "TEXT" << endl;
      arquivo << "  100"  << endl;
      arquivo << "AcDbText" << endl;
      arquivo << "  8" << endl;
      arquivo << lay << endl;
      arquivo << "  62" << endl;
      arquivo << colo << endl;
      
   
   
    //TPZGeoNode *no = &(Gmesh->NodeVec()[control]);

    //TPZGeoEl *gel = Gmesh->ElementVec()[control];
    //TPZGeoNode *node = gel->NodePtr(control);
    //TPZGeoNode *geonode = &(Gmesh->NodeVec()[control]);


    TPZGeoNode *node = gel->NodePtr(control);
     
    TPZGeoNode *geonode = &(Gmesh->NodeVec()[node->Id()]);

    int z,k;
     for (int z=0; z<3 ; z++){
       node_coord[z] = geonode->Coord(z);
       k=(z*10)+10;
       arquivo << k << endl;
       arquivo << node_coord[z] << endl;
    }

   
   arquivo << "  40" << endl;
   arquivo << "0.1" << endl;
   arquivo << "  1" << endl;
   arquivo << w << endl;
   w++ ;
  }


}
  

 
   arquivo << "  0" << endl;
   arquivo << "ENDSEC" << endl;
   arquivo << "  0" << endl;
   arquivo << "EOF" << endl;

}
















/* Inicio da gera��o do arquivo DXF que gera imagens 
   das componentes filhas reduzidas */




















 void dxfdrawsep (TPZGeoMesh *Gmesh, TPZRefPattern * refPatt)
{
 int nnodes = Gmesh->NNodes();
 int node;
 TPZVec <double> node_coord(3);
 TPZVec <double> node_coord2(3);

 //escrever os elementos;
 
 string nome ;
  
 int el,n;
 int nelem = Gmesh->ElementVec().NElements();
 int meshindex;

 
 //string filename =  refPatt->GetFileRefPatt() ;
 //ofstream arquivo(filename.c_str());


 ofstream arquivo("prismunif2.dxf");
 arquivo << "  0" << endl ;
 arquivo << "SECTION" << endl ;
 arquivo << "  2" << endl;
 arquivo << "ENTITIES" << endl;


 for (el=0;el<nelem;el++){
    
   //ofstream arquivo(nome.c_str());
 
   //In�cio da declara��o de polyfacemesh

   meshindex=el;
   polyfacemesh (arquivo,meshindex);
  
   //Fim da declara��o de polyfacemesh
   
   TPZGeoEl *gel = Gmesh->ElementVec()[el];
   if (!gel) continue;
   
   // Obtendo Center Point

   int side = gel->NSides()-1;
   TPZVec <REAL> centerER (gel->Dimension(),0.);
   gel->CenterPoint(side,centerER);
   TPZVec <REAL> center (3,0.);
   gel->X(centerER,center);

   // Fim do trabalho com Center Point

   int ncorner = gel->NCornerNodes();

    arquivo << "  71" << endl;
    arquivo << ncorner << endl;

   // Identificando o n�mero de faces do elemento.
  
    int sidedim, faces=0;
    int nsides = gel->NSides();
    
    for (int controlsides=0 ; controlsides<nsides ; controlsides++){
      sidedim = gel->SideDimension(controlsides);
      if (sidedim==2){
	faces++;
      }
    }

   // Fim da identifica��o do n�mero de faces.

    arquivo << "  72" << endl;
    arquivo << faces << endl;
   
   for (n=0;n<ncorner;n++){
     arquivo << "  0" << endl;
     arquivo << "VERTEX" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbEntity" << endl;
     arquivo << "  8" << endl;
     arquivo << "0" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbVertex" << endl;
     arquivo << "  100" << endl;
     arquivo << "AcDbPolyFaceMeshVertex" << endl;
     
     
     TPZGeoNode *node = gel->NodePtr(n);
     
     TPZGeoNode *geonode = &(Gmesh->NodeVec()[node->Id()]);


      // Modificacao de coordenadas
   
     double modulo_vec ;
     if (el==0){

       node_coord[0] = geonode->Coord(0) ;
       node_coord[1] = geonode->Coord(1) ;
       node_coord[2] = geonode->Coord(2) ;

     }
     else {

       node_coord[0] = ( ( (geonode->Coord(0))-center[0] ) /2 ) + center[0]  ;
       node_coord[1] = ( ( (geonode->Coord(1))-center[1] ) /2 ) + center[1]  ;
       node_coord[2] = ( ( (geonode->Coord(2))-center[2] ) /2 ) + center[2]  ;
       


     }




       // Fim da modificacao de coordenadas  


     int i,k;
     for (int i=0; i<3 ; i++){
       //node_coord[i] = geonode->Coord(i);
        k=(i*10)+10;
        arquivo << k << endl;
        arquivo << node_coord[i] << endl;        
     }
     arquivo << "  70" << endl;
     arquivo << "192" << endl;



     
   }



   // Constru��o de faces

   facemaker (arquivo , faces , gel , el , Gmesh);

   // Fim da constru��o de faces


   arquivo << "  0" << endl;
   arquivo << "SEQEND" << endl;

 }








 // Inicio da insercao de indice dos pontos





   int f=0;

   int ncorners = Gmesh->NNodes();

   for (int control2=0;control2<ncorners;control2++){

    arquivo << "  0" << endl;
    arquivo << "TEXT" << endl;
    arquivo << "  100"  << endl;
    arquivo << "AcDbText" << endl;
    arquivo << "  8" << endl;
    arquivo << "Indices" << endl;
    arquivo << "  62" << endl;
    arquivo << "5" << endl;
   
   
    TPZGeoNode *no = &(Gmesh->NodeVec()[control2]);


   int z,k;
     for (int z=0; z<3 ; z++){
       node_coord[z] = no->Coord(z);
       k=(z*10)+10;
       arquivo << k << endl;
       arquivo << node_coord[z] << endl;
    }

   
   arquivo << "  40" << endl;
   arquivo << "0.1" << endl;
   arquivo << "  1" << endl;
   arquivo << f << endl;
   f++ ;
  }


   //ncorners = Gmesh->NNodes();
   int numerel = Gmesh->ElementVec().NElements();
   int lay = 0 , colo = 0 ;

   for (int y=1;y<numerel;y++){

     int w=0;
     lay++ ;
     colo++ ;

     TPZGeoEl *gel = Gmesh->ElementVec()[y];
   
     int numNodes = gel->NNodes();
   
     for (int control=0;control<numNodes;control++){

     

      arquivo << "  0" << endl;
      arquivo << "TEXT" << endl;
      arquivo << "  100"  << endl;
      arquivo << "AcDbText" << endl;
      arquivo << "  8" << endl;
      arquivo << lay << endl;
      arquivo << "  62" << endl;
      arquivo << colo << endl;
      
   
   



    TPZGeoNode *node = gel->NodePtr(control);
     
    TPZGeoNode *geonode = &(Gmesh->NodeVec()[node->Id()]);




   // Obtendo Center Point

   int side = gel->NSides()-1;
   TPZVec <REAL> centerER (gel->Dimension(),0.);
   gel->CenterPoint(side,centerER);
   TPZVec <REAL> center (3,0.);
   gel->X(centerER,center);

   // Fim do trabalho com Center Point
   

     if (el==0){

       node_coord2[0] = geonode->Coord(0) ;
       node_coord2[1] = geonode->Coord(1) ;
       node_coord2[2] = geonode->Coord(2) ;

     }
     else {

       node_coord2[0] = ( ( (geonode->Coord(0))-center[0] ) /2 ) + center[0]  ;
       node_coord2[1] = ( ( (geonode->Coord(1))-center[1] ) /2 ) + center[1]  ;
       node_coord2[2] = ( ( (geonode->Coord(2))-center[2] ) /2 ) + center[2]  ;
       


     }



    int z,k;
     for (int z=0; z<3 ; z++){
       //node_coord[z] = geonode->Coord(z);
       k=(z*10)+10;
       arquivo << k << endl;
       arquivo << node_coord2[z] << endl;
    }

   
   arquivo << "  40" << endl;
   arquivo << "0.1" << endl;
   arquivo << "  1" << endl;
   arquivo << w << endl;
   w++ ;
  }


}
  

 
   arquivo << "  0" << endl;
   arquivo << "ENDSEC" << endl;
   arquivo << "  0" << endl;
   arquivo << "EOF" << endl;

}


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






void siderecog (TPZGeoMesh *Gmesh){

  // Para reconhecer o lado dividido ao meio, basta trabalharmos com o primeiro elemento da malha (pai).
  TPZGeoEl *elemento = Gmesh->ElementVec()[0];
  // Acima � gerado o GelEl do elemento pai .
  
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



  // No loop abaixo testamos todos os sides, e naqueles que s�o veificados como sendo arestas, fazemos o teste de alinhamento dos pontos .

  for (int controle=0;controle<nsides;controle++){
    dimensao = elemento->SideDimension(controle);
    if (dimensao==1){

    lado++;
    no1index = elemento->SideNodeIndex(controle,0);
    no2index = elemento->SideNodeIndex(controle,1);
    
   
    TPZGeoNode *no1 = elemento->NodePtr(no1index);
    TPZGeoNode *no2 = elemento->NodePtr(no2index);
    for (int i=0;i<3;i++){
      no1coord[i] = no1->Coord(i);
      no2coord[i] = no2->Coord(i);
      nomediocoord[i] = nomedio->Coord(i);
    }

    // Calcularemos o determinante da matriz contendo as coordenadas de 
    // tr�s pontos para checar alinhamento (os dois extremos do lado e o n� m�dio)
    
    detcheck = ((nomediocoord[1]*no1coord[2])+(nomediocoord[2]*no2coord[1])+(no1coord[1]*no2coord[2])-(no1coord[2]*no2coord[1])-(nomediocoord[2]*no1coord[1])-(nomediocoord[1]*no2coord[2]))*((nomediocoord[1]*no1coord[2])+(nomediocoord[2]*no2coord[1])+(no1coord[1]*no2coord[2])-(no1coord[2]*no2coord[1])-(nomediocoord[2]*no1coord[1])-(nomediocoord[1]*no2coord[2]))+((nomediocoord[2]*no1coord[0])+(nomediocoord[0]*no2coord[2])+(no1coord[2]*no2coord[0])-(no1coord[0]*no2coord[2])-(nomediocoord[0]*no1coord[2])-(nomediocoord[2]*no2coord[0]))*((nomediocoord[2]*no1coord[0])+(nomediocoord[0]*no2coord[2])+(no1coord[2]*no2coord[0])-(no1coord[0]*no2coord[2])-(nomediocoord[0]*no1coord[2])-(nomediocoord[2]*no2coord[0]))+((nomediocoord[0]*no1coord[1])+(nomediocoord[1]*no2coord[0])+(no1coord[0]*no2coord[1])-(no1coord[1]*no2coord[0])-(nomediocoord[1]*no1coord[0])-(nomediocoord[0]*no2coord[1]))*((nomediocoord[0]*no1coord[1])+(nomediocoord[1]*no2coord[0])+(no1coord[0]*no2coord[1])-(no1coord[1]*no2coord[0])-(nomediocoord[1]*no1coord[0])-(nomediocoord[0]*no2coord[1]));
    




    // Caso o alinhamento seja verificado, guardamos algumas refer�ncias abaixo.
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
  cout << "Estes lados s�o os sides de �ndice :" << endl;
    for (int p=0;p<considesnumber;p++){
      cout << neighbourindex[p] << endl;

    }
  
  
  

}





int main ()
{

  TPZGeoMesh ownermesh;
  std::string filename("/home/ic/luis/Documents/ic/catalogacao/meuRefPattern/Prisma_Triangular_1");
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Triang_Unif.rpt");
  TPZRefPattern *reftetra = new TPZRefPattern (&ownermesh,filename);
  //TPZRefPattern *reftetra = new TPZRefPattern ("/home/pos/cesar/RefPattern/Piram_Rib_Side_7.rpt");
  TPZGeoMesh *Gmesh = reftetra->Mesh();
  dxfdraw(Gmesh);  
  dxfdrawsep(Gmesh, reftetra);

  //Teste da fun��o de reconhecimento de lados

  siderecog(Gmesh);
  //Fim do teste

  
  return 0 ;
}



