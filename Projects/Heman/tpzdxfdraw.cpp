//
// C++ Implementation: tpzdxfdraw
//
// Description: 
//
//
// Author: Luís Guilherme Mello Décourt <gdecourt@gmail.com>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzdxfdraw.h"
#include "pzgmesh.h"
#include "pzgeoel.h"

using namespace std;

namespace pz_dxf {

TPZDXFDraw::TPZDXFDraw(std::string &file,TPZGeoMesh *gmesh) : fDXFFile(file.c_str())
{
  fMesh = gmesh;
}


TPZDXFDraw::~TPZDXFDraw()
{
}


void TPZDXFDraw::PolyFaceMesh (int meshindex){
    fDXFFile << "  0"      << std::endl;
    fDXFFile << "POLYLINE" << std::endl;
    fDXFFile << "  100"    << std::endl;
    fDXFFile << "AcDbEntity" << std::endl;
    fDXFFile << "  8"      << std::endl;
    fDXFFile << "Mesh"     << meshindex << std::endl;
    fDXFFile << "  100"    << std::endl;
    fDXFFile << "AcDbPolyFaceMesh" << std::endl;
    fDXFFile << "  66"     << std::endl;
    fDXFFile << "1"        << std::endl;
    fDXFFile << "  10"     << std::endl;
    fDXFFile << "0.0"      << std::endl;
    fDXFFile << "  20"     << std::endl;
    fDXFFile << "0.0"      << std::endl;
    fDXFFile << "  30"     << std::endl;
    fDXFFile << "0.0"      << std::endl;
    fDXFFile << "  70"     << std::endl;
    fDXFFile << "64"       << std::endl;
}


void TPZDXFDraw::FaceMaker (int faces ,TPZGeoEl *gel , int el){
  int dxfnodeindex, sidenodes , nodeindex , finalnodeindex  ;
  int nsides = gel->NSides();
  int sidedim ;
 
  for (int controlsides=0 ; controlsides<nsides ; controlsides++){
    sidedim = gel->SideDimension(controlsides);
      
    if (sidedim==2){
      dxfnodeindex = 70 ;
  
      fDXFFile << "  0" << std::endl;
      fDXFFile << "VERTEX" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbEntity" << std::endl;
      fDXFFile << "  8" << std::endl;
      fDXFFile << "Faces"  << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbFaceRecord" << std::endl;
      fDXFFile << "  10" << std::endl;
      fDXFFile << "0.0" << std::endl;
      fDXFFile << "  20" << std::endl;
      fDXFFile << "0.0" << std::endl;
      fDXFFile << "  30" << std::endl;
      fDXFFile << "0.0" << std::endl;
      fDXFFile << "  70" << std::endl;
      fDXFFile << "128" << std::endl;
  
      //int numerodenos = gel->NNodes() ;
      sidenodes = gel->NSideNodes(controlsides) ; 
      
  
      for (int c = 0 ; c<sidenodes ; c++){
      //for (int c = 0 ; c<numerodenos ; c++){  
        dxfnodeindex++ ; 
        nodeindex = gel->SideNodeLocIndex(controlsides,c);
        //nodeindex = gel->NodeIndex(c);
        finalnodeindex = nodeindex + 1 ;
      
        fDXFFile << "  " << dxfnodeindex << std::endl;
        fDXFFile << finalnodeindex << std::endl ;
      }
    }
  }
}


void TPZDXFDraw:: DXFDraw (){
  //int nnodes = fMesh->NNodes();
  //int node;
  TPZVec <double> node_coord(3);
  
  string nome ;
  
  int meshindex = 0 ;
    
  int el,n;
  int nelem = fMesh->ElementVec().NElements();
  
  fDXFFile << "  0" << std::endl ;
  fDXFFile << "SECTION" << std::endl ;
  fDXFFile << "  2" << std::endl;
  fDXFFile << "ENTITIES" << std::endl;
  
  for (el=0;el<nelem;el++){
    //ofstream fDXFFile(nome.c_str());
    meshindex = el ; 
    PolyFaceMesh (meshindex);
    TPZGeoEl *gel = fMesh->ElementVec()[el];
    if (!gel) continue;
    int ncorner = gel->NCornerNodes();
    fDXFFile << "  71" << std::endl;
    fDXFFile << ncorner << std::endl;
    int faces = 0 , nsides , controlsides = 0 , sidedim ;
  
    // Identificando o nmero de faces do elemento.
  
    nsides = gel->NSides();
    for (controlsides=0 ; controlsides<nsides ; controlsides++){
      sidedim = gel->SideDimension(controlsides);
      if (sidedim==2){
        faces++;
      }
    }

    // Fim da identifica�o do nmero de faces.
    fDXFFile << "  72" << std::endl;
    fDXFFile << faces << std::endl;
  
    for (n=0;n<ncorner;n++){
      fDXFFile << "  0" << std::endl;
      fDXFFile << "VERTEX" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbEntity" << std::endl;
      fDXFFile << "  8" << std::endl;
      fDXFFile << "Vertices" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbVertex" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbPolyFaceMeshVertex" << std::endl;
      
      TPZGeoNode *node = gel->NodePtr(n);
      
      TPZGeoNode *geonode = &(fMesh->NodeVec()[node->Id()]);
      int i,k;
  
      for (i=0; i<3 ; i++){
        node_coord[i] = geonode->Coord(i);
        k = (i*10)+10 ;       
        fDXFFile << k << std::endl;
        fDXFFile << node_coord[i] << std::endl;        
      }
      fDXFFile << "  70" << std::endl;
      fDXFFile << "192" << std::endl;
    }
  
    // Defini�o do tipo de elemento e chamada de fun�o construtora de faces
    FaceMaker (faces , gel , el);
  
    // Fim da constru�o de faces
    fDXFFile << "  0" << std::endl;
    fDXFFile << "SEQEND" << std::endl;
  }

  // Inicio da insercao de indice dos pontos
  
  int f=0;
  
  //int numerel = fMesh->ElementVec().NElements();
  
  //for (int control=0;control<numerel;control++){
  
  // TPZGeoEl *gel = fMesh->ElementVec()[control];
  //if (!gel) continue;
  //int ncorner = gel->NCornerNodes();
  
  int ncorners = fMesh->NNodes();

  for (int control2=0;control2<ncorners;control2++){

    fDXFFile << "  0" << std::endl;
    fDXFFile << "TEXT" << std::endl;
    fDXFFile << "  100"  << std::endl;
    fDXFFile << "AcDbText" << std::endl;
    fDXFFile << "  8" << std::endl;
    fDXFFile << "Indices" << std::endl;
    fDXFFile << "  62" << std::endl;
    fDXFFile << "5" << std::endl;
  
  
    TPZGeoNode *no = &(fMesh->NodeVec()[control2]);

    //TPZGeoNode *node = gel->NodePtr(control2);
    
    //TPZGeoNode *geonode = &(fMesh->NodeVec()[node->Id()]);
    int z,k;
    for (z=0; z<3 ; z++){
      k=z;
      node_coord[z] = no->Coord(z);
      k=(z*10)+10;
      fDXFFile << k << std::endl;
      fDXFFile << node_coord[z] << std::endl;
    }
    fDXFFile << "  40" << std::endl;
    fDXFFile << "0.1" << std::endl;
    fDXFFile << "  1" << std::endl;
    fDXFFile << f << std::endl;
    f++ ;
  }
   
  //ncorners = fMesh->NNodes();
  int numerel = fMesh->ElementVec().NElements();
  int lay = 0 , colo = 0 ;

  for (int y=0;y<numerel;y++){
    int w=0;
    lay++ ;
    colo++ ;

    TPZGeoEl *gel = fMesh->ElementVec()[y];
  
    int numNodes = gel->NNodes();
  
    for (int control=0;control<numNodes;control++){
      fDXFFile << "  0" << std::endl;
      fDXFFile << "TEXT" << std::endl;
      fDXFFile << "  100"  << std::endl;
      fDXFFile << "AcDbText" << std::endl;
      fDXFFile << "  8" << std::endl;
      fDXFFile << lay << std::endl;
      fDXFFile << "  62" << std::endl;
      fDXFFile << colo << std::endl;
   
      //TPZGeoNode *no = &(fMesh->NodeVec()[control]);

      //TPZGeoEl *gel = fMesh->ElementVec()[control];
      //TPZGeoNode *node = gel->NodePtr(control);
      //TPZGeoNode *geonode = &(fMesh->NodeVec()[control]);
      
      
      TPZGeoNode *node = gel->NodePtr(control);
        
      TPZGeoNode *geonode = &(fMesh->NodeVec()[node->Id()]);
      
      int z,k;
      for (z=0; z<3 ; z++){
          node_coord[z] = geonode->Coord(z);
          k=(z*10)+10;
          fDXFFile << k << std::endl;
          fDXFFile << node_coord[z] << std::endl;
      }
      
      
      fDXFFile << "  40" << std::endl;
      fDXFFile << "0.1" << std::endl;
      fDXFFile << "  1" << std::endl;
      fDXFFile << w << std::endl;
      w++ ;
    }
  }
  fDXFFile << "  0" << std::endl;
  fDXFFile << "ENDSEC" << std::endl;
  fDXFFile << "  0" << std::endl;
  fDXFFile << "EOF" << std::endl;
}


/* Inicio da gera�o do arquivo DXF que gera imagens 
   das componentes filhas reduzidas */

void TPZDXFDraw::DXFDrawSep (){
//  int nnodes = fMesh->NNodes();
//  int node;
  TPZVec <double> node_coord(3);
  TPZVec <double> node_coord2(3);
  
  //escrever os elementos;
  
  string nome ;
    
  int el,n;
  int nelem = fMesh->ElementVec().NElements();
  int meshindex;
  
  //string filename =  refPatt->GetFileRefPatt() ;
  //ofstream fDXFFile(filename.c_str());
  fDXFFile << "  0" << std::endl ;
  fDXFFile << "SECTION" << std::endl ;
  fDXFFile << "  2" << std::endl;
  fDXFFile << "ENTITIES" << std::endl;
  
  
  for (el=0;el<nelem;el++){
      
    //ofstream fDXFFile(nome.c_str());
  
    //In�io da declara�o de polyfacemesh
  
    meshindex=el;
    PolyFaceMesh (meshindex);
    
    //Fim da declara�o de polyfacemesh
    
    TPZGeoEl *gel = fMesh->ElementVec()[el];
    if (!gel) continue;
    
    // Obtendo Center Point
  
    int side = gel->NSides()-1;
    TPZVec <REAL> centerER (gel->Dimension(),0.);
    gel->CenterPoint(side,centerER);
    TPZVec <REAL> center (3,0.);
    gel->X(centerER,center);
  
    // Fim do trabalho com Center Point
  
    int ncorner = gel->NCornerNodes();
  
    fDXFFile << "  71" << std::endl;
    fDXFFile << ncorner << std::endl;
  
    // Identificando o nmero de faces do elemento.
    
    int sidedim, faces=0;
    int nsides = gel->NSides();
      
    for (int controlsides=0 ; controlsides<nsides ; controlsides++){
      sidedim = gel->SideDimension(controlsides);
      if (sidedim==2){
        faces++;
      }
    }
  
    // Fim da identifica�o do nmero de faces.
  
    fDXFFile << "  72" << std::endl;
    fDXFFile << faces << std::endl;
  
    for (n=0;n<ncorner;n++){
      fDXFFile << "  0" << std::endl;
      fDXFFile << "VERTEX" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbEntity" << std::endl;
      fDXFFile << "  8" << std::endl;
      fDXFFile << "0" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbVertex" << std::endl;
      fDXFFile << "  100" << std::endl;
      fDXFFile << "AcDbPolyFaceMeshVertex" << std::endl;
      
      
      TPZGeoNode *node = gel->NodePtr(n);
      
      TPZGeoNode *geonode = &(fMesh->NodeVec()[node->Id()]);
        
      // Modificacao de coordenadas
    
      //double modulo_vec ;
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
      for (i=0; i<3 ; i++){
        //node_coord[i] = geonode->Coord(i);
        k=(i*10)+10;
        fDXFFile << k << std::endl;
        fDXFFile << node_coord[i] << std::endl;        
      }
      fDXFFile << "  70" << std::endl;
      fDXFFile << "192" << std::endl;
    }
    // Constru�o de faces
    FaceMaker (faces , gel , el);
    // Fim da constru�o de faces
    fDXFFile << "  0" << std::endl;
    fDXFFile << "SEQEND" << std::endl;
  }
  // Inicio da insercao de indice dos pontos
  int f=0;

  int ncorners = fMesh->NNodes();

  for (int control2=0;control2<ncorners;control2++){

    fDXFFile << "  0" << std::endl;
    fDXFFile << "TEXT" << std::endl;
    fDXFFile << "  100"  << std::endl;
    fDXFFile << "AcDbText" << std::endl;
    fDXFFile << "  8" << std::endl;
    fDXFFile << "Indices" << std::endl;
    fDXFFile << "  62" << std::endl;
    fDXFFile << "5" << std::endl;

    TPZGeoNode *no = &(fMesh->NodeVec()[control2]);

    int z,k;
    for (z=0; z<3 ; z++){
      node_coord[z] = no->Coord(z);
      k=(z*10)+10;
      fDXFFile << k << std::endl;
      fDXFFile << node_coord[z] << std::endl;
    }
    
    fDXFFile << "  40" << std::endl;
    fDXFFile << "0.1" << std::endl;
    fDXFFile << "  1" << std::endl;
    fDXFFile << f << std::endl;
    f++ ;
  }

  //ncorners = fMesh->NNodes();
  int numerel = fMesh->ElementVec().NElements();
  int lay = 0 , colo = 0 ;

  for (int y=1;y<numerel;y++){
    int w=0;
    lay++ ;
    colo++ ;

    TPZGeoEl *gel = fMesh->ElementVec()[y];
  
    int numNodes = gel->NNodes();
  
    for (int control=0;control<numNodes;control++){
      fDXFFile << "  0" << std::endl;
      fDXFFile << "TEXT" << std::endl;
      fDXFFile << "  100"  << std::endl;
      fDXFFile << "AcDbText" << std::endl;
      fDXFFile << "  8" << std::endl;
      fDXFFile << lay << std::endl;
      fDXFFile << "  62" << std::endl;
      fDXFFile << colo << std::endl;
      TPZGeoNode *node = gel->NodePtr(control);
      TPZGeoNode *geonode = &(fMesh->NodeVec()[node->Id()]);
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
      for (z=0; z<3 ; z++){
        //node_coord[z] = geonode->Coord(z);
        k=(z*10)+10;
        fDXFFile << k << std::endl;
        fDXFFile << node_coord2[z] << std::endl;
      }
      fDXFFile << "  40" << std::endl;
      fDXFFile << "0.1" << std::endl;
      fDXFFile << "  1" << std::endl;
      fDXFFile << w << std::endl;
      w++ ;
    }
  }
  fDXFFile << "  0" << std::endl;
  fDXFFile << "ENDSEC" << std::endl;
  fDXFFile << "  0" << std::endl;
  fDXFFile << "EOF" << std::endl;
}


};
