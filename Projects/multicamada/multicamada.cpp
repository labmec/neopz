
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzplaca.h"
#include "pzintel.h"

#include "TPZGeoElement.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzelgc3d.h"
#include "pzbndcond.h"
#include "pztempmat.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelcc3d.h"
#include "pzelgt3d.h"
#include "pzelct3d.h"
#include "pzelgpi3d.h"
#include "pzelcpi3d.h"
#include "pzelgpr3d.h"
#include "pzelcpr3d.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
#include "TPZConservationLaw.h"
#include "TPZConsLawTest.h"
#include "TPZEulerConsLaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
#include "pzreal.h"
#include "pzdxmesh.h"
#include "pzmatorthotropic.h"
#include "TPZMulticamadaOrtho.h"
#include "TPZPlacaOrthotropic.h"
//#include "TPZSurfaceInteg.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <string.h>
#include <math.h>
using namespace std;

void ProcessamentoLocal(TPZGeoMesh &gmesh,ostream &out);
void CriaNos(int num, TPZGeoMesh &geomesh, double list [20][3]);
void CriaElementos(int numelem,int ncon,TPZGeoMesh &geomesh, int list[20][8] );
void CriaCondCont(TPZGeoMesh &gmesh);
void ExecutaAnalysis(TPZCompMesh &cmesh,TPZMaterial *mat);
void PosProcessamento(TPZAnalysis &an);
void Print(TPZGeoMesh *geomesh);
void CriaCondCont2(TPZGeoMesh &gmesh);
//static double NosPlaca[8][3] = {{




/**
 * -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <-
 */

int main(){
  int ord, quantplacas;
  char cout;
  //double hight, dx, dy;


  //chamando construtor da classe TPZMulticamadaOrtho
  //TPZMulticamadaOrthotropic *multcam = new TPZMulticamadaOrthotropic(*geomesh, *compmesh);
// mult (0.2, 1.5, 1.8, gmesh, compmesh);
//  TPZPlacaOrthotropic *placa = new TPZPlacaOrthotropic;
//  placa (0.2);

//criando objetos das classes computacional e geométrica
  TPZGeoMesh *geomesh = new TPZGeoMesh;
  TPZCompMesh *compmesh = new TPZCompMesh(geomesh);
  //  TPZMaterial *mat = new TPZMaterial(1);
  TPZMulticamadaOrthotropic *multcam = new TPZMulticamadaOrthotropic(2.3, 2.0, 1.5, geomesh, compmesh);
                        //  cout << "Entre ordem : 1,2,3,4,5 : -> ";
                        //   cin >> ord;
  //  ord = 1; cout << endl;
                        //TPZCompEl::gOrder = ord;
  //TPZMulticamadaOrthotropic multcam2(2.3, 1.2, 1.5, geomesh, compmesh);
 //TPZInterpolatedElement *cel = new TPZInterpolatedElement(*geomesh, TPZGeoEl *reference, int &index); 
 //TPZPlacaOrthotropic *placa = new TPZPlacaOrthotropic(cel); //, 1.2, 1.5);
                        //  geomesh->NodeVec().Resize(12);

  multcam->RQPlacas();
 TPZVec<REAL> co(3,0.);//zerando as 3 coordenadas
 geomesh->NodeVec().Resize(12);
 geomesh->NodeVec()[0].Initialize(co,*geomesh);//primeiro elemento[0] recebe as coord (0 0 0)
 co[0] = 1.;//atribui o valor 1 ao primeiro coord co[0],(co[1], co[2] não alteram o valor)
 geomesh->NodeVec()[1].Initialize(co,*geomesh);//o elemento [1] da lista recebe as coord(1 0 0)sendo que co[0] foi alterada
 co[1] = 1.;//alterando valor da posicão co[1]para passar ao elemento NodeVec[2] (junto com as coordenadas co[0], co[2])     
 geomesh->NodeVec()[2].Initialize(co,*geomesh);//1 1 0
 co[0] = 0.;
 geomesh->NodeVec()[3].Initialize(co,*geomesh);//0 1 0
 co[1] = 0.;
 co[2] = .2;
 geomesh->NodeVec()[4].Initialize(co,*geomesh);//0 0 .2
 co[0] = 1.;
 geomesh->NodeVec()[5].Initialize(co,*geomesh);//1 0 .2
 co[1] = 1.;
 geomesh->NodeVec()[6].Initialize(co,*geomesh);//1 1 .2
 co[0] = 0.;
 geomesh->NodeVec()[7].Initialize(co,*geomesh);//0 1 .2
 TPZVec<int> index(8);
 int i, indice;
 for (i=0; i<8; i++){
   index[i] = i;
 }
 // geomesh->Print();
 
 TPZGeoEl*  geo1 =  geomesh->CreateGeoElement(ECube, index, 0, indice);//criando elemento geo, pois os 8 nós da placa estão criados

 // TPZPlacaOrthotropic placa(*cel);//colocar indice para placa ???????????????????????
 co[1] = 0.;
 co[2] = .3;
 geomesh->NodeVec()[8].Initialize(co,*geomesh);//0 0 .3
 co[0] = 1.;
 geomesh->NodeVec()[9].Initialize(co,*geomesh);//1 0 .3
 co[1] = 1.;
 geomesh->NodeVec()[10].Initialize(co,*geomesh);//1 1 .3
 co[0] = 0.;
 geomesh->NodeVec()[11].Initialize(co,*geomesh);//0 1 .3
 for (i=0; i<8; i++){
   index[i] = i+4;
 }
 TPZGeoEl*  geo2 = geomesh->CreateGeoElement(ECube, index, 1, indice);

 //geomesh->Print();

 //chamar outra placa para inclui-la na multicamada

  TPZFMatrix naxes(3,3);
  //  REAL E = 200000.0, G = 173913.0, v = 0.15;
  //REAL E = 10000.0, G = 4347.83, v = 0.15;//teste
  //REAL E = 15.133e6, G = 2.59398e7, v = 0.471;//Pinus Pinaster ait.
 //  REAL E = 69.0e6, G = 5.14378e6, v = 0.330;//Aluminio
//   REAL eppx=E, eppy=E, eppz=E, vxy=v, vyz=v,vzx=v, gxy=G, gyz=G, gzx=G;
//  int matindex = 1;

  naxes(0,0) =  1.0;    naxes(0,1) =  0.0;    naxes(0,2) =  0.0;
  naxes(1,0) =  0.0;    naxes(1,1) =  1.0;    naxes(1,2) =  0.0;
  naxes(2,0) =  0.0;    naxes(2,1) =  0.0;    naxes(2,2) =  1.0;

  //criando objeto para a classe material
  TPZMaterial *orto;
   orto = new TPZMatOrthotropic(1, naxes, 69.e06, 69.e06, 69.e06, 0.33, 0.33, 0.33, 5.15378e06, 5.15378e06, 5.15378e06) ;
  compmesh->InsertMaterialObject(orto);  //inserindo o material na classe computacional
 TPZMaterial *orto2 = new TPZMatOrthotropic(2, naxes, 50, 50, 69.e06, 0.33, 0.33, 0.33, 5.15378e06, 5.15378e06, 5.15378e06) ;
  compmesh->InsertMaterialObject(orto2);  //inserindo o material na classe computacional

// montagem da conectividade 
 geomesh->BuildConnectivity2();
 
 geomesh->Print();
 
 

/**
 * dica: devera ler as carasteristica de cada placa a partir
 * de um arquivo de entrada e criar o material ortotropico
 * com aqueles dados (classe TPZMatOrthotropic), logo chamar o
 * método AddPlacaOrth() para adicionar esta placa na 
 * lista de placas, variável fPlacaOrth
 */ 
   
  //cria condi¢ões de contorno
 //  CriaCondCont(*geomesh);
 CriaCondCont2(*geomesh);


  //cria elementos computacionais
  compmesh->AutoBuild();

 // criando objeto multicamada

 TPZMulticamadaOrthotropic multcam2(0.2, 1.7, 2.0, geomesh, compmesh); // passando altura, dx e dy

 TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(geo1->Reference());
 TPZPlacaOrthotropic placa1 (intel);
 multcam->AddPlacaOrtho(&placa1);// usa-se & ao invés de * porque placa não foi criado com new
 // multcam->Print();


 //multcam.AddPlacaOrtho(*placa1);
 intel = dynamic_cast<TPZInterpolatedElement *>(geo2->Reference());
 TPZPlacaOrthotropic placa2 (intel);
 multcam->AddPlacaOrtho(&placa2);
 // criando o objeto para o  material da placa
 //orto = new TPZMatOrthotropic(1, naxes, 69.e06, 69.e06, 69.e06, 0.33, 0.33, 0.33, 5.15378e06, 5.15378e06, 5.15378e06) ;
 //multcam.AddPlacaOrtho(*placa2);
 multcam->Print();


 // TPZVec <elastortot> fPlacaOrth(0.2,orto);
  
  //ajusta elementos no contorno
  compmesh->AdjustBoundaryElements();

  //analysis do problema
 //  ExecutaAnalysis(*compmesh,orto);

  //arquivo de saida de dados
  ofstream data("mesh.out");
  geomesh->Print(data);
  compmesh->Print(data);
  data.flush();
  data.close();

  delete compmesh;
  delete geomesh;
  return 0;
}


void CriaCondCont2(TPZGeoMesh &gmesh){
  //as CC de faces locais não estão contidas nas CC -3
  //essa diferença é necessária para os elementos global-local
  //cout << "main::CriaCondContLocal CC -1 e -2 ja devem existir\n";
  //StopTime();

  int indicematerial = 1;
  TPZMaterial *orthotropic = gmesh.Reference()->FindMaterial(indicematerial);
  if(!orthotropic){
    cout << "main::CriaCondCont5El material nao existe, CC nao criadas\n";
    cout << "\t\tindice material pedido : " << indicematerial << endl;
    return;
  }

  //malha computacional
  TPZCompMesh *cmesh = gmesh.Reference();
  //valor das CC
  TPZFMatrix val1(3,3,0.),val2(3,1,0.);
  //BIG number
  REAL big = 1.e12;
  TPZBndCond *bc;
  val1(1,1) = big;
  val2(1,0) = big * (1.0);
  
  //orthotropic->TPZMaterial(-1);
  bc = orthotropic->CreateBC(-1,2,val1,val2);
  cmesh->InsertMaterialObject(bc);

  val1(1,1) = big;
  val2(1,0) = big * (-1.0);
  //orthotropic->TPZMaterial(-2);
  bc = orthotropic->CreateBC(-2,2,val1,val2);//bc->SetForcingFunction(BCSolution);
  cmesh->InsertMaterialObject(bc);
 
  
  TPZGeoEl *geo1 = gmesh.ElementVec()[0];
  
  
  TPZGeoElBC(geo1,22,-1,gmesh);//o local está contido no elemento central

  TPZGeoElBC(geo1,24,-2,gmesh);
  
  TPZGeoEl *geo2 = gmesh.ElementVec()[1];
  
  
  TPZGeoElBC(geo2,22,-1,gmesh);//o local está contido no elemento central

  TPZGeoElBC(geo2,24,-2,gmesh);

}
