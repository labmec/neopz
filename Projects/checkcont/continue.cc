// -*- c++ -*-
/** 
    Este programa testa as possiveis combinacões da definicão dos elementos por
    meio dos seus id's globais em uma malha formada por um elemento simples com 
    arestas e faces superpostas. Os elementos testados são:
    triangulo, quadrilatero, tetraedro, piramide, prisma e hexaedro.
    O objetivo é testar as transformacões que podem se dar entre os diferentes 
    lados dos elementos. As ordens de interpolacão dos lados dos elementos desta malha
    são diferentes entre sim, para diferentes lados. Divisões não são efetuadas.
    As ordens para os diferentes lados do elemento devem ser entradas por arquivos

    Nota: basta utilizar orden uniforme para os lados do elemento todo pois a ordem 
    de cada lado do elemento sempre deve ser conforme com a ordem do lado do elemento vizinho
*/
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgnode.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelgc3d.h"
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
#include "pzmattest3d.h"
#include "TPZRefCube.h"
#include "pzvec.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;


static int test = 0;
void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f);
void CC(TPZGeoMesh *firstmesh);
TPZMaterialTest3D *mat;
void Beep(int batida);
void AvisoAudioVisual(int nbat);
void GeometricTest(TPZGeoMesh *geomesh,ofstream &out);
void ComputationalTest(TPZCompMesh *cmesh,ofstream &out);
void ComputationalTestSide(TPZCompMesh *cmesh, ofstream &out);
void SideParameterToElement1d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementT2d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementQ2d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementT3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementPi3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementPr3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void SideParameterToElementC3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
void LerMalha(TPZGeoMesh *geomesh,TPZCompMesh *& compmesh);
void AutomaticDivide(TPZCompMesh &cmesh);
void Divide(TPZCompMesh *compmesh);
void NeighbourhoodCycles(TPZGeoMesh &gmesh,ofstream &out);
int FirstSideShape(TPZCompMesh *mesh,TPZInterpolatedElement *cel,TPZConnect *connect);
void SideParameterToElement(int nsides,int side,TPZVec<REAL> &pint);
void ShapeLocId(TPZInterpolatedElement *cel,TPZInterpolatedElement *neighref,int locid,int neighside,int &neighlocshape);
void ComputationalTestInnerSide(TPZCompMesh *cmesh, ofstream &out);
void TesteShapeFunctions(TPZGeoMesh *firstmesh);
void Test_Continuity(TPZGeoMesh *firstmesh,ofstream &outgm);

//static int permuta2[2][2] = {{1,0},{0,1}};// 2! = 2

static int permuta3[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};// 3! = 6

static int permuta4[24][4] = {{0,1,2,3},{0,1,3,2},{0,2,1,3},{0,2,3,1},{0,3,1,2},{0,3,2,1},
			      {1,0,2,3},{1,0,3,2},{1,2,0,3},{1,2,3,0},{1,3,0,2},{1,3,2,0},
			      {2,0,1,3},{2,0,3,1},{2,1,0,3},{2,1,3,0},{2,3,0,1},{2,3,1,0},
			      {3,0,1,2},{3,0,2,1},{3,1,0,2},{3,1,2,0},{3,2,0,1},{3,2,1,0}};// 4! = 24

static double triangulo[3][3] = { {0.,0.,-1.},{2.,0.,0.},{0.,0.,1.} };//0 1 2

static double quadrilatero[4][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.} };//0 1 2 3

static double tetraedro[4][3] = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.} };//0 1 2 3 4

static double prisma[6][3] = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},
		             {0.,0.,1.},{1.,0.,1.},{0.,1.,1.} };//0 1 2 3 4 5

static double piramide[5][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},{0.5,0.5,1.} };//0 1 2 3 4

static double hexaedro[8][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
		             {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.} };//0 1 2 3 4 5 6 7

void ConstroeNodes(int typel,TPZGeoMesh *geomesh);
void IdsDefinitions(int typeel,TPZVec<int> &nodes,int &ntype,int &numel);
void LargestNeighbour(TPZGeoElSide THIS,TPZGeoElSide &largeel, int targetdimension);


int main() {


  //AvisoAudioVisual(5);
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   //Arquivos de saida
   ofstream outgm("MGEOMETRICA.OUT");
   ofstream outcm("MCOMPUTACIONAL.OUT");
   TPZCompMesh *secondmesh;// = new TPZCompMesh(firstmesh);
   //malha geométrica
   LerMalha(firstmesh,secondmesh);
   //CC : condicões de contorno
   //CC(elg,firstmesh,dir);
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
   //ordem de interpolacao
   int ord;
   cout << "Enter order 1,2,3,4,5,... : \n";
   cin >> ord;
   TPZCompEl::gOrder = ord;
   //construção malha computacional
   secondmesh->AutoBuild();
   if(0){
     firstmesh->Print(outgm);
     outgm.flush();
     secondmesh->Print(outcm);
     outcm.flush();
   }

   //PRIMEIRO TESTE: MALHA
   if(0) return 0;

   if(1) Test_Continuity(firstmesh,outgm);

   cout << "\nEnd # \n";
   int sim;
   //cin >> sim;
   delete secondmesh;
   delete firstmesh;
   return 0;
}

//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN

void Test_Continuity(TPZGeoMesh *firstmesh,ofstream &outgm){

   TPZGeoEl *geoel = firstmesh->ElementVec()[0];
   int side,ncorners = geoel->NNodes(),nsides = geoel->NSides();
   if(geoel->Dimension() == 3) nsides--;//tira o volume
   for(side=ncorners;side<nsides;side++){
     TPZGeoElSide geoelside(geoel,side);
     TPZGeoElSide neigh = geoelside.Neighbour();
     while(neigh.Exists() && neigh != geoelside){
       if(neigh.Element()->Dimension() == geoel->SideDimension(side)) break;
       neigh = neigh.Neighbour();
     }
     if(neigh.Element() == geoelside.Element()){
       cout << "\nmain:: neighbourhooh cycle element error, BC element not exists!\n";
       continue;
     }
     TPZGeoEl *bcel = neigh.Element();
     if(geoel->SideDimension(side) == 1) {//arestas
       
       TPZGeoNode node0 = firstmesh->NodeVec()[bcel->NodeIndex(0)];
       int nodeid0 = node0.Id();
       TPZGeoNode node1 = firstmesh->NodeVec()[bcel->NodeIndex(1)];
       int nodeid1 = node1.Id();
       node0.SetNodeId(nodeid1);
       node1.SetNodeId(nodeid0);
       firstmesh->Print(outgm);
       TesteShapeFunctions(firstmesh);
     }
     int i,j;
     if(geoel->SideDimension(side) == 2) {//faces

       TPZVec<TPZGeoNode *> nodes;
       TPZVec<int> nodesid;
       if(geoel->NSides() == 7 || geoel->NSides() == 15 || geoel->SideNodeIndex(side,3) < 0){//face quadrilateral
	 nodes.Resize(3);
	 nodesid.Resize(3);
	 for(i=0;i<3;i++) nodes[i] = &firstmesh->NodeVec()[bcel->NodeIndex(i)];
	 for(j=0;j<6;j++){	   
	   for(i=0;i<3;i++)  nodesid[i] = nodes[i]->Id();
	   for(i=0;i<3;i++)  nodes[i]->SetNodeId(nodesid[permuta3[j][i]]);
	   firstmesh->Print(outgm);
	   TesteShapeFunctions(firstmesh);
	 }
       }
       else {//faces quadrilaterais
	 nodes.Resize(4);
	 nodesid.Resize(4);
	 for(i=0;i<4;i++) nodes[i] = &firstmesh->NodeVec()[bcel->NodeIndex(i)];
	 for(j=0;j<24;j++){	   
	   for(i=0;i<4;i++)  nodesid[i] = nodes[i]->Id();
	   for(i=0;i<4;i++)  nodes[i]->SetNodeId(nodesid[permuta4[j][i]]);
	   firstmesh->Print(outgm);
	   TesteShapeFunctions(firstmesh);
	 }
       }//if
     }//dimension = 2
   }//for sides do elemento 3D
}

void GeometricTest(TPZGeoMesh *geomesh, ofstream &out) {//quando a malha computacional tem elementos de diferentes níveis

   int nelg = geomesh->ElementVec().NElements();
   int el,nsides;
   TPZGeoEl *gel;
   for(el=0;el<nelg;el++) {
      gel = geomesh->ElementVec()[el];
      if(!gel) continue;
      nsides = gel->NSides();
      int ncorners = 0;//gel->NNodes();
      for(int side=ncorners;side<nsides;side++) {
         TPZGeoElSide largel,gelside(gel,side);
         for(int target = gelside.Dimension();target < 3;target++) {
	   LargestNeighbour(gelside,largel,target);
	   if(!largel.Element() || largel.Element()->Id()==gel->Id()) continue;
	   TPZIntPoints *intrule = gel->CreateSideIntegrationRule(side,2);
	   if(!intrule) {
	     cout << "\nRegra de integracao nula\n";
	     continue;
	   }
	   while(gel!=largel.Element()) {
	     REAL weight;
	     int dim = gelside.Dimension();
	     TPZTransform t(dim);
	     gelside.SideTransform3(largel,t);
	     int numint = intrule->NPoints();
	     for(int it=0; it<numint; it++) {
	       TPZVec<REAL> par(3,0.),pointl(3,0.),result1(3,0.),result2(3,0.);
	       intrule->Point(it,par,weight);
	       t.Apply(par,pointl);
	       SideParameterToElement(nsides,side,par);
	       int nsideslargel = largel.Element()->NSides();
	       SideParameterToElement(nsideslargel,largel.Side(),pointl);
	       gel->X(par,result1);
	       largel.Element()->X(pointl,result2);
	       REAL result,error = 0.;;
	       for(int i=0; i<3; i++) {
		 result = result1[i] - result2[i];
		 //error += pow( result , 2.0);//error da 4 com result = 0
		 error += result*result;
	       }
	       if(sqrt(error) > 1.e-10) {
		 cout << " error = " << error << endl;
		 out << " error = " << error << endl;
	       }
	     }
	     out << gel->Id() << "/" << side << " -> ";
	     out << largel.Element()->Id() << "/" << largel.Side() << endl;
	     if(largel.Element()->Level()!=gel->Level()) break;
	     largel = largel.Neighbour();
	   }
	   delete intrule;
         }
      }
   }
}



void ComputationalTest(TPZCompMesh *cmesh, ofstream &out) {

   int nelems = cmesh->ElementVec().NElements();
   int nsides;
   TPZInterpolatedElement *cel;
   TPZGeoEl *gel;
   for(int el=0;el<nelems;el++) {
      cel = (TPZInterpolatedElement *) cmesh->ElementVec()[el];
      if(!cel) continue;
      gel = cel->Reference();
      if(!gel) continue;
      nsides = gel->NSides();
      for(int side=0;side<nsides;side++) {
         TPZGeoEl *neighbour =  gel->Neighbour(side).Element();
         if(!neighbour) continue;
         TPZIntPoints *sideintrule = gel->CreateSideIntegrationRule(side,cel->SideOrder(side)*2);
         if(!sideintrule) {
            cout << "\nRegra de integracao nula\n";
            continue;
         }
         TPZGeoElSide neighcycle(neighbour,gel->Neighbour(side).Side());
         while(gel!=neighcycle.Element()) {//!!!!!!!!!!!!!!!
            neighbour = neighcycle.Element();
            TPZInterpolatedElement *neighref = (TPZInterpolatedElement *) neighbour->Reference();
            if(!neighref) break;
            int nsideshapes = cel->NSideShapeF(side);
            int neighbourside = neighcycle.Side();
            int nneighsideshapes = neighref->NSideShapeF(neighbourside);
            if(nsideshapes != nneighsideshapes) {
               cout << "\nnNumero de shapes incompatibles\n";
                out << "\nnNumero de shapes incompatibles\n";
               continue;
            }
            int dimside = gel->SideDimension(side);
            int neighdimside = neighbour->SideDimension(neighbourside);
            if(dimside!=neighdimside) {

               cout << "\nLados de dimensoes distintas\n";

                out << "\nLados de dimensoes distintas\n";

            }

            int nsideconnect = cel->NSideConnects(side);

            int nneighsideconnet = neighref->NSideConnects(neighbourside);

            if(nsideconnect!=nneighsideconnet) {

               cout << "\nNumero de connects distintos\n";

                out << "\nNumero de connects distintos\n";

            }

            //aqui nsideshape = nneighsideshape e dimside = neighdimside

            int nshape = cel->NShapeF();

            int nneighshape = neighref->NShapeF();

            int dim = cel->Dimension();

            int neighdim = neighref->Dimension();
	    //suficiente para ordem 5 do cubo (o mais grande)
            REAL phistore[220],dphistore[660],phineighstore[220],dphineighstore[660];	
            TPZFMatrix phi(nshape,1,phistore,220),phineigh(nneighshape,1,phineighstore,220);
	    //phiel(nshape,1,0.),phineigh(nneighshape,1,0.)
            TPZFMatrix dphi(dim,nshape,dphistore,660),dphineigh(neighdim,nneighshape,dphineighstore,660);
	    //dphiel(dim,nshape,0.),dphineigh(neighdim,nneighshape,0.)

            //calcula a transformação entre lados vizinhos

            TPZGeoElSide geoside(gel,side),neighside(neighbour,neighbourside);

            TPZTransform t(dimside);

            geoside.SideTransform3(neighside,t);//t = geoside.NeighbourSideTransform(neighside);

            int numintp = sideintrule->NPoints();

            REAL weight,error;

            int locid,cindex,locshape,neighlocshape,nconshp;

            for(int i=0;i<nsideconnect;i++) {

               locid = cel->SideConnectLocId(i,side);//id local do connect i do side

               cindex = cel->ConnectIndex(locid);//index do connect

               TPZConnect *connect = &(cmesh->ConnectVec()[cindex]);

               locshape = FirstSideShape(cmesh,cel,connect);

               neighlocshape = FirstSideShape(cmesh,neighref,connect);

               nconshp = cel->NConnectShapeF(locid);

               for(int ip=0;ip<numintp;ip++) {

                  TPZVec<REAL> point(3,0.),pint(3,0.),pont(3,0.);

                  sideintrule->Point(ip,point,weight);

                  pint = point;

                  SideParameterToElement(gel->NSides(),side,pint);

                  cel->Shape(pint,phi,dphi);

                  t.Apply(point,pont);

                  SideParameterToElement(neighbour->NSides(),neighbourside,pont);

                  neighref->Shape(pont,phineigh,dphineigh);

                  for(int k=0;k<nconshp;k++) {

                     error = fabs(phi(locshape+k,0)-phineigh(neighlocshape+k,0));

                     if(error > 1.e-10) {

                        cout << "\nerro no valor da funcao ";
                        cout << "\nid do elemento " << gel->Id();
                        cout << "\nside do elemento " << side;
                        cout << "\nshape local do elemento "<< (locshape+k);
                        cout << "\nid do vizinho " << neighbour->Id();
                        cout << "\nside do vizinho " << neighbourside;
                        cout << "\nshape local do vizinho " << (neighlocshape+k) << endl;
                         out << "\nerro no valor da funcao ";
                         out << "\nid do elemento " << gel->Id();
                         out << "\nside do elemento " << side;
                         out << "\nshape local do elemento " << (locshape+k);
                         out << "\nid do vizinho " << neighbour->Id();
                         out << "\nside do vizinho " << neighbourside;
                         out << "\nshape local do vizinho " << (neighlocshape+k) << endl;
			 //AvisoAudioVisual(5);
                     }
                  }
               }
            }
            out << gel->Id() << "/" << side << " -> ";
            out << neighbour->Id() << "/" << neighbourside << endl;
            neighcycle = neighcycle.Neighbour();
      	}//while
         delete sideintrule;
      }
   }
}



void ComputationalTestSide(TPZCompMesh *cmesh, ofstream &out) {

   int nelems = cmesh->ElementVec().NElements();

   int nsides;

   TPZInterpolatedElement *cel;

   TPZGeoEl *gel;

   for(int el=0;el<nelems;el++) {

      cel = (TPZInterpolatedElement *) cmesh->ElementVec()[el];

      if(!cel) continue;

      gel = cel->Reference();

      if(!gel) continue;

		nsides = gel->NSides();

      for(int side=0;side<nsides;side++) {

         TPZGeoEl *neighbour =  gel->Neighbour(side).Element();

         if(!neighbour) continue;

         TPZIntPoints *sideintrule = gel->CreateSideIntegrationRule(side,cel->SideOrder(side)*2);

         if(!sideintrule) {

            cout << "\nRegra de integracao nula\n";

            continue;

         }

         TPZGeoElSide neighcycle(neighbour,gel->Neighbour(side).Side());

         while(gel!=neighcycle.Element()) {

            neighbour = neighcycle.Element();

            TPZInterpolatedElement *neighref = (TPZInterpolatedElement *) neighbour->Reference();

            if(!neighref) break;

            int nsideshapes = cel->NSideShapeF(side);

            int neighbourside = neighcycle.Side();

            int nneighsideshapes = neighref->NSideShapeF(neighbourside);

            if(nsideshapes != nneighsideshapes) {

               cout << "\nnNumero de shapes incompatibles\n";

                out << "\nnNumero de shapes incompatibles\n";

               continue;

            }

            int dimside = gel->SideDimension(side);

            int neighdimside = neighbour->SideDimension(neighbourside);

            if(dimside!=neighdimside) {
               cout << "\nLados de dimensoes distintas\n";
                out << "\nLados de dimensoes distintas\n";
            }
            int nsideconnect = cel->NSideConnects(side);
            int nneighsideconnet = neighref->NSideConnects(neighbourside);
            if(nsideconnect!=nneighsideconnet) {
               cout << "\nNumero de connects distintos\n";
                out << "\nNumero de connects distintos\n";
            }

            //aqui nsideshape = nneighsideshape e dimside = neighdimside

            if(!neighdimside) neighdimside = 1;
            if(!dimside) dimside = 1;
            REAL phistore[50],dphistore[100],phineighstore[50],dphineighstore[100];
	    //suficiente para o side do quadrilatero de ordem 6
            TPZFMatrix phi(nsideshapes,1,phistore,50),phineigh(nneighsideshapes,1,phineighstore,50);
	    //phi(nsideshapes,1,0.),phineigh(nneighsideshapes,1,0.)
            TPZFMatrix dphi(dimside,nsideshapes,dphistore,100),dphineigh(neighdimside,nneighsideshapes,dphineighstore,100);
	    //dphi(dimside,nsideshapes,0.),dphineigh(neighdimside,nneighsideshapes,0.)

            //calcula a transformação entre lados vizinhos

            TPZGeoElSide geoside(gel,side),neighside(neighbour,neighbourside);

            TPZTransform t(dimside);

            t = geoside.NeighbourSideTransform(neighside);//geoside.SideTransform2(neighside,t);

            int numintp = sideintrule->NPoints();

            REAL weight,error;

            int nsideshp,locid,neighlocshape,locshape;

            for(int i=0;i<nsideconnect;i++) {

               locid = cel->SideConnectLocId(i,side);

               nsideshp = cel->NConnectShapeF(locid);

               locshape = 0;

               for(int l=0;l<i;l++) locshape += cel->NConnectShapeF(cel->SideConnectLocId(l,side));

               ShapeLocId(cel,neighref,locid,neighbourside,neighlocshape);

               for(int ip=0;ip<numintp;ip++) {

                  TPZVec<REAL> point(3,0.),pint(3,0.),pont(3,0.);

                  sideintrule->Point(ip,point,weight);

                  pint = point;

                  cel->SideShapeFunction(side,pint,phi,dphi);

                  t.Apply(point,pont);

                  neighref->SideShapeFunction(neighbourside,pont,phineigh,dphineigh);

                  for(int k=0;k<nsideshp;k++) {

                     error = fabs(phi(locshape+k,0)-phineigh(neighlocshape+k,0));

                     if(error > 1.e-10) {
                        cout << "\nerro no valor da funcao ->\n";
                        cout << "\nid do elemento " << gel->Id();
                        cout << "\nside do elemento " << side;
                        cout << "\nshape local dos elementos " << (locshape+k) << endl;
                        cout << "\nid do vizinho " << neighbour->Id();
                        cout << "\nside do vizinho " << neighbourside;
                        cout << "\nshape local dos elementos " << (neighlocshape+k) << endl;
                         out << "\nerro no valor da funcao ";
                         out << "\nid do elemento " << gel->Id();
                         out << "\nside do elemento " << side;
                         out << "\nshape local dos elementos " << (locshape+k) << endl;
                         out << "\nid do vizinho " << neighbour->Id();
                         out << "\nside do vizinho " << neighbourside;
                         out << "\nshape local do elemento " << (neighlocshape+k) << endl;
			 //AvisoAudioVisual(5);
                     }
                  }
               }
            }
            out << gel->Id() << "/" << side << " -> ";
            out << neighbour->Id() << "/" << neighbourside << endl;
            neighcycle = neighcycle.Neighbour();
      	}//while
         delete sideintrule;
      }
   }
}


void ComputationalTestInnerSide(TPZCompMesh *cmesh, ofstream &out) {

   int nelems = cmesh->ElementVec().NElements();

   int nsides;

   TPZInterpolatedElement *cel;

   TPZGeoEl *gel;

   for(int el=0;el<nelems;el++) {

      cel = (TPZInterpolatedElement *) cmesh->ElementVec()[el];

      if(!cel) continue;

      gel = cel->Reference();

      if(!gel) continue;

		nsides = gel->NSides();

      for(int side=0;side<nsides;side++) {

         TPZGeoEl *neighbour =  gel->Neighbour(side).Element();

         if(!neighbour) continue;

         TPZIntPoints *sideintrule = gel->CreateSideIntegrationRule(side,cel->SideOrder(side)*2);

         if(!sideintrule) {

            cout << "\nRegra de integracao nula\n";

            continue;

         }

         TPZGeoElSide neighcycle(neighbour,gel->Neighbour(side).Side());

         while(gel!=neighcycle.Element()) {

            neighbour = neighcycle.Element();

            TPZInterpolatedElement *neighref = (TPZInterpolatedElement *) neighbour->Reference();

            if(!neighref) break;

            int nsideshapes = cel->NSideShapeF(side);

            int neighbourside = neighcycle.Side();

            int nneighsideshapes = neighref->NSideShapeF(neighbourside);

            if(nsideshapes != nneighsideshapes) {

               cout << "\nnNumero de shapes incompatibles\n";

                out << "\nnNumero de shapes incompatibles\n";

               continue;

            }

            int dimside = gel->SideDimension(side);

            int neighdimside = neighbour->SideDimension(neighbourside);

            if(dimside!=neighdimside) {

               cout << "\nLados de dimensoes distintas\n";

                out << "\nLados de dimensoes distintas\n";

            }

            int nsideconnect = cel->NSideConnects(side);

            int nneighsideconnet = neighref->NSideConnects(neighbourside);

            if(nsideconnect!=nneighsideconnet) {

               cout << "\nNumero de connects distintos\n";

                out << "\nNumero de connects distintos\n";

            }

            //aqui nsideshape = nneighsideshape e dimside = neighdimside

            int nshape = cel->NShapeF();

            int dim = cel->Dimension();

            if(!neighdimside) neighdimside = 1;

            REAL phistore[220],dphistore[660],phineighstore[50],dphineighstore[100];
	    //suficientes para ordem 5 do cubo (o mais grande)
            TPZFMatrix phi(nshape,1,phistore,220),phineigh(nneighsideshapes,1,phineighstore,50);
	    //phi(nshape,1,0.),phineigh(nneighsideshapes,1,0.)
            TPZFMatrix dphi(dim,nshape,dphistore,660),dphineigh(neighdimside,nneighsideshapes,dphineighstore,100);
	    //dphi(dim,nshape,0.),dphineigh(neighdimside,nneighsideshapes,0.)

            //calcula a transformação entre lados vizinhos

            TPZGeoElSide geoside(gel,side),neighside(neighbour,neighbourside);

            TPZTransform t(dimside);

            t = geoside.NeighbourSideTransform(neighside);//geoside.SideTransform2(neighside,t);

            int numintp = sideintrule->NPoints();

            REAL weight,error;

            int nsideshp,locid,neighlocshape,locshape,cindex;

            for(int i=0;i<nsideconnect;i++) {

               locid = cel->SideConnectLocId(i,side);//id local do connect i do side

               cindex = cel->ConnectIndex(locid);//index do connect

               TPZConnect *connect = &(cmesh->ConnectVec()[cindex]);

               locshape = FirstSideShape(cmesh,cel,connect);

               nsideshp = cel->NConnectShapeF(locid);

               ShapeLocId(cel,neighref,locid,neighbourside,neighlocshape);

               for(int ip=0;ip<numintp;ip++) {

                  TPZVec<REAL> point(3,0.),pint(3,0.),pont(3,0.);

                  sideintrule->Point(ip,point,weight);

                  pint = point;

                  SideParameterToElement(gel->NSides(),side,pint);

                  cel->Shape(pint,phi,dphi);

                  t.Apply(point,pont);

                  neighref->SideShapeFunction(neighbourside,pont,phineigh,dphineigh);

                  for(int k=0;k<nsideshp;k++) {

                     error = fabs(phi(locshape+k,0)-phineigh(neighlocshape+k,0));

                     if(error > 1.e-10) {
                        cout << "\nerro no valor da funcao ";
                        cout << "\nid do elemento " << gel->Id();
                        cout << "\nside do elemento " << side;
                        cout << "\nshape local do elemento " << (locshape+k);
                        cout << "\nid do vizinho " << neighbour->Id();
                        cout << "\nside do vizinho " << neighbourside;
                        cout << "\nshape local do vizinho " << (neighlocshape+k) << endl;
                         out << "\nerro no valor da funcao ";
                         out << "\nid do elemento " << gel->Id();
                         out << "\nside do elemento " << side;
                         out << "\nshape local do elemento " << (locshape+k);
                         out << "\nid do vizinho " << neighbour->Id();
                         out << "\nside do vizinho " << neighbourside;
                         out << "\nshape local do vizinho " << (neighlocshape+k) << endl;
			 //AvisoAudioVisual(5);
                     }
                  }
               }
            }
            out << gel->Id() << "/" << side << " -> ";
            out << neighbour->Id() << "/" << neighbourside << endl;
            neighcycle = neighcycle.Neighbour();
      	}//while
         delete sideintrule;
      }
   }
}


void ShapeLocId(TPZInterpolatedElement *cel,TPZInterpolatedElement *neighref,int locid,int neighside,int &neighlocshape) {



	int index = cel->ConnectIndex(locid);

   TPZCompMesh *cmesh = cel->Mesh();

   TPZConnect *connect = &(cmesh->ConnectVec()[index]);

   //TPZConnect connect = &(cmesh->ConnectVec()[index]);//esta forma serve

   //TPZConnect &connect = cmesh->ConnectVec()[index];//esta forma não serve-----|

   //int seq = connect.SequenceNumber(),neighseq,;                            // |

   int nsideconn = neighref->NSideConnects(neighside),i;                      // |

   neighlocshape = 0;                                                         // |

   for(i=0;i<nsideconn;i++) {                                                 // |

	   locid = neighref->SideConnectLocId(i,neighside);                        // |

      index = neighref->ConnectIndex(locid);                                  // |

      //TPZConnect connect = cmesh->ConnectVec()[index];//da erro ao atualizar <-

      TPZConnect *actconnect = &(cmesh->ConnectVec()[index]);

      //neighseq = connect.SequenceNumber();

      //if(neighseq == seq) break;

      if(connect == actconnect) break;

      neighlocshape += neighref->NConnectShapeF(locid);

   }

}



int FirstSideShape(TPZCompMesh *cmesh,TPZInterpolatedElement *cel,TPZConnect *connect) {



   int seq = connect->SequenceNumber();

   int ncon = cel->NConnects();

   int locshape = 0,index;

   TPZConnect *elconn;

   for(int i=0;i<ncon;i++) {

   	index = cel->ConnectIndex(i);

      elconn = &(cmesh->ConnectVec()[index]);

      if(seq == elconn->SequenceNumber()) break;

      locshape += cel->NConnectShapeF(i);

   }

   return locshape;

}



void LerMalha(TPZGeoMesh *geomesh,TPZCompMesh *&compmesh) {

   TPZFMatrix xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);

   TPZMaterialTest3D *mat1d,*mat2d,*mat3d;//TPZMat1dLin *mat1d;//TPZMat2dLin *mat2d;

   int nnode,nel,ncorners,typeel;

   TPZVec<REAL> coord(3);

   cout << "Entre tipo de elemento\n"
        << "0: triangulo\n"
        << "1: quadrilatero\n"
        << "2: tetraedro\n"
        << "3: piramide\n"
        << "4: prisma\n"
        << "5: hexaedro\n";

   cin >> typeel;

   ConstroeNodes(typeel,geomesh);

   if(typeel==0) {ncorners = 3; nel = 05;}//até faces superpostas
   if(typeel==1) {ncorners = 4; nel = 06;}
   if(typeel==2) {ncorners = 4; nel = 11;}
   if(typeel==3) {ncorners = 5; nel = 14;}
   if(typeel==4) {ncorners = 6; nel = 15;}
   if(typeel==5) {ncorners = 8; nel = 19;}
   
   compmesh = new TPZCompMesh(geomesh);
   int numel = 0;

   for(int el=0;el<nel;el++)  {

      int mat=1,ntype;

      TPZVec<int> nodes(0);//faz resize no próximo passo

      IdsDefinitions(typeel,nodes,ntype,numel);//cosntroe elemento, arestas e faces

      switch(ntype) {//tipo de elemento

         case 10://unidimensional ; elg1d =

            new TPZGeoEl1d(nodes,mat,*geomesh);

            mat1d = new TPZMaterialTest3D(mat);//mat1d = new TPZMat1dLin(mat);

            compmesh->InsertMaterialObject(mat1d);

            mat1d->SetMaterial(xk);//mat1d->SetMaterial(xk,xb,xc,xf);

            break;

         case 0://triângulo ; elgt2d =

         	new TPZGeoElT2d(nodes,mat,*geomesh);

            mat2d = new TPZMaterialTest3D(mat);//mat2d = new TPZMat2dLin(mat);

            compmesh->InsertMaterialObject(mat2d);//mat2d->SetMaterial(xk,xc,xf);

            mat2d->SetMaterial(xk);

            break;

         case 1://quadrilátero ; elgq2d =

            new TPZGeoElQ2d(nodes,mat,*geomesh);

            mat2d = new TPZMaterialTest3D(mat);//mat2d = new TPZMat2dLin(mat);

            compmesh->InsertMaterialObject(mat2d);//mat2d->SetMaterial(xk,xc,xf);

            mat2d->SetMaterial(xk);

            break;

         case 2://tetraedro ; elgt3d =

            new TPZGeoElT3d(nodes,mat,*geomesh);

            mat3d = new TPZMaterialTest3D(mat);

            compmesh->InsertMaterialObject(mat3d);

            mat3d->SetMaterial(xk);

            break;

         case 3://pirâmide ; elgpi3d =

            new TPZGeoElPi3d(nodes,mat,*geomesh);

            mat3d = new TPZMaterialTest3D(mat);

	    compmesh->InsertMaterialObject(mat3d);

            mat3d->SetMaterial(xk);

            break;

         case 4://prismae ; elgpi3d =

            new TPZGeoElPr3d(nodes,mat,*geomesh);

            mat3d = new TPZMaterialTest3D(mat);

	    compmesh->InsertMaterialObject(mat3d);

            mat3d->SetMaterial(xk);

            break;

         case 5://cubo ; elgc3d =

            new TPZGeoElC3d(nodes,mat,*geomesh);

            mat3d = new TPZMaterialTest3D(mat);

   	    compmesh->InsertMaterialObject(mat3d);

            mat3d->SetMaterial(xk);

	    break;

         default:

         	for(int i=0;i<300;i++)

            	cout << "\nmain::LerMalha -> Elemento nao conhecido\n";

            cout << "\nChao\n";

            exit(1);

            delete mat1d;

            delete mat2d;

            delete mat3d;

      }

   }

}



void SideParameterToElement1d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

  if(side == 0) point[0] = -1.;

  else if(side == 1) point[0] = 1;

  else {

    point[0] = par[0];

  }

}



void SideParameterToElementT2d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

  REAL aux;

  //point[0] = point[1] = 0.;//não pode

  switch(side) {

  case 0:

    return;

  case 1:

    point[0] = 1.;

    point[1] = 0.;

    return;

  case 2:

  	 point[0] = 0.;

    point[1] = 1.;

    return;

  case 3:

    point[0] = (1.+par[0])*.5;

    point[1] = 0.;

	 return;

  case 4:

	 aux = par[0];

    point[0] = (1.-aux)*.5;

	 point[1] = (1.+aux)*.5;

    return;

  case 5:

    point[1] = (1.-par[0])*.5;

    point[0] = 0.;

    return;

  case 6:

    point[0] = par[0];

    point[1] = par[1];

    return;

  default:

    PZError << "TPZCompElT2d::SideParameterToElement. Bad paramenter side.\n";

  }

}



void SideParameterToElementQ2d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

  switch(side) {//0 a 8 quadrilatero

     case 0:

       point[0]=point[1]=-1.;

       return;

     case 1:

       point[0]= 1.;

       point[1]=-1.;

       return;

     case 2:

       point[0]=point[1]=1.;

       return;

     case 3:

       point[0]=-1.;

       point[1]=1.;

       return;

     case 4:

       point[0] = par[0];

       point[1] = -1.0;

       return;

     case 5:

       point[1] = par[0];

       point[0] = 1.;

       return;

     case 6:

       point[0] = -par[0];

       point[1] = 1.0;

       return;

     case 7:

       point[1] = -par[0];

       point[0] = -1.;

       return;

     case 8:

       point[0]=par[0];//ponto

       point[1]=par[1];//de integração

       return;

     default:

       PZError << "TPZCompElQ2d::SideParameterToElement. Bad paramenter side.\n";

       return;

  }

}



void SideParameterToElementT3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

   //transforma ponto da regra de integração do lado para o elemento todo

   if(side>-1 && side<4) {//cantos

     point[0] = TPZCompElT3d::MasterCoord[side][0];//8x3

     point[1] = TPZCompElT3d::MasterCoord[side][1];

     point[2] = TPZCompElT3d::MasterCoord[side][2];

     return;

   }

   REAL aux;

   //point[0] = point[1] = point[2] = 0.;//não pode

   switch(side) {

  	  //lados

     case  4:

        point[0] = (1.+par[0])*.5;

        point[1] = 0.;

        point[2] = 0.;

        return;

     case  5:

        aux = par[0];

        point[0] = (1.-aux)*.5;

        point[1] = (1.+aux)*.5;

        point[2] = 0.;

        return;

     case 6:

        point[1] = (1.-par[0])*.5;

        point[0] = 0.;

        point[2] = 0.;

        return;

     case 7:

        point[1] = 0.;

        point[2] = (1.+par[0])*.5;

        point[0] = 0.;

        return;

     case 8:

        aux = par[0];

        point[0] = (1.-aux)*.5;

        point[1] = 0.;

        point[2] = (1.+aux)*.5;

        return;

     case 9:

        aux = par[0];

        point[0] = 0.;

        point[1] = (1.-aux)*.5;

        point[2] = (1.+aux)*.5;

        return;

     //faces

     case 10:

        point[0] = par[0];

        point[1] = par[1];

        point[2] = 0.;

        return;

     case 11:

        point[0] = par[0];

        point[2] = par[1];

        point[1] = 0.;

        return;

     case 12:

     	  aux = par[0];

        point[2] = par[1];

        point[0] = 1.-aux-par[1];

        point[1] = aux;

        return;

     case 13:

     	  aux = par[0];

        point[0] = 0.;

        point[2] = par[1];

        point[1] = aux;

        return;

     //interior

     case 14:

        point[0] = par[0];//ponto

        point[1] = par[1];//de integração

        point[2] = par[2];

       return;

     default:

       PZError << "TPZCompElT3d::SideParameterToElement. Bad paramenter side.\n";

   }

}



void SideParameterToElementPi3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

   //transforma ponto da regra de integração do lado para o elemento todo

   //basta fazer isto só para o elemento mestre, a função X(par1,par2) corrige para o elemento deformado

   if(side>-1 && side<5) {//cantos

     point[0] = TPZCompElPi3d::MasterCoord[side][0];//8x3

     point[1] = TPZCompElPi3d::MasterCoord[side][1];

     point[2] = TPZCompElPi3d::MasterCoord[side][2];

     return;

   }

   REAL aux;

   switch(side) {

  	  //lados

     case  5:

        point[0] = par[0];

        point[1] = -1.;

        point[2] =  0.;

        return;

     case  6:

     	  point[1] = par[0];

        point[0] = 1.;

        point[2] = 0.;

        return;

     case 7:

        point[0] = -par[0];

        point[1] =  1.;

        point[2] =  0.;

        return;

     case 8:

        point[1] = -par[0];

        point[0] = -1.;

        point[2] =  0.;

        return;

     case 9:

        point[1] = .5*(par[0]-1.);

        point[2] = .5*(par[0]+1.);

        point[0] = .5*(par[0]-1.);

        return;

     case 10:

        point[1] = .5*(par[0]-1.);

        point[2] = .5*(par[0]+1.);

        point[0] =-.5*(par[0]-1.);

        return;

     case 11:

        point[1] =-.5*(par[0]-1.);

        point[2] = .5*(par[0]+1.);

        point[0] =-.5*(par[0]-1.);

        return;

     case 12:

        point[1] =-.5*(par[0]-1.);

        point[2] = .5*(par[0]+1.);

        point[0] = .5*(par[0]-1.);

        return;

     //faces

     case 13:

        point[0] = par[0];//não precisa declarar

        point[1] = par[1];//não precisa declarar

        point[2] = 0.;

        return;

     case 14:

        point[0] = 2.*par[0]+par[1]-1.;

        point[2] = par[1];

        point[1] = par[1]-1.;

        return;

     case 15:

        aux = par[0];

        point[0] = -par[1]+1.;

        point[2] = par[1];

        point[1] = 2.*aux+par[1]-1.;

        return;

     case 16:

        aux = par[1];

        point[1] = -par[1]+1.;

        point[2] = aux;

        point[0] = 2.*par[0]+aux-1.;

        return;

     case 17:

        aux = par[0];

        point[0] = par[1]-1.;

        point[2] = par[1];

        point[1] = 2.*aux+par[1]-1.;

        return;

     //interior

     case 18:

        point[0] = par[0];//o proprio ponto

        point[1] = par[1];//de integração

        point[2] = par[2];//

       return;

     default:

       PZError << "TPZCompElPi3d::SideParameterToElement. Bad paramenter side.\n";

   }

}



void SideParameterToElementPr3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

   //transforma ponto da regra de integração do lado para o elemento todo

   //basta fazer isto só para o elemento mestre, a função X(par1,par2) corrige para o elemento deformado

   if(side>-1 && side<6) {//cantos

     point[0] = TPZCompElPr3d::MasterCoord[side][0];//6x3

     point[1] = TPZCompElPr3d::MasterCoord[side][1];

     point[2] = TPZCompElPr3d::MasterCoord[side][2];

     return;

   }

   REAL aux;

   switch(side) {

     case  6:

        point[0] = (1.+par[0])*.5;

        point[1] =  0.;

        point[2] = -1.;

        return;

     case  7:

        aux = par[0];

        point[0] = (1.-aux)*.5;

        point[1] = (1.+aux)*.5;

        point[2] = -1.;

        return;

     case 8:

        point[1] = (1.-par[0])*.5;//(0,pt,-1)

        point[0] =  0.;//do canto 2 para o canto 0

        point[2] = -1.;

        return;

     case 9:

        point[2] = par[0];

        point[0] = 0.;

        point[1] = 0.;

        return;

     case 10:

        point[2] = par[0];

        point[0] = 1.;

        point[1] = 0.;

        return;

     case 11:

        point[2] = par[0];

        point[0] = 0.;

        point[1] = 1.;

        return;

     case 12:

        point[0] = .5*(1.+par[0]);

        point[1] = 0.;

        point[2] = 1.;

        return;

     case 13:

        point[1] = .5*(1.+par[0]);// do canto 4 para o canto 5

        point[0] = .5*(1.-par[0]);

        point[2] = 1.;

        return;

     case 14:

        point[1] = .5*(1.-par[0]);//-1<=pt[0]<=1 => 1>=point[1]>=0

        point[0] = 0.;           //no sentido da aresta 14

        point[2] = 1.;

        return;

     //faces

     case 15:

        point[0] = par[0];

        point[1] = par[1];

        point[2] = -1.;

        return;

     case 16:

        point[0] = .5*(par[0]+1.);

        point[2] = par[1];

        point[1] = 0.;

        return;

     case 17:

     	  aux = par[0];

        point[0] = .5*(1.-aux);

        point[2] = par[1];

        point[1] = .5*(1.+aux);

        return;

     case 18:

        point[2] = par[1];//regra para quadrilatero [-1,1]x[-1,1]

        point[1] = .5*(par[0]+1.);//-1<=par[0]<=1 => 0<=point[1]<=1

        point[0] = 0.;

       return;

     case 19:

        point[0] = par[0];

        point[1] = par[1];

        point[2] = 1.;

        return;

     //interior

     case 20:

        point[0] = par[0];//o proprio ponto

        point[1] = par[1];//de integração

        point[2] = par[2];//

       return;

     default:

       PZError << "TPZCompElPi3d::SideParameterToElement. Bad paramenter side.\n";

   }

}

void SideParameterToElementC3d(int side,TPZVec<REAL> &par,TPZVec<REAL> &point) {

   //transforma ponto da regra de integração do lado para o elemento todo

   //basta fazer isto só para o elemento mestre, a função X(par1,par2) corrige para o elemento deformado

   if(side>-1 && side<8) {//cantos

     point[0] = TPZCompElC3d::MasterCoord[side][0];//8x3

     point[1] = TPZCompElC3d::MasterCoord[side][1];

     point[2] = TPZCompElC3d::MasterCoord[side][2];

     return;

   }

   REAL aux;

   switch(side) {

  	  //lados

     case  8:

        point[0] = par[0];

        point[1] = -1.;

        point[2] = -1.;

        return;

     case  9:

     	  point[1] = par[0];

        point[0] = 1.;

        point[2] = -1.;

        return;

     case 10:

        point[0] = -par[0];

        point[1] =  1.;

        point[2] = -1.;

        return;

     case 11:

        point[1] = -par[0];

        point[0] = -1.;

        point[2] = -1.;

        return;

     case 12:

        point[2] = par[0];

        point[0] = -1.;

        point[1] = -1.;

        return;

     case 13:

        point[2] = par[0];

        point[0] = 1.;

        point[1] = -1.;

        return;

     case 14:

        point[2] = par[0];

        point[0] = 1.;

        point[1] = 1.;

        return;

     case 15:

        point[2] = par[0];

        point[0] = -1.;

        point[1] = 1.;

        return;

     case 16:

        point[0] = par[0];

        point[1] = -1.;

        point[2] =  1.;

        return;

     case 17:

        point[1] = par[0];

        point[0] = 1.;

        point[2] = 1.;

        return;

     case 18:

        point[0] = -par[0];

        point[1] = 1.;

        point[2] = 1.;

        return;

     case 19:

        point[1] = -par[0];

        point[0] = -1.;

        point[2] = 1.;

        return;

     //faces

     case 20:             //seguindo o sentido exterior a face deve ser

        point[0] = par[0];//-par[0];

        point[1] = par[1];//

        point[2] = -1.;

        return;

     case 21:

        point[0] = par[0];

        point[2] = par[1];

        point[1] = -1.;

        return;

     case 22:

        aux = par[1];

        point[1] = par[0];

        point[0] = 1.;

        point[2] = aux;

        return;

     case 23:

        point[0] = par[0];//-par[0];

        point[2] = par[1];

        point[1] = 1.;

        return;

     case 24:

        aux = par[1];

        point[1] = par[0];

        point[0] = -1.;

        point[2] = aux;//-par[0];

        return;

     case 25:

        point[0] = par[0];

        point[1] = par[1];

        point[2] = 1.;

        return;

     //interior

     case 26:

        point[0] = par[0];//ponto

        point[1] = par[1];//de integração

        point[2] = par[2];

       return;

     default:

       PZError << "TPZCompElC3d::SideParameterToElement. Bad paramenter side.\n";

   }

}





void CC(TPZGeoMesh *firstmesh) {

   TPZBndCond *bc;

   REAL big  = 1.e12;

   TPZFMatrix val1(1,1,big), val2(1,1,0.);

   int side,ncc,type;

   ifstream in("CCon.out");

   in >> ncc;

   int iel=0,idc=0;

   while(ncc) {

      TPZGeoEl *elgi = firstmesh->ElementVec()[iel++];

      if(elgi) {

         for(int k=0;k<ncc;k++) {

            in >> side >> type;

            TPZGeoElBC(elgi,side,--idc,*firstmesh);

            bc = mat->CreateBC(idc,type,val1,val2);

            firstmesh->Reference()->InsertMaterialObject(bc);

            if(type==0) bc->SetForcingFunction(BCDirichlet);

         }

      }

      in >> ncc;

   }

}



void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f) {

   REAL x = point[0];

   REAL y = point[1];

   REAL z = point[2];

   //f[0]  = pow(x,p)*pow(y,q)*pow(z,r);

   f[0]  = pow(x,1)*pow(y,1)*pow(z,1);

}



void Divide(TPZCompMesh *compmesh) {



   TPZVec<int> csub;

   int n1=1;

   while(n1) {

	   cout << "Id do elemento geometrico a dividir ? : ";

      cin >> n1;

      if(n1 < 0) break;

      int nelc = compmesh->ElementVec().NElements();

      int el;

      TPZCompEl *cpel;

      for(el=0;el<nelc;el++) {

         cpel = compmesh->ElementVec()[el];

         if(cpel && cpel->Reference()->Id() == n1) break;

      }

      cout << "\nNao interpola ! : Divide(el,csub,0)\n";

      compmesh->Divide(el,csub,0);

      n1 = 1;

   }

}



void NeighbourhoodCycles(TPZGeoMesh &gmesh,ofstream &out) {

  TPZAdmChunkVector<TPZGeoEl *> gm = gmesh.ElementVec();
  int i,j,k,nsides,count;
  int nel = gm.NElements();
  TPZGeoElSide neigh,neighkeep;
  for(i=0;i<nel;i++) {
    TPZGeoEl *gel = gm[i];
    if(!gel) continue;
    nsides = gel->NSides();    
    out << "Elemento geometrico de id = " << gel->Id() << endl;    
    for(j=0;j<nsides;j++) {      
      count = 1;      
      TPZGeoElSide gelside(gel,j);      
      if(gelside.Exists()) neigh = gelside.Neighbour();      
      if(!neigh.Exists()) {	
	out << "Lado sem vizinhanza : " << j << endl;	
	continue; 	
      }
      int nsel = gel->NSideNodes(j);
      int nsneigh = neigh.Element()->NSideNodes(neigh.Side());
      if(nsel!=nsneigh) {
	cout << "\nSides incompatíveis\n";
	continue;
      }


      TPZVec<REAL> nodel(3,0.),nodneigh(3,0.);
      gel->CenterPoint(j,nodel);
      int sidegel = gel->WhichSide(nodel);
      neigh.Element()->CenterPoint(neigh.Side(),nodneigh);
      int sideneigh = neigh.Element()->WhichSide(nodneigh);


      if(sidegel!=gelside.Side() || sideneigh!=neigh.Side()) {
	cout << "\nSides incompativeis\n";
	continue;
      }
      while(neigh.Element() && neigh!=gelside) {
	count++;//>1
	neigh = neigh.Neighbour();
	if(count > 1001) {
	  cout << "Erro 1 : Ciclo degenerado ou com mais de mil elementos, continua o ciclo (0/1)? : ";
	  cin >> count;
	}
	if(count == 0) {
	  out << "Erro 1 : Ciclo degenerado ou com mais de mil elementos\n";
	  out << gelside.Element()->Id() << "/" << gelside.Side() << " ";
	  neigh = gelside.Neighbour();
	  for(k=0;k<30;k++) {
	    out << neigh.Element()->Id() << "/" << neigh.Side() << " ";
	    neigh = neigh.Neighbour();
	  }
	  out << endl;
	}
      }//while
      if(!neigh.Element()) {
	out << "\nErro 2 : Ciclo interrupto\n";
	cout << "Erro 2 : Ciclo interrupto, continua (1)! : ";
	cin >> k;
      }
    }//for lados
  }//for elementos
}





void SideParameterToElement(int nsides,int side,TPZVec<REAL> &pint) {

   switch(nsides) {//ncorners

      case 3://linear : 2

         SideParameterToElement1d(side,pint,pint); break;

      case 7://triângulo : 3

         SideParameterToElementT2d(side,pint,pint); break;

      case 9://quadrilátero : 4

         SideParameterToElementQ2d(side,pint,pint); break;

      case 15://tetraedro : 4

         SideParameterToElementT3d(side,pint,pint); break;

      case 19://pirâmide : 5

         SideParameterToElementPi3d(side,pint,pint); break;

      case 21://prisma : 6

         SideParameterToElementPr3d(side,pint,pint); break;

      case 27://cubo : 8

         SideParameterToElementC3d(side,pint,pint);

         break;

      default:

         cout << "\nmain::SideParameterToElement elemento nao conhecido\n";

   }

}



ofstream divide("divide.out");

void AutomaticDivide(TPZCompMesh &cmesh) {



   int maxlevel;

   cout << "\nNivel maximo a ser dividido : ";

   cin >> maxlevel;

   TPZVec<int> csub;

   int el,elid,numeldiv=0;

   TPZCompEl *cpel;

   TPZGeoEl *gel;

   TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh.ElementVec();

   int level = 0,nelc;

   cout << "\nInterpolar a divisao (0/1)? ";

   int nao;

   cin >> nao;

   while(level < (maxlevel+1)) {

      el = 0;

      nelc = elementvec.NElements();

      numeldiv = 0;

      while(el < nelc) {

         cpel = elementvec[el];

         if(cpel) gel = cpel->Reference();

         else gel = 0;

         if(cpel && gel && gel->Level() < (level+1)) {

            elid = gel->Id();

         	cmesh.Divide(el,csub,nao);

            numeldiv++;

         }

         el++;

      }

      divide << "nivel dividido " << level << endl;

      divide << "id do ultimo elemento dividido deste nivel : " << elid << endl;

      level++;

   }

   divide << "\nNumero de elementos divididos no ultimo nivel : " << numeldiv;

   //verificando a divisao

   nelc = elementvec.NElements();

   el = 0;

   while(el < nelc) {

      cpel = elementvec[el];

      if(cpel) gel = cpel->Reference();

      else gel = 0;

      if(cpel && gel) {

         if(gel->Level()!=level)

	         divide << "\n\nexistem elementos de niveis diferentes : \n";

      }

      el++;

   }

}


void ConstroeNodes(int typeel,TPZGeoMesh *geomesh){

  TPZVec<REAL> coord(3,0.);
  if(typeel == 0){//triangulo
     geomesh->NodeVec().Resize(3);
     for(int inode=0;inode<3;inode++) {
       coord[0] = triangulo[0][inode];
       coord[1] = triangulo[1][inode];
       coord[2] = triangulo[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
  if(typeel == 1){//quadrilatero 
     geomesh->NodeVec().Resize(4);
     for(int inode=0;inode<4;inode++) {
       coord[0] = quadrilatero[0][inode];
       coord[1] = quadrilatero[1][inode];
       coord[2] = quadrilatero[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
  if(typeel == 2){//tetraedro
     geomesh->NodeVec().Resize(4);
     for(int inode=0;inode<4;inode++) {
       coord[0] = tetraedro[0][inode];
       coord[1] = tetraedro[1][inode];
       coord[2] = tetraedro[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
  if(typeel == 3){//piramide
     geomesh->NodeVec().Resize(5);
     for(int inode=0;inode<5;inode++) {
       coord[0] = piramide[0][inode];
       coord[1] = piramide[1][inode];
       coord[2] = piramide[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
  if(typeel == 4){//prisma
     geomesh->NodeVec().Resize(6);
     for(int inode=0;inode<6;inode++) {
       coord[0] = prisma[0][inode];
       coord[1] = prisma[1][inode];
       coord[2] = prisma[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
  if(typeel == 5){//hexaedro
     geomesh->NodeVec().Resize(8);
     for(int inode=0;inode<8;inode++) {
       coord[0] = hexaedro[0][inode];
       coord[1] = hexaedro[1][inode];
       coord[2] = hexaedro[2][inode];
       geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
     }
     return;
   }
}

void IdsDefinitions(int typeel,TPZVec<int> &nodes,int &ntype,int &numel){

  int i;
  if(typeel == 0){//triangulo
    switch(numel){
    case 0:
    case 1:
      ntype = 0;//triangulo
      nodes.Resize(3);
      for(i=0;i<3;i++) nodes[i] = i;
      break;
    case 2:
      ntype = 10;//linha
      nodes.Resize(2);
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 2;
      nodes[1] = 0;
    }    
    numel++;
    return;
  }
  if(typeel == 1){//quadrilatero
    switch(numel){
    case 0:
    case 1:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      for(i=0;i<4;i++) nodes[i] = i;
      break;
    case 2:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 3;
      break;
    case 5:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 0;
    }
    numel++;
    return;
  }
  if(typeel == 2){//tetraedro
    switch(numel){
    case 0:
      nodes.Resize(4);
      ntype = 2;//tetraedro
      for(i=0;i<4;i++) nodes[i] = i;
      break;
    case 1:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 2:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 0;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 3;
      break;
    case 5:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 3;
      break;
    case 6:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 3;
      break;
    case 7:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      break;
    case 8:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 3;
      break;
    case 9:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 2;
      nodes[2] = 3;
      break;
    case 10:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 1;
      nodes[1] = 2;
      nodes[2] = 3;
    }
    numel++;
    return;
  }
  if(typeel == 3){//piramide
    switch(numel){
    case 0:
      nodes.Resize(5);
      ntype = 3;//piramide
      for(i=0;i<5;i++) nodes[i] = i;
      break;
    case 1:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 2:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 3;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 0;
      break;
    case 5:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 4;
      break;
    case 6:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 4;
      break;
    case 7:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 4;
      break;
    case 8:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 4;
      break;
    case 9:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      nodes[3] = 3;
      break;
    case 10:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 4;
      break;
    case 11:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 1;
      nodes[1] = 2;
      nodes[2] = 4;
      break;
    case 12:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 3;
      nodes[1] = 2;
      nodes[2] = 4;
      break;
    case 13:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 3;
      nodes[2] = 4;
    }
    numel++;
    return;
  }
  if(typeel == 4){//prisma
    switch(numel){
    case 0:
      nodes.Resize(6);
      ntype = 4;//prisma
      for(i=0;i<6;i++) nodes[i] = i;
      break;
    case 1:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 2:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 0;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 3;
      break;
    case 5:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 4;
      break;
    case 6:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 5;
      break;
    case 7:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 4;
      break;
    case 8:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 4;
      nodes[1] = 5;
      break;
    case 9:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 5;
      nodes[1] = 3;
      break;
    case 10:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      break;
    case 11:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 4;
      nodes[3] = 3;
      break;
    case 12:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 1;
      nodes[1] = 2;
      nodes[2] = 5;
      nodes[3] = 4;
      break;
    case 13:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 2;
      nodes[2] = 5;
      nodes[3] = 3;
      break;
    case 14:
      nodes.Resize(3);
      ntype = 0;//triangulo
      nodes[0] = 3;
      nodes[1] = 4;
      nodes[2] = 5;
    }
    numel++;
    return;
  }
  if(typeel == 5){//hexaedro
    switch(numel){
    case 0:
      nodes.Resize(8);
      ntype = 5;//prisma
      for(i=0;i<8;i++) nodes[i] = i;
      break;
    case 1:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 1;
      break;
    case 2:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 2;
      break;
    case 3:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 3;
      break;
    case 4:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 0;
      break;
    case 5:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 0;
      nodes[1] = 4;
      break;
    case 6:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 1;
      nodes[1] = 5;
      break;
    case 7:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 2;
      nodes[1] = 6;
      break;
    case 8:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 3;
      nodes[1] = 7;
      break;
    case 9:
      nodes.Resize(2);
      ntype = 10;//linha
      nodes[0] = 4;
      nodes[1] = 5;
      break;
    case 10:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 5;
      nodes[1] = 6;
      break;
    case 11:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 6;
      nodes[1] = 7;
      break;
    case 12:
      nodes.Resize(2);
      ntype = 10;
      nodes[0] = 7;
      nodes[1] = 4;
      break;
    case 13:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 2;
      nodes[3] = 3;
      break;
    case 14:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 1;
      nodes[2] = 5;
      nodes[3] = 4;
      break;
    case 15:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 1;
      nodes[1] = 2;
      nodes[2] = 6;
      nodes[3] = 5;
      break;
    case 16:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 3;
      nodes[1] = 2;
      nodes[2] = 6;
      nodes[3] = 7;
      break;
    case 17:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 0;
      nodes[1] = 3;
      nodes[2] = 7;
      nodes[3] = 4;
      break;
    case 18:
      nodes.Resize(4);
      ntype = 1;//quadrilatero
      nodes[0] = 4;
      nodes[1] = 5;
      nodes[2] = 6;
      nodes[3] = 7;
    }
    numel++;//contagem dos elementos criados
  }
}


void LargestNeighbour(TPZGeoElSide THIS,TPZGeoElSide &largeel, int targetdimension) {

  TPZGeoElSide fatherkeep(THIS.Element(),THIS.Side());

  if(THIS.Dimension() > targetdimension) {
      largeel = TPZGeoElSide();//fatherkeep;
      return;
  }

  TPZGeoElSide father = THIS.Father2();

  while(father.Element()) {//father.Dimension() = 1,2
    fatherkeep = father;
    father = father.Father2();
  }

  father = fatherkeep;
  if(THIS.Dimension() < targetdimension){
    TPZGeoEl *gel = THIS.Element();
    if(gel->NSides() == 6){
      father = ( THIS.Element())->HigherDimensionSides(THIS.Side(),targetdimension);
    }

  }
  if(!father.Element()) {
    largeel = THIS;//fatherkeep;
    return;
  }

  while(father.Element()) {//father.Dimension() = 1,2
    fatherkeep = father;
    father = father.Father2();
  }

  largeel = fatherkeep.Neighbour();
  if(!largeel.Element()) largeel = fatherkeep;
}

void TesteShapeFunctions(TPZGeoMesh *firstmesh){

  int sim = 1;
  TPZCompMesh *secondmesh = firstmesh->Reference();

  cout << "\nNeighbourhood test (0/1)! \n";
  //cin >> sim;
  if(sim) {
    ofstream neighbourhood("neighbourhood.out");
    NeighbourhoodCycles(*firstmesh,neighbourhood);
  }
  cout << "Geometric test (0/1)! \n";
  //cin >> sim;
  if(sim) {
    ofstream geometrictest("GeometricTest.out");
    GeometricTest(firstmesh,geometrictest);
  }
  cout << "Computational test inner-inner (0/1)! \n";
  //cin >> sim;
  if(sim) {
    ofstream computationaltest("ComputationalTest.out");
    ComputationalTest(secondmesh,computationaltest);
  }
  cout << "Computational test side-side (0/1)! \n";
  //cin >> sim;
  if(sim) {
    ofstream computationaltest("ComputationalSideTest.out");
    ComputationalTestSide(secondmesh,computationaltest);
  }
  cout << "Computational test inner-side (0/1)! \n";
  //cin >> sim;
  if(sim) {
    ofstream computationaltest("ComputationalInnerSideTest.out");
    ComputationalTestInnerSide(secondmesh,computationaltest);
  }
  cout << "\n     MAIN::TesteShapeFunctions teste = " << ++test << endl;
}

void Beep(int batida){

  clock_t start,end;
  start = clock();
  end = start;
  REAL tempo = 80.0;
  while( ((end - start)/tempo) < 10000) end = clock();
  int cont = 0;
  cout << "\n\t\t\tERRO DETECTADO :   \a" << batida;
  cout << "\n\n";
}

void AvisoAudioVisual(int nbat){

  int step = 0;
  while(step++ < nbat) Beep(step);
}
