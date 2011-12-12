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
#include "pzvec.h"
#include "pztrnsform.h"
#include "pzsolve.h"

#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzmattest3d.h"
#include <stdio.h>
#include <time.h>
#include <cstdlib>
#include <iostream>

using namespace std;

REAL Norm(REAL *vec);
void FirstFace(int &face0,int &nfaces,TPZGeoEl *gel);
REAL AngleFaceFi(TPZGeoEl *gel,int f,int node);
void ComputationalTest(TPZCompMesh *cmesh,ofstream &out);
void LerMalha(const char *malha,TPZGeoMesh *geomesh,TPZCompMesh *& compmesh);
void AutomaticDivide(TPZCompMesh &cmesh,int actuallevel);
void Divide(TPZCompMesh *compmesh);
void AngularTest(TPZGeoMesh *geomesh,ofstream &out);
void ElementName(TPZCompEl *cel,ofstream &out);
TPZMaterialTest3D *mat;

int main() {

   //malha geometrica
   TPZGeoMesh *gmesh = new TPZGeoMesh;
   //Arquivos de saida
	ofstream outgm("outgm.txt");
   ofstream outcm("outcm.txt");
   TPZCompMesh *compmesh;
   //malha geom�trica
   LerMalha("malha.txt",gmesh,compmesh);
   //ordem de interpolacao
   int ord;
   cout << "Enter order 1,2,3,4,5,... : \n";
   cin >> ord;
   TPZCompEl::SetgOrder(ord);
   //constru��o malha computacional
   compmesh->AutoBuild();
   //Divide computacional
   cout << "\nDivisao manual (0/1)? : ";
   int nao; cin >> nao;
   if (nao) Divide(compmesh);
   //AutomaticDivide(*compmesh,5);
   gmesh->Print(outgm);
   outgm.flush();
   compmesh->Print(outcm);
   outcm.flush();
   cout << "\n\nCalculo dos angulos nas faces do elemento\n\n";
   ofstream angle("anglemesh.txt");
   AngularTest(gmesh,angle);
   int sim;
   cout << "\nEnd # \n";
   cin >> sim;
   delete compmesh;
   delete gmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN
void AngularTest(TPZGeoMesh *geomesh,ofstream &out) {

   int el,maxlevel,level=0;
   REAL minangle=4.0,maxangle=0.0;
   cout << "\nMaximo nivel a ser atingido (nivel dividido+1) : ";
   cin >> maxlevel;
   TPZGeoEl *gel;
   TPZVec<REAL> angles(250000,4.0);//angulo t �: |t|<Pi
   while(level < maxlevel+1){
      REAL levelminangle=4.0,levelmaxangle=0.0;
      int nelg = geomesh->ElementVec().NElements();
      for(el=0;el<nelg;el++) {//percorre elementos geometricos
         gel = geomesh->ElementVec()[el];
         if(!gel) continue;
         TPZCompEl *gelcomp =(TPZCompEl *) gel->Reference();
         if(!gelcomp) {
            if(gel->Level() == level) cout << "\nmain::gelcomp nulo\n";
            continue;
         }
         int face0,nfaces;
         FirstFace(face0,nfaces,gel);//o side da primeira face
         out << endl;
         ElementName(gelcomp,out);
         out << "  Element id " << gel->Id() << "  level " << level << endl;
         for(int f=0;f<nfaces;f++) {//percorre faces
            int facei = f+face0;
            int nnodeface = gel->NSideNodes(facei);//num n�s da face
            for(int node=0;node<nnodeface;node++) {//percorre n�s das faces
               out << "face " << f << " : node local " << node;
               REAL angulo = AngleFaceFi(gel,f,node);
               if(angulo<levelminangle) levelminangle = angulo;
               if(angulo>levelmaxangle) levelmaxangle = angulo;
               if(angulo < 0.0) {cout << "\nFace degenerada!\n"; out << "\nFace degenerada!\n";}
               out << " : angulo " << angulo << endl;
               //guardando os angulos distintos
               int i=0;
               while(i<250000 && (4.0-angles[i]) > 0.5){//4-|t| > 0.5 onde |t|<Pi
                  if( fabs(angles[i] - angulo) < 1.e-12 ) break;
                  i++;
               }
               if(i>249999) cout << "\n\nmain::angles capacidade esgotada\n\n";
               if( fabs(4.0 - angles[i]) < 1.e-12 ) angles[i] = angulo;//angles.Push(angulo);
            }
         }
      }
      if(level < maxlevel) AutomaticDivide(*geomesh->Reference(),level);
      out << "\nMinimo angulo do nivel " << level << " : " << levelminangle;
      out << "\nMaximo angulo do nivel " << level << " : " << levelmaxangle << endl;
      if(levelminangle<minangle) minangle = levelminangle;
      if(levelmaxangle>maxangle) maxangle = levelmaxangle;
      level++;
   }
   out << "\nAngulos diferentes da malha\n";
   int i=0;
   while( (4.0-angles[i]) > 0.5 ) 
   {
     out << i << ": " << angles[i] << endl;
     i++;
   }
   out << "\nMinimo angulo da malha : " << minangle;
   out << "\nMaximo angulo da malha : " << maxangle;
}

REAL AngleFaceFi(TPZGeoEl *gel,int f,int node){

  TPZGeoElSide gels(gel,f);
  if(gels.Dimension() != 2) return 0.;
  int nsidenodes = gels.NSideNodes();
   REAL u[3],v[3];
   TPZVec<int> nodes(4);
   int i,index[3];
   index[0] = gels.SideNodeIndex(node%nsidenodes);
   index[1] = gels.SideNodeIndex((node+1)%nsidenodes);
   index[2] = gels.SideNodeIndex((node+nsidenodes-1)%nsidenodes);
   TPZGeoMesh *gm = gel->Mesh();
   TPZAdmChunkVector<TPZGeoNode> &nodevec = gm->NodeVec();
   for(i=0;i<3;i++){
    u[i] = nodevec[index[1]].Coord(i) - nodevec[index[0]].Coord(i);
    v[i] = nodevec[index[2]].Coord(i) - nodevec[index[0]].Coord(i);
   }
   REAL tetha = acos( (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]) / Norm(u) / Norm(v) );
   return tetha;
}

REAL Norm(REAL vec[3]){
	return (sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));
}

void FirstFace(int &face0,int &nfaces,TPZGeoEl *gel){

   int nsides = gel->NSides();
   switch(nsides) {//ncorners
      case 3://linear : 2
         nfaces = 0;
         face0 = -1;//=side
         break;
      case 7://tri�ngulo : 3
         nfaces = 1;
         face0 =  6;//=side
         break;
      case 9://quadril�tero : 4
         nfaces = 1;
         face0 =  8;//=side
         break;
      case 15://tetraedro : 4
         nfaces = 4;
         face0 = 10;//=side
         break;
      case 19://pir�mide : 5
         nfaces = 5;
         face0 = 13;//=side
         break;
      case 21://prisma : 6
         nfaces = 5;
         face0 = 15;//=side
         break;
      case 27://cubo : 8
         nfaces = 6;
         face0 = 20;//=side
         break;
      default:
         cout << "\nmain::FirstFace elemento nao conhecido\n";
   }
}

void ElementName(TPZCompEl *cel,ofstream &out) {

  switch(cel->Type()){
      case 0:
      	out << "Ponto";
         break;
      case 1:
         out << "Linear";
         break;
      case 2:
      	out << "Triangulo";
         break;
      case 3:
      	out << "Quadrilatero";
         break;
      case 4:
      	out << "Tetraedro";
         break;
      case 5:
      	out << "Piramide";
         break;
      case 6:
      	out << "Prisma";
         break;
      case 7:
         out << "Cubo";
         break;
      default:
      	cout << "Nao tratado";
   }
}

void LerMalha(const char *malha,TPZGeoMesh *geomesh,TPZCompMesh *&compmesh) {

   ifstream grid(malha);
   TPZFMatrix xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
   TPZAutoPointer<TPZMaterial> mat1d;
   int matindex = 1;
   TPZMaterialTest3D *mattest = new TPZMaterialTest3D(1);
   mattest->SetMaterial(xk);
   mat1d = mattest;
   int nnode,nel,ncorners;
   TPZVec<REAL> coord(3);
   grid >> nnode >> nel;
   geomesh->NodeVec().Resize(nnode);
   for(int inode=0;inode<nnode;inode++) {
      grid >> coord[0];
      grid >> coord[1];
      grid >> coord[2];
      geomesh->NodeVec()[inode].Initialize(coord,*geomesh);
   }
   compmesh = new TPZCompMesh(geomesh);
   compmesh->InsertMaterialObject(mat1d);
   for(int el=0;el<nel;el++)  {
      int ntype;
      grid >> ntype;
      ncorners = ntype;
      if(ntype == 7) ncorners = 4;
      TPZVec<int> nodes(ncorners);
      int index;
      for(int incid=0;incid<ncorners;incid++) grid >> nodes[incid];
      geomesh->CreateGeoElement((MElementType) ntype,nodes,matindex,index,1);
   }
      //montagem de conectividades entre elementos
   geomesh->BuildConnectivity();

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

ofstream divide("divide.txt");
void AutomaticDivide(TPZCompMesh &cmesh,int actuallevel) {

   int maxlevel;
   //cout << "\nNivel maximo a ser dividido : ";
   //cin >> maxlevel;
   maxlevel = actuallevel;
   TPZVec<int> csub;
   int el,elid=-1,numeldiv=0;
   TPZCompEl *cpel;
   TPZGeoEl *gel;
   TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh.ElementVec();
   int level = 0,nelc;
   static int nao = -1;
   if(nao == -1){
      cout << "\nInterpolar a divisao (0/1)? ";
      cin >> nao;
   }
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
      divide << "\nnivel dividido " << level << endl;
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
