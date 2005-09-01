//#include "pzmetis.h"
//#include "pztrnsform.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgnode.h"
#include "pznonlinanalysis.h"
#include "pztempmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pzsolve.h"
#include "pzintel.h"
#include "pzelg1d.h"
#include "pzelc1d.h"
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
#include "pzmathyperelastic.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream.h>
struct TPZGeoElBC;


//#define NOTDEBUG

void force(TPZVec<REAL> &point, TPZVec<REAL> &f);
void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f);
void BCNeumann(TPZVec<REAL> &point, TPZVec<REAL> &f);
void ExactSol(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv);
void PostProcess(TPZGeoMesh &gmesh,ostream &out);
void CC(TPZGeoMesh *gmesh);
void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *& compmesh);
void Divide(TPZCompMesh *compmesh);
void AutomaticDivide(TPZCompMesh &cmesh,int maxlevel);
TPZMaterial *mat;
using namespace std;
int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   //malha computacional
   TPZCompMesh *secondmesh;
   LerMalha("malha.txt",firstmesh,secondmesh);
   //Arquivos de saida
	ofstream outgm("outgm.dat");
   ofstream outcm("outcm.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
   //malha computacional
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");
   //forcingfunction
   int nmat = secondmesh->MaterialVec().NElements(),k;
   for(k=0;k<nmat;k++) {
      mat = secondmesh->MaterialVec()[k];
      mat->SetForcingFunction(force);
   }//para diferentes tipos de elementos o material pode ser o mesmo
   //CC : condicões de contorno
   CC(firstmesh);
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,... : ";
   cin >> ord;
   TPZCompEl::gOrder = ord;
   //construção malha computacional
   secondmesh->AutoBuild();
   //Divisão
   cout << "\nDivide niveis completos (0/1)? ";
   cin >> ord;
   if(!ord){
      cout << "\nDivide elementos (0/1)? ";
      cin >> ord;
      if(ord) Divide(secondmesh);
   } else {
      cout << "\nEntre nivel maximo a ser dividido : ";
      int maxlevel;
      cin >> maxlevel;
   	AutomaticDivide(*secondmesh,maxlevel);
	}
   //imprime malha geometrica antes de AdjustBoundaryElement
   firstmesh->Print(outgm);
   outgm.flush();
   secondmesh->Print(outcm);
   outcm.flush();
   secondmesh->AdjustBoundaryElements();
   //imprime malha geometrica depois de AdjustBoundaryElement
   firstmesh->Print(outgm);
   outgm.flush();
   secondmesh->Print(outcm);
   outcm.flush();
   //analysis
	secondmesh->InitializeBlock();
   TPZHyperElasticAnalysis *an = new TPZHyperElasticAnalysis(secondmesh,outcm);
   int numeq = secondmesh->NEquations();
   //TPZVec<int> skyline;
   //secondmesh->Skyline(skyline);
	//TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
   TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
   an->SetMatrix(stiff);
   an->Solver().SetDirect(ELU);//ELDLt , ECholesky
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   /*TPZVec<char *> scalnames(2);
   scalnames[0] = "Solution";
   scalnames[1] = "Derivate";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "plot.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);*/
   int numiter;
   REAL tol;
   cout << "\nNumero de iteracoes? : ";
   cin >> numiter;
   cout << "\nTolerancia? : 1.0e-15";
   //cin >> tol;
   tol = 1.0e-15;
   an->SetExact(ExactSol);
   //TPZSkylMatrix *stifffix = new TPZSkylMatrix(numeq,skyline);
   an->IterativeProcess(outcm,tol,numiter);//,*stifffix
   PostProcess(*firstmesh,outcm);
   int fim;
   cout << "\n~>Fim ";
   cin >> fim;
   delete an;
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN
static int mattype1=-1,mattype2=-1,key = 1;
static REAL pi = 2.*asin(1.);
static REAL tetha;
//solução exata e sua derivada
void ExactSol(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv) {

   if(mattype2==-1) {
   	cout << "\nEntre BC material (0/1/2/3/4) :";
      cin >> mattype2;
      cout << endl;
   }
   val[0] = 0.0;
   val[1] = 0.0;
   val[2] = 0.0;
   deriv.Zero();
   if(mattype2 == 1){//mu=702, lambda=400 <=>  a=100, b=-1602, c=-1/2
      REAL fact = sqrt(2.0)-1.0;
      val[0] = fact*point[0];
      val[1] = fact*point[1];
      val[2] = fact*point[2];
      deriv(0,0) = fact;//du/dx
      deriv(1,1) = fact;//dv/dy
      deriv(2,2) = fact;//dw/dz
   }else
   if(mattype2 == 2){//a = 1 , b = -2 <=> mu = 0 , lambda = 4
      val[0] = 0.5*point[1]*point[1];
      deriv(1,0) = point[1];//du/dy
   }else
   if(mattype2 == 3){//a = 1 , b = -2 <=> mu = 0 , lambda = 4
      val[0] = 0.5*point[0]*point[0];
      deriv(0,0) = point[0];//du/dy
   }else
   if(mattype2==4){//mu=3/2, lambda=4 => a=1, b=-7/2, c=1/2,
      val[0] = point[0]*(cos(tetha)-1.) - point[2]*sin(tetha);
      val[1] = 0.0;
      val[2] = point[0]*sin(tetha) + point[2]*(cos(tetha)-1.);
      deriv(0,0) = -1. + cos(tetha);//du/dx
      deriv(2,0) = -sin(tetha);//du/dz
      deriv(0,2) =  sin(tetha);//dw/dx
      deriv(2,2) = -1. + cos(tetha);//dw/dz
   }
}

void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f) {

   TPZVec<REAL> normal(3,0.);
   if(mattype1==-1) {//mattype1
   	cout << "\nEntre BC material (0/1/2/3/4) :";
      cin >> mattype1;
      cout << endl;
      mattype2 = mattype1;
   }
   f[0] = 0.0;//=u1
   f[1] = 0.0;//=u2
   f[2] = 0.0;//=u3
   if(mattype1==1){//mu=702, lambda=400
      REAL fact = sqrt(2.0)-1.0;
      f[0] = fact*point[0];
      f[1] = fact*point[1];
      f[2] = fact*point[2];
   }else
   if(mattype1 == 2){//a = 1 , b = -2 <=> mu = 0 , lambda = 4
      f[0] = 0.5*point[1]*point[1];
   } else
   if(mattype1 == 3){//a = -c , b = 0, c = -1/2 <=> mu = 1 , lambda = -2
      f[0] = 0.5*point[0]*point[0];
   } else
   if(mattype1==4){//mu=3/2, lambda=4 => a=1, b=-7/2, c=1/2,
      if(key==1) {
         REAL grau;
         cout << "\nEntre grau em (0,90) : ";
         cin >> grau;
         tetha = grau*pi/180.;
         key = 0;
      }
      f[0] = point[0]*(cos(tetha)-1.) - point[2]*sin(tetha);
      f[1] = 0.0;
      f[2] = point[0]*sin(tetha) + point[2]*(cos(tetha)-1.);
   }
   TPZHyperElasticAnalysis han;
   han.BoundCondStep_n(point,0,f,normal);
}

void BCNeumann(TPZVec<REAL> &point, TPZVec<REAL> &f) {

   TPZHyperElasticAnalysis han;
   f[0] =  0.0;//=du1/dn  : U(x,y,z) = (u1,u2,u3)
   f[1] =  0.0;//=du2/dn          ui = ui(x,y,z)
   f[2] =  0.0;//=du3/dn      dui/dn = n1*dui/dx+n2*dui/dy+n3*dui/dz
   REAL fact;
   TPZVec<REAL> normal(3,0.);
   int var=0;
   if(mattype1==1) fact = sqrt(2.0)-1.0;//mattype = 1
   if(point[0] == 1.){//x = 1 : F22 do cubo unitario :  d/dx
      normal[0] = 1.0;
      var = 1;//f[0] : du1/dn
	   if(mattype1==1) f[0] =  fact;//du1/dx : n = (1,0,0)
      if(mattype1==3) f[0] =  point[0];//du1/dx
      if(mattype1==4){
      	f[0] =  cos(tetha);//du1/dx
         han.BoundCondStep_n(point,1,f,normal);
			f[2] =  sin(tetha);//du3/dx
         han.BoundCondStep_n(point,3,f,normal);
         return;
      }
   }else
   if(point[0] == 0.){//x = 0 : F24 do cubo unitario : -d/dx
      normal[0] = -1.0;
      var = 1;//f[0] : du1/dn
	   if(mattype1==1) f[0] =  -fact;//-du1/dx  : n = (-1,0,0)
      if(mattype1==3) f[0] =  -point[0];//-du1/dx
      if(mattype1==4){
      	f[0] = -cos(tetha);//-du1/dx
         han.BoundCondStep_n(point,1,f,normal);
			f[2] = -sin(tetha);//-du3/dx
         han.BoundCondStep_n(point,3,f,normal);
         return;
      }
   }else
   if(point[1] == 1.){//y = 1 : F21 do cubo unitario :  d/dy
      normal[1] = 1.0;
	   if(mattype1==1) {
         var = 2;//f[1] : du2/dn
      	f[1] =  fact;//du2/dy : n = (0,1,0)
      }
      if(mattype1==2) {
         var = 1;//f[0] : du1/dn
      	f[0] =  point[1];//du1/dy : n = (0,1,0)
      }
   }else
   if(point[1] == 0.){//y = 0 : F23 do cubo unitario : -d/dy :
      normal[1] = -1.0;
	   if(mattype1==1) {
         var = 2;//f[1] : du2/dn
      	f[1] =  -fact;//-du2/dy : n = (0,-1,0)
      }
      if(mattype1==2) {
      	var = 1;//f[0] : du1/dn
         f[0] =  -point[1];//-du1/dy
      }
   }else
   if(point[2] == 1.){//z = 1 : F25 do cubo unitario : d/dz
      normal[2] = 1.0;
      var = 3;//f[2] : du3/dn
	   if(mattype1==1) f[2] =  fact;//du3/dz : n = (0,0,1)
      if(mattype1==4){
      	f[0] = -sin(tetha);//du1/dz
         han.BoundCondStep_n(point,1,f,normal);
			f[2] =  cos(tetha);//du3/dz
         han.BoundCondStep_n(point,3,f,normal);
         return;
      }
   }else
   if(point[2] == 0.){//z = 0 : F20 do cubo unitario :  -d/dz
      normal[2] = -1.0;
      var = 3;//f[2] : du3/dn
	   if(mattype1==1) f[2] =  -fact;//-du3/dz : n = (0,0,-1)
      if(mattype1==4){
      	f[0] =  sin(tetha);//-du1/dz
         han.BoundCondStep_n(point,1,f,normal);
			f[2] = -cos(tetha);//-du3/dz
         han.BoundCondStep_n(point,3,f,normal);
         return;
      }
   }
   han.BoundCondStep_n(point,var,f,normal);
}

//força de corpo : b0
void force(TPZVec<REAL> &/*point*/, TPZVec<REAL> &f) {
   f[0] =  .0;//força
   f[1] =  .0;//de corpo
   f[2] =  .0;//nula
}

void CC(TPZGeoMesh *gmesh) {
   TPZBndCond *bc;
   REAL big  = 1.e12;
   TPZFMatrix val1(3,1,big), val2(3,1,0.);
   int side,ncc,type,elid,matindex;
   ifstream in("ccon.txt");
   in >> elid >> matindex >> ncc;
   int idc=0,el;
   int nelg = gmesh->ElementVec().NElements();
   while(elid > -1) {
      TPZGeoEl *elgi;
      for(el=0;el<nelg;el++) {
         elgi = gmesh->ElementVec()[el];
         if(elgi && elgi->Id() == elid) break;
      }
      mat = gmesh->Reference()->MaterialVec()[matindex];
      if(elgi) {
         for(int k=0;k<ncc;k++) {
            in >> side >> type;
            TPZGeoElBC(elgi,side,--idc,*gmesh);
            bc = mat->CreateBC(idc,type,val1,val2);
            gmesh->Reference()->InsertMaterialObject(bc);
            if(type==0) bc->SetForcingFunction(BCDirichlet);
            if(type==1) bc->SetForcingFunction(BCNeumann);
         }
      }
      in >> elid >> matindex >> ncc;
   }
}

void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *&compmesh) {

	ifstream grid(malha);
   //material
   REAL e,mu,nu,lambda,coef1,coef2,coef3;
   grid >> e >> mu >> nu >> lambda >> coef1 >> coef2 >> coef3;
   TPZFMatrix xk(3,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
   TPZMat1dLin *mat1d;
   TPZMat2dLin *mat2d;
   TPZMaterialTest3D *mat3d;
   TPZMatHyperElastic *mat4d;
   int nnode,nel,ncorners;//,yaxisindex;
   //TPZGeoEl1d *elg1d;TPZGeoElT2d *elgt2d;TPZGeoElQ2d *elgq2d;TPZGeoElC3d *elgc3d;TPZGeoElT3d *elgt3d;TPZGeoElPi3d *elgpi3d;TPZGeoElPr3d *elgpr3d
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
   for(int el=0;el<nel;el++)  {
      int mat,ntype;
      grid >> ntype >> mat;
      ncorners = ntype;
      if(ntype == 7) ncorners = 4;
      TPZVec<int> nodes(ncorners);
      for(int incid=0;incid<ncorners;incid++) grid >> nodes[incid];
      switch(ntype) {//tipo de elemento
         case 2://unidimensional ; elg1d =
            new TPZGeoEl1d(nodes,mat,*geomesh);
            if(mat==1) {
               mat1d = new TPZMat1dLin(mat);
               compmesh->InsertMaterialObject(mat1d);
               if(mat==1) mat1d->SetMaterial(xk,xb,xc,xf);
            }
            if(mat==2) {
            	mat2d = new TPZMat2dLin(mat);
               compmesh->InsertMaterialObject(mat2d);
               mat2d->SetMaterial(xk,xc,xf);
            }
            break;
         case 3://triângulo ; elgt2d =
         	new TPZGeoElT2d(nodes,mat,*geomesh);
            if(mat==2) {
            	mat2d = new TPZMat2dLin(mat);
               compmesh->InsertMaterialObject(mat2d);
               mat2d->SetMaterial(xk,xc,xf);
            }
            if(mat==3) {
            	mat3d = new TPZMaterialTest3D(mat);
               compmesh->InsertMaterialObject(mat3d);
               mat3d->SetMaterial(xk);
            }
            break;
         case 4://quadrilátero ; elgq2d =
            new TPZGeoElQ2d(nodes,mat,*geomesh);
            if(mat==2) {
            	mat2d = new TPZMat2dLin(mat);
               compmesh->InsertMaterialObject(mat2d);
               mat2d->SetMaterial(xk,xc,xf);
            }
            if(mat==3) {
            	mat3d = new TPZMaterialTest3D(mat);
               compmesh->InsertMaterialObject(mat3d);
               mat3d->SetMaterial(xk);
            }
            break;
         case 7://tetraedro ; elgt3d =
            new TPZGeoElT3d(nodes,mat,*geomesh);
            mat4d = new TPZMatHyperElastic(mat,e,mu,nu,lambda,coef1,coef2,coef3);
            compmesh->InsertMaterialObject(mat4d);
            mat4d->SetMaterial(xk);
            break;
         case 5://pirâmide ; elgpi3d =
            new TPZGeoElPi3d(nodes,mat,*geomesh);
            mat4d = new TPZMatHyperElastic(mat,e,mu,nu,lambda,coef1,coef2,coef3);
				compmesh->InsertMaterialObject(mat4d);
            mat4d->SetMaterial(xk);
            break;
         case 6://prisma ; elgpr3d =
            new TPZGeoElPr3d(nodes,mat,*geomesh);
            mat4d = new TPZMatHyperElastic(mat,e,mu,nu,lambda,coef1,coef2,coef3);
				compmesh->InsertMaterialObject(mat4d);
            mat4d->SetMaterial(xk);
            break;
         case 8://cubo ; elgc3d =
            new TPZGeoElC3d(nodes,mat,*geomesh);
            mat4d = new TPZMatHyperElastic(mat,e,mu,nu,lambda,coef1,coef2,coef3);
				compmesh->InsertMaterialObject(mat4d);
            mat4d->SetMaterial(xk);
	         break;
         default:
         	for(int i=0;i<300;i++)
            	cout << "\nmain::LerMalha -> Elemento nao conhecido\n";
            cout << "\nChao\n";
            exit(1);
            delete mat1d;
            delete mat2d;
            delete mat3d;
            delete mat4d;
      }
   }
   grid.close();
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
      //cout << "\nNao interpola ! : Divide(el,csub,0)\n";
      compmesh->Divide(el,csub,0);
      n1 = 1;
   }
}

ofstream divide("divide.dat");
void AutomaticDivide(TPZCompMesh &cmesh,int maxlevel) {

   TPZVec<int> csub;
   int el,elid,numeldiv=0;
   TPZCompEl *cpel;
   TPZGeoEl *gel;
   TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh.ElementVec();
   int level = 0,nelc;
   //cout << "\nInterpolar a divisao (0/1)? ";
   int nao = 0;
   //cin >> nao;
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

void PostProcess(TPZGeoMesh &gmesh,ostream &out) {

 int numpoints;
 int nel = gmesh.Reference()->ElementVec().NElements();
 for(int iel=0;iel<nel;iel++) {
 	if(!gmesh.Reference()->ElementVec()[iel]) continue;
 	int elemtype = gmesh.Reference()->ElementVec()[iel]->Type();
   if(elemtype==0) continue;
   TPZCompEl *el = gmesh.Reference()->ElementVec()[iel];
   TPZGeoEl *gel = el->Reference();
   if(el && gel) {
      TPZGeoEl1d  *el1d=0;
      TPZGeoElT2d *elt2d=0;
      TPZGeoElQ2d *elq2d=0;
      TPZGeoElT3d *elt3d=0;
      TPZGeoElPi3d *elpi3d=0;
      TPZGeoElPr3d *elpr3d=0;
      TPZGeoElC3d  *elc3d=0;
      if(elemtype==1) el1d   = (TPZGeoEl1d    *) gel;
      if(elemtype==2) elt2d  = (TPZGeoElT2d   *) gel;
      if(elemtype==3) elq2d  = (TPZGeoElQ2d   *) gel;
      if(elemtype==4) elt3d  = (TPZGeoElT3d   *) gel;
      if(elemtype==5) elpi3d = (TPZGeoElPi3d  *) gel;
      if(elemtype==6) elpr3d = (TPZGeoElPr3d  *) gel;
      if(elemtype==7) elc3d  = (TPZGeoElC3d   *) gel;
   	out << "Elemento " << el->Reference()->Id() << endl;;
   	TPZManVector<REAL> sol(3);
      TPZVec<REAL> csi(3,0.),x(3);
      ifstream in("pontos.txt");
      in >> numpoints;
      for(int p=0;p<numpoints;p++) {
         in >> csi[0] >> csi[1] >>  csi[2];
         gel->X(csi,x);
         if(elemtype==1) el1d->Reference()->Solution(csi,0,sol);
         if(elemtype==2) elt2d->Reference()->Solution(csi,0,sol);
         if(elemtype==3) elq2d->Reference()->Solution(csi,0,sol);
         if(elemtype==4) elt3d->Reference()->Solution(csi,0,sol);
         if(elemtype==5) elpi3d->Reference()->Solution(csi,0,sol);
         if(elemtype==6) elpr3d->Reference()->Solution(csi,0,sol);
         if(elemtype==7) elc3d->Reference()->Solution(csi,0,sol);
         out << "solucao em x y z  = " << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
         out << "               u  = " << sol[0] << endl;
         out << "               v  = " << sol[1] << endl;
         out << "               w  = " << sol[2] << endl;
      }
      in.close();
   }
 }
}
/*
   //comparação dos connects de lados viz da malha
   int i,side;
   TPZCompEl *cel;
   for (i=0;i<secondmesh->NElements();i++){
      cel = secondmesh->ElementVec()[i];
      if(!cel) continue;
      if(!cel->Reference()){
         cout << "\nReferencia NULA ?!\n";
         continue;
      }
      TPZCompElSide vizcelside;
      TPZGeoElSide vizgelside;
      for(side=0;side<cel->NConnects();side++){
         TPZCompElSide celside(cel,side);
         TPZGeoElSide gelside = celside.Reference();
         vizgelside = gelside.Neighbour();
         if(!vizgelside.Element()){
         	//cout << "\nVizinho NULO ?!\n";
            //cout << "\nElemento/lado : " << gelside.Element()->Id() << "/" << gelside.Side();
            continue;
         }
         vizcelside = vizgelside.Reference();
         if(!vizcelside.Element()){
         	//cout << "\nReferencia NULA ?!\n";
            continue;
         }
         while(vizcelside.Element() != celside.Element()){
            if(celside.Element()->ConnectIndex(side) != vizcelside.Element()->ConnectIndex(vizcelside.Side()))
            	cout << "\nErro de vizinhanza achado\n";
            vizgelside = vizcelside.Reference();
            do{
               vizgelside = vizgelside.Neighbour();
               if(!vizgelside.Element()){
                  cout << "\nVizinho NULO ?!\n";
               }
            } while(!vizgelside.Reference().Element());
            vizcelside = vizgelside.Reference();
         }
      }
   }
*/
