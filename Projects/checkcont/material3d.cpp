//#include "pzmetis.h"
//#include "pztrnsform.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgnode.h"
#include "pzanalysis.h"
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
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream.h>
struct TPZGeoElBC;
//#define NOTDEBUG

int eq;//tipo de equação
REAL p,q,r;//expoentes da solução exata  u =  xp*yq*zr
void force(TPZVec<REAL> &point, TPZVec<REAL> &f);
void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f);
void BCNeumann(TPZVec<REAL> &point, TPZVec<REAL> &f);
void ExactSol(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv);
void InverseX(TPZVec<REAL> x,TPZVec<REAL> &sol);
void PostProcess(TPZGeoMesh &gmesh,ostream &out);
void CC(TPZGeoMesh *firstmesh);
void NormalVectorLin(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorTr(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorQ(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorTe(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorPi(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorPr(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void NormalVectorC(REAL **n,REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s);
void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *& compmesh);
void LocalAxes(TPZVec<REAL> point,REAL &x,REAL &y,REAL &z);
void Divide(TPZCompMesh *compmesh);
void AutomaticDivide(TPZCompMesh &cmesh,int maxlevel);
TPZMaterial *mat;

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
   //TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");
   //forcingfunction
   int nmat = secondmesh->MaterialVec().NElements(),k;
   for(k=0;k<nmat;k++) {
      mat = secondmesh->MaterialVec()[k];
      mat->SetForcingFunction(force);
   	//if(mat->Id()==1) mat->SetForcingFunction(force);
      //if(mat->Id()==2) mat->SetForcingFunction(force);
      //etc.
   }//para diferentes tipos de elementos o material pode ser o mesmo
   //escolha da equacao a discretizar
   cout << "0 : Projecao L2 ; 1 : Equacao de Poisson\n";
   cin >> eq;
 //  TPZMat1dLin::eq1 = eq;//caso o material
 //  TPZMat2dLin::eq2 = eq;//correspondente
   TPZMaterialTest3D::eq3 = eq;//foi considerado
   cout << "Expoentes positivos p,q,r da solucao exata  u =  xp*yq*zr\n";
   cin >> p >> q >> r;
   //CC : condicões de contorno
   CC(firstmesh);
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,... : ";
   cin >> ord;
   TPZCompEl::gOrder = ord;
   //construção malha computacional
   secondmesh->AutoBuild();
   //designando o axisindex
   /*int nelc = secondmesh->ElementVec().NElements();
   int yaxisindex;
   for(k=0;k<nelc;k++) {
      TPZCompEl *cel = secondmesh->ElementVec()[k];
      if(!cel) continue;
      if(cel->Dimension()==2) break;
      if(cel->Reference()->NSides()==3) {
         TPZCompEl1d *c1d = (TPZCompEl1d *) cel;
         TPZGeoEl1d *g1d = (TPZGeoEl1d *) c1d->Reference();
         if(g1d->fYAxisIndex!=-1) continue;
         cout << "\nEntre index de yaxis : ";
         cin >> yaxisindex;
      	g1d->SetAxisIndex(yaxisindex);
      }
   }*/
   //Divisão
   cout << "\nDivide (0/1)?\n";
   cin >> ord;
   if(ord) Divide(secondmesh);
   firstmesh->Print(outgm);
   outgm.flush();
   secondmesh->Print(outcm);
   outcm.flush();
   secondmesh->AdjustBoundaryElements();
   //analysis
	secondmesh->InitializeBlock();
   TPZAnalysis an(secondmesh,outcm);
   int numeq = secondmesh->NEquations();
   TPZVec<int> skyline;
   secondmesh->Skyline(skyline);
	TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
   //TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
   an.SetMatrix(stiff);
   an.Solver().SetDirect(ELDLt);//ELU , ECholesky
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   /*TPZVec<char *> scalnames(2);
   scalnames[0] = "Solution";
   scalnames[1] = "Derivate";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "plot.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);*/
   an.Run(outcm);
   an.Print("FEM SOLUTION ",outcm);
   firstmesh->Print(outgm);
   outgm.flush();
   //erro da aproximacao FEM
   REAL estimate=0.,H1_error=0.,L2_error=0.;
   REAL fluxstore[9];
   TPZManVector<REAL> flux(0,fluxstore,9);
   secondmesh->EvaluateError(ExactSol,H1_error,L2_error,estimate);
	outcm << "\n\n\n";
   outcm << "Energy    : " << H1_error << endl;
   outcm << "L2 Error  : " << L2_error << endl;
   outcm << "Semi Norm : " << estimate << endl;
   outcm.flush();
   PostProcess(*firstmesh,outcm);
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN
void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *&compmesh) {

	ifstream grid(malha);
   TPZFMatrix xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
   TPZMat1dLin *mat1d;
   TPZMat2dLin *mat2d;
   TPZMaterialTest3D *mat3d;
   int nnode,nel,ncorners,yaxisindex;
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
            //cout << "\nfYAxisIndex (n/-1)? : ";
            //cin >> yaxes;
            grid  >> yaxisindex;//no final do arquivo malh.txt
            new TPZGeoEl1d(nodes,mat,*geomesh,yaxisindex);
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
            if(mat==3) {
            	mat3d = new TPZMaterialTest3D(mat);
               compmesh->InsertMaterialObject(mat3d);
               mat3d->SetMaterial(xk);
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
            mat3d = new TPZMaterialTest3D(mat);
            compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 5://pirâmide ; elgpi3d =
            new TPZGeoElPi3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
				compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 6://prisma ; elgpr3d =
            new TPZGeoElPr3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
				compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 8://cubo ; elgc3d =
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
      cout << "\nNao interpola ! : Divide(el,csub,0)\n";
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

//solução exata e sua derivada
void ExactSol(TPZVec<REAL> &point,TPZVec<REAL> &val,TPZFMatrix &deriv) {
   if(p<0. || q<0. || r<0.) {
   	cout << "Expoente nao permitido\n";
   	exit(-1);
   }
   REAL x = point[0];
   REAL y = point[1];
   REAL z = point[2];
   LocalAxes(point,x,y,z);
   val[0]  = 0.;
   deriv(0,0) = 0.;
   deriv(1,0) = 0.;
   deriv(2,0) = 0.;
   BCDirichlet(point,val);
   if(p) {
   	if(p==1) deriv(0,0)  = 1.;
      else     deriv(0,0)  = p*pow(x,p-1.);
      if(q)    deriv(0,0) *= pow(y,q);
      if(r)    deriv(0,0) *= pow(z,r);
   }
   if(q) {
   	if(q==1) deriv(1,0)  = 1.;
      else     deriv(1,0)  = q*pow(y,q-1.);
      if(p)    deriv(1,0) *= pow(x,p);
      if(r)    deriv(1,0) *= pow(z,r);
   }
   if(r) {
   	if(r==1) deriv(2,0)  = 1.;
      else     deriv(2,0)  = r*pow(z,r-1.);
      if(p)    deriv(2,0) *= pow(x,p);
      if(q)    deriv(2,0) *= pow(y,q);
   }
   /////////////////////////////////////
   //REAL derivv1 = deriv(1,0);
   //REAL derivv2 = -.7071067811865*deriv(0,0);
   //deriv(0,0) = derivv1;
   //deriv(1,0) = derivv2;
}

//carga ou fonte
void force(TPZVec<REAL> &point, TPZVec<REAL> &f) {
   //p,q,r > 0
   f[0] = 0.;
   if(p<1. && q<1. && r<1.) {
   	cout << "Expoente nao permitido\n";
      exit(-1);
   }
   REAL x = point[0];
   REAL y = point[1];
   REAL z = point[2];
   LocalAxes(point,x,y,z);
   if(eq==0) {//projeção L2 : solução = u = fonte = f =
      BCDirichlet(point,f);
   } else if(eq==1) {//Poisson : solução u : f = -Lap(u)
      REAL d1=0.,d2=0.,d3=0.;
      if(p>1.) {
         if(p==2) d1  = -2.;
         else     d1  = -p*(p-1.)*pow(x,p-2.);
         if(q)    d1 *= pow(y,q);
         if(r)    d1 *= pow(z,r);
      }
      if(q>1.) {
         if(q==2) d2  = -2.;
         else     d2  = -q*(q-1.)*pow(y,q-2.);
         if(p)    d2 *= pow(x,p);
         if(r)    d2 *= pow(z,r);
      }
      if(r>1.) {
         if(r==2) d3  = -2.;
         else     d3  = -r*(r-1.)*pow(z,r-2.);
         if(p)    d3 *= pow(x,p);
         if(q)    d3 *= pow(y,q);
      }
      f[0] = d1+d2+d3;
      ////////////////////////////////
      //f[0] = .5*d1 + d2;
   }
}

void LocalAxes(TPZVec<REAL> point,REAL &x,REAL &y,REAL &z) {

   static int key=1;
   static REAL axis[4][3];
   if(key==1) {
      ifstream axes("axes.txt");
   	for(int i=0;i<4;i++) for(int j=0;j<3;j++) axes >> axis[i][j];
      key = 0;
      axes.close();
   }
   x = axis[0][0]*point[0]+axis[0][1]*point[1]+axis[0][2]*point[2]+axis[0][3];
   y = axis[1][0]*point[0]+axis[1][1]*point[1]+axis[1][2]*point[2]+axis[1][3];
   z = axis[2][0]*point[0]+axis[2][1]*point[1]+axis[2][2]*point[2]+axis[2][3];
}

void BCDirichlet(TPZVec<REAL> &point, TPZVec<REAL> &f) {
   REAL x = point[0];
   REAL y = point[1];
   REAL z = point[2];
   LocalAxes(point,x,y,z);
   f[0] = 0.;
   if(p)         f[0]  = pow(x,p);
   if(q) {
      if(!p)     f[0]  = pow(y,q);
      else       f[0] *= pow(y,q);
   }
   if(r) {
      if(p || q) f[0] *= pow(z,r);
      else       f[0]  = pow(z,r);
   }
}

void NormalVectorLin(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {
   if(fabs(c+1.)<1.e-10 && fabs(e)<1.e-10 && fabs(s)<1.e-10) n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
   if(fabs(c-1.)<1.e-10 && fabs(e)<1.e-10 && fabs(s)<1.e-10) n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
   if(fabs(c)<1.e-10 && fabs(e-1.)<1.e-10 && fabs(s)<1.e-10) n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
}


void NormalVectorTr(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {
//triangulo : c:csi , e:eta ; s:zeta
   if (fabs(e) < 1.e-12) {//aresta 3
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(e+c-1.) < 1.e-12) {//aresta 4
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(c) < 1.e-12) {//aresta 5
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(s) < 1.e-12) {//face 6
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   }
   cout << "\nmain::NormalVectorTr error\n";
}

void NormalVectorQ(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {

   if (fabs(e+1.) < 1.e-12) {//aresta 4
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(c-1.) < 1.e-12) {//aresta 5
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(e-1.) < 1.e-12) {//aresta 6
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(c+1.) < 1.e-12) {//aresta 7
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   } else
   if (fabs(s) < 1.e-12) {//face 8
      n1 = n[0][4]; n2 = n[1][4]; n3 = n[2][4];
      return;
   }
   cout << "\nmain::NormalVectorQ error\n";
}

void NormalVectorPi(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {

   if (fabs(s) < 1.e-12) {//face 13
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(s-1.-e) < 1.e-12) {//face 14
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(s-1.+c) < 1.e-12) {//face 15
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(s-1.+e) < 1.e-12) {//face 16
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   } else
   if (fabs(s-1.-c) < 1.e-12) {//face 17
      n1 = n[0][4]; n2 = n[1][4]; n3 = n[2][4];
      return;
   }
   cout << "\nmain::NormalVectorPi error\n";
}

void NormalVectorTe(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {

   if (fabs(s) < 1.e-12) {//face 10
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(e) < 1.e-12) {//face 11
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(s-1.+c+e) < 1.e-12) {//face 12
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(c) < 1.e-12) {//face 13
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   }
   cout << "\nmain::NormalVectorTe error\n";
}

void NormalVectorPr(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {

   if (fabs(s+1.) < 1.e-12) {//face 15
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(e) < 1.e-12) {//face 16
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(e-1.+c) < 1.e-12) {//face 17
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(c) < 1.e-12) {//face 18
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   } else
   if (fabs(s-1.) < 1.e-12) {//face 19
      n1 = n[0][4]; n2 = n[1][4]; n3 = n[2][4];
      return;
   }
   cout << "\nmain::NormalVectorPi error\n";
}

void NormalVectorC(REAL n[3][6],REAL &n1,REAL &n2,REAL &n3,REAL c,REAL e,REAL s) {

   if (fabs(s+1.) < 1.e-12) {//face 20
      n1 = n[0][0]; n2 = n[1][0]; n3 = n[2][0];
      return;
   } else
   if (fabs(e+1.) < 1.e-12) {//face 21
      n1 = n[0][1]; n2 = n[1][1]; n3 = n[2][1];
      return;
   } else
   if (fabs(c-1.) < 1.e-12) {//face 22
      n1 = n[0][2]; n2 = n[1][2]; n3 = n[2][2];
      return;
   } else
   if (fabs(e-1.) < 1.e-12) {//face 23
      n1 = n[0][3]; n2 = n[1][3]; n3 = n[2][3];
      return;
   } else
   if (fabs(c+1.) < 1.e-12) {//face 24
      n1 = n[0][4]; n2 = n[1][4]; n3 = n[2][4];
      return;
   } else
   if (fabs(s-1.) < 1.e-12) {//face 25
      n1 = n[0][5]; n2 = n[1][5]; n3 = n[2][5];
      return;
   }
   cout << "\nmain::NormalVectorC error\n";
}

void BCNeumann(TPZVec<REAL> &point, TPZVec<REAL> &f) {

	//point é o ponto deformado
   //pinv é o ponto de integracao
   //o valor da C.C. é a primeira derivada
	REAL pinvstore[3];
	TPZManVector<REAL> pinv(3,pinvstore,3);
   REAL x = point[0];//ponto real
   REAL y = point[1];//ou deformado
   REAL z = point[2];
   LocalAxes(point,x,y,z);//ponto no sistema local
   TPZVec<REAL> local(3);
   local[0] = x;
   local[1] = y;
   local[2] = z;
	InverseX(local,pinv);
   REAL c = pinv[0];//ponto no
   REAL e = pinv[1];//elemento
   REAL s = pinv[2];//de referencia
   f[0] = 0.;//derivada normal
   REAL n1,n2,n3;//cordenadas da normal a face
   int i,j;
   static REAL n[3][6];//maximo de faces para o cubo
   static int key=1,elemtype;
   if(key==1) {
      ifstream normal("normal.txt");
      normal >> elemtype;//tetraedro ou piramide, etc
      int nvec=0;
      if(elemtype==2) nvec = 3;//linear
      if(elemtype==3 || elemtype==7) nvec = 4;//triangulo e tetraedro
      if(elemtype==4 || elemtype==5 || elemtype==6) nvec = 5;//quadrilatero ou piramide ou prisma
      if(elemtype==8) nvec = 6;//cubo
      for(j=0;j<nvec;j++) for(i=0;i<3;i++) normal >> n[i][j];//1 vetor por coluna
      normal.close();
      key = 0;
   }
	if(elemtype==2) NormalVectorLin(n,n1,n2,n3,c,e,s);//linear
	if(elemtype==3) NormalVectorTr(n,n1,n2,n3,c,e,s);//triangulo
	if(elemtype==4) NormalVectorQ(n,n1,n2,n3,c,e,s);//quadrilatero
   if(elemtype==7) NormalVectorTe(n,n1,n2,n3,c,e,s);//tetraedro
	if(elemtype==5) NormalVectorPi(n,n1,n2,n3,c,e,s);//piramide
   if(elemtype==6) NormalVectorPr(n,n1,n2,n3,c,e,s);//prisma
   if(elemtype==8) NormalVectorC(n,n1,n2,n3,c,e,s);//cubo
   REAL d1=0.;
   REAL d2=0.;
   REAL d3=0.;
   if(p) {
   	if(p==1) d1  = n1;
      else     d1  = p*pow(x,p-1.)*n1;
      if(q)    d1 *= pow(y,q);
      if(r)    d1 *= pow(z,r);
   }
   if(q) {
   	if(q==1) d2  = n2;
      else     d2  = q*pow(y,q-1.)*n2;
      if(p)    d2 *= pow(x,p);
      if(r)    d2 *= pow(z,r);
   }
   if(r) {
   	if(r==1) d3  = n3;
      else     d3  = r*pow(z,r-1.)*n3;
      if(p)    d3 *= pow(x,p);
      if(q)    d3 *= pow(y,q);
   }
   f[0] = d1+d2+d3;
}

void CC(TPZGeoMesh *firstmesh) {
   TPZBndCond *bc;
   REAL big  = 1.e12;
   TPZFMatrix val1(1,1,big), val2(1,1,0.);
   int side,ncc,type,elid,matindex;
   ifstream in("ccon.txt");
   in >> elid >> matindex >> ncc;
   int idc=0,el;
   int nelg = firstmesh->ElementVec().NElements();
   while(elid > -1) {
      TPZGeoEl *elgi;
      for(el=0;el<nelg;el++) {
         elgi = firstmesh->ElementVec()[el];
         if(elgi && elgi->Id() == elid) break;
      }
      mat = firstmesh->Reference()->MaterialVec()[matindex];
      if(elgi) {
         for(int k=0;k<ncc;k++) {
            in >> side >> type;
            TPZGeoElBC(elgi,side,--idc,*firstmesh);
            bc = mat->CreateBC(idc,type,val1,val2);
            firstmesh->Reference()->InsertMaterialObject(bc);
            if(type==0) bc->SetForcingFunction(BCDirichlet);
            if(type==1) bc->SetForcingFunction(BCNeumann);
         }
      }
      in >> elid >> matindex >> ncc;
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
   	TPZManVector<REAL> sol(1);
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
         out << "solucao em x    = "<< x[0] << ' ' << x[1] << ' ' << x[2] << endl;
         REAL x1,y1,z1;
         LocalAxes(x,x1,y1,z1);
         out << "solucao local x = "<< x1 << ' ' << y1 << ' ' << z1 << endl;
         out << "              u = "<< sol[0] << endl;
      }
      in.close();
   }
 }
}

ifstream inverse("inverse.txt");//coeficientes da inversa
void InverseX(TPZVec<REAL> p,TPZVec<REAL> &sol) {

   int i,j;
   static REAL T[3][4];
   static int key = 1;
   if(key==1) {
      for(i=0;i<3;i++) for(j=0;j<4;j++) inverse >> T[i][j];
      key = 0;
      inverse.close();
   }
   for(i=0;i<3;i++) sol[i] =  T[i][0]*p[0]+T[i][1]*p[1]+T[i][2]*p[2]+T[i][3];
         //sol é ponto no elemento de referencia
}
/*
      /*if(elemtype==5)
      	for(j=0;j<5;j++) for(i=0;i<3;i++) normal >> n[i][j];//1 vetor por coluna
      if(elemtype==4 || elemtype==3)
      	for(j=0;j<4;j++) for(i=0;i<3;i++) normal >> n[i][j];//1 vetor por coluna
   //designando o axisindex
   int nelbc = firstmesh->BCElementVec().NElements();
   int yaxisindex;
   for(k=0;k<nelbc;k++) {
      TPZGeoElBC *gelbc = &(firstmesh->BCElementVec()[k]);
      if(gelbc->fElement->NNodes()==2) {
         TPZGeoEl1d *g1d = (TPZGeoEl1d *) gelbc->fElement;
         if(g1d->fYAxisIndex!=-1) continue;
         cout << "\nEntre index de yaxis : ";
         cin >> yaxisindex;
      	g1d->SetAxisIndex(yaxisindex);
      }
   }
*/
/*
CC pontual
   /////////////////////////////
   //TPZConnect *cn = &secondmesh->ConnectVec()[8];
   //TPZFMatrix val1(1,1,1.e12),val2(1,1,-2.);
   //TPZBndCond *bnc = secondmesh->MaterialVec()[0]->CreateBC(8,0, val1,val2);
   //secondmesh->BCConnectVec().Resize(1);
   //secondmesh->BCConnectVec()[0] = TPZConnectBC(cn,bnc);//TPZConnectBC(TPZConnect *nd,TPZBndCond *bc)
   /////////////////////////////

*/
