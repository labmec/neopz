#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream.h>

#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzmat2dlin.h"
#include "pzelasmat.h"
#include "pzanalysis.h"
#include "pzmetis.h"

#include <stdio.h>//Cedric
#include <time.h>//Cedric
#include "pzelct2d.h"//Cedric

#define NOTDEBUG
int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh);
void DivideAny(TPZCompMesh &cmsh);
void CoarsenAny(TPZCompMesh &cmsh);
void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel);
void force(TPZVec<REAL> &x, TPZVec<REAL> &f);
void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f);
void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f);

int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : Nós e Elementos");
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;//00
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;//10
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   coord[1] = 1.0;//11
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;//01
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);//triangulo
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;//ex. anterior ; 021
   //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //orientacao local de um segundo elemento superposto
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;//302
   TPZGeoElT2d *elg1 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
  	//teste de divisao geometrica : 1 elemento
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");
   //material
   int matindex;
   TPZFMatrix k(1,1,100.),f(1,1,0.),c(1,2,0.);
   TPZMat2dLin * mat = new TPZMat2dLin(1);
   mat->SetMaterial(k,c,f);
   //TPZElasticityMaterial * mat = new TPZElasticityMaterial(1,1.e5,.25,0.,0.);
   matindex = secondmesh->InsertMaterialObject(mat);
   //forcingfunction
   int carga;
   cout << "Teste tipo funcao/derivada (0/1) : \n";
   //cin >> carga;
   carga = 0;
   if(carga==0) {
   	mat->SetForcingFunction(force);
   } else {//default
   	mat->SetForcingFunction(derivforce);

   }
   //CC : condicao de contorno
   int bcindex;
   cout << "CC (0/1)\n";
   cin >> bcindex;
//	bcindex = 0;
   if(bcindex) {
        TPZGeoElBC(elg1,4,-1,*firstmesh);
        TPZFMatrix val1(1,1,0.), val2(1,1,1.);
        TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
        bcindex = secondmesh->InsertMaterialObject(bc);
        if(carga==0) bc->SetForcingFunction(bcforce);
        else bc->SetForcingFunction(bcforce);
   }
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem geral 1,2,3,4,5 : \n";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   cmesh.SetDefaultOrder(ord);
   //construção malha computacional
   secondmesh->AutoBuild();
   //redistribuicao de ordem aos lados do elemento
	int nel = secondmesh->ElementVec().NElements();
   int intinpose[] = {1,1,1,1,1,1,1,1,1,1};//{4,4,4,4,4,4,4};
   srand(3);
   cout << "Ordem por elemento : \n";
   for(int el=0;el<nel;el++) {
   	cout << "\nElemento : " << el << ' ';
   	cin >> intinpose[el];
   }
   for(int index=0;index<nel;index++) {
   	//int ran = random(3)+2;
      TPZInterpolatedElement *el;
      el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
      el->PRefine(intinpose[index]);
      //el->PRefine(ran);
   }
   int n1,n2;
   cout << "Cyclerefinamnets ? (0/1) \n";
   cin >> n1;
   if(n1) {
        //numeros minimos e maximos de elementos que a malha nao devera ultrapassar apos agrupamento e divisao, respectivamente
        cout << "Parametros do CycleRefinements(n1,n2) n1 < n2 : ";
        cin >> n1 >> n2;
        CycleRefinements(*secondmesh,n1,n2);
   }
   cout << "Divide computacional ? (-1/1) \n";
   cin >> n1;
   while(n1>=0) {
        TPZVec<int> csub;
        cout << "Elemento a dividir : ";
        cin >> n1;
        if(n1<0) break;
        secondmesh->Divide(n1,csub,1);
   }
   //analysis
	secondmesh->InitializeBlock();
//   secondmesh->Print(outcm1);
   TPZAnalysis an(secondmesh,outcm1);
   int numeq = secondmesh->NEquations();
   secondmesh->Print(outcm1);
   firstmesh->Print(outgm1);
   outcm1.flush();
   outgm1.flush();
   TPZVec<int> skyline;
   secondmesh->Skyline(skyline);
   int nskyl = skyline.NElements();
   int sk;
   cout << "Skyline profile";
   for(sk=0; sk<nskyl; sk++) {
   	if((sk%20==0)) cout << endl;
   	cout << skyline[sk] << ' ';
   }
	//TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
   TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
   an.SetMatrix(stiff);
   an.Solver().SetDirect(ELU);
//   an.Solver().SetDirect(ELDLt);
   //an.Solver().SetDirect(ECholesky);
   //an.Solver().SetJacobi(4, 1E-8, 0);
   //an.Solver().SetSOR(4, 0.2, 1E-8, 0);
   //an.Solver().SetSSOR(6, 1.3, 1E-8,0);
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   an.Run(outcm2);

   cout << " '.plt' ou '.pos' ou '.dx'? (0/1/2) : ";
   int arq;
   cin >> arq;
   TPZVec<char *> scalnames(3);
   TPZVec<char *> vecnames(1);
   scalnames[0] = "SigmaX";
   scalnames[1] = "SigmaY";
   scalnames[2] = "TauXY";
   vecnames[0]  = "Displacement";
   char *plotfile =  "plot.pos";
   if(arq==2) plotfile =  "plot.dx";
   if(arq==0) plotfile =  "plot.plt";
   if(arq == 0) {
        if(! strcmp(mat->Name(),"TPZMat2dLin")) {
            scalnames.Resize(1);
            vecnames.Resize(0);
            scalnames[0] = "state";
        }
        plotfile =  "plot.plt";
        an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   } else {
        TPZVec<char *> scalnames(3);
        TPZVec<char *> vecnames(1);
        scalnames[0] = "SigmaX";
        scalnames[1] = "SigmaY";
        scalnames[2] = "TauXY";
        vecnames[0]  = "Displacement";
        an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   }
   cout << "Resolucao : ";
   int res;
   cin >> res;
   an.PostProcess(res);

   //an->Print("Saida para viga.dat\n",out);
   an.Print("FEM SOLUTION ",outcm1);

//   firstmesh->Print(outgm1);
//   outgm1.flush();
///////////////////////////////////////////////////////
   delete secondmesh;
   delete firstmesh;
   return 0;
}
static int init=0,n,m;
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {
   if(init == 0) {
  		init = 1;
      cout<<"\nMaterial::TPZMat2dLin -> Teste com monomios xn*ym. Entre expoentes n,m\n";
      cin>>n>>m;
	}
   int r = f.NElements();
   int ic;
   for(ic=0; ic< r; ic++) {
      if (!n && !m) f[ic]  = 1.;
      if ( n      ) f[ic]  = 1.*pow(x[0], n);
      if ( n &&  m) f[ic] *= pow(x[1], m);
      if (!n &&  m) f[ic]  = 1.*pow(x[1], m);
   }
}

void force(TPZVec<REAL> &x, TPZVec<REAL> &f) {
   if(init == 0) {
  		init = 1;
      cout<<"\nMaterial::TPZMat2dLin -> Teste com monomios xn*ym. Entre expoentes n,m\n";
      cin>>n>>m;
	}
   int r = f.NElements();
   int ic;
   for(ic=0; ic< r; ic++) {
      if (!n && !m) f[ic]  = 100.;
      if ( n      ) f[ic]  = 100.*pow(x[0], n);
      if ( n &&  m) f[ic] *= pow(x[1], m);
      if (!n &&  m) f[ic]  = 100.*pow(x[1], m);
   }
}

void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {
	static int init=0,n,m;
   if(init == 0) {
  		init = 1;
      cout<<"\nMaterial::TPZMat2dLin -> Solucao  f =  xn*ym. Entre expoentes nao negativos n,m \n";
      cin>>n>>m;
	}
   int r = f.NElements();
   int ic;
   for(ic=0; ic< r; ic++) {
   	if(n>0 && m>0) {
      	f[ic]  = pow(x[0], n) * pow(x[1], m);
      	f[ic] += n*pow(x[0], n-1)*pow(x[1], m) + m*pow(x[0], n)*pow(x[1], m-1);
      } else
   	if(n==0 && m>0) {
      	f[ic]  = pow(x[1], m);
      	f[ic] += m*pow(x[1], m-1);
      } else
   	if(n>0 && m==0) {
      	f[ic]  = pow(x[0], n);
      	f[ic] += n*pow(x[0], n-1);
      } else {//outros casos
      	f[ic] = 0;
      }
   }
}

int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh) {
	if(!cel) return 0;
	TPZGeoEl *ref = cel->Reference();
   TPZGeoEl *father = ref->Father();
   if(!father) return 0;
   int numsub = father->NSubElements();
   TPZVec<TPZCompEl *> compsubel(numsub);
   int returncode = 1;
   int is;
   for(is=0; is<numsub; is++) {
		TPZGeoEl *gsub = father->SubElement(is);
      compsubel[is] = gsub->Reference();
      if(!gsub->Reference()) returncode = 0;

   }
   if(!returncode) return returncode;
   TPZManVector<int> subindex;
   int index;
   numsub = compsubel.NElements();
   subindex.Resize(numsub);
   for(is=0; is<numsub; is++) subindex[is] = compsubel[is]->Index();
   if(returncode) {
   	cout << "Coarsening " << compsubel[0]->Reference()->Id() << ' ' << compsubel[1]->Reference()->Id() << endl;
   	mesh.Coarsen(subindex,index);
      cout << "Created " << mesh.ElementVec()[index]->Reference()->Id() << endl;
   }
   return returncode;
}

//ofstream locfile("CycleRefinements");

void DivideAny(TPZCompMesh &cmesh) {
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
   int nelem = vec.NElements();
	int elindex = rand()%nelem;
   TPZCompEl *cel = vec[elindex];
   while(!cel) {
   	elindex = rand()%nelem;
   	cel = vec[elindex];
   }
   TPZManVector<int> subindex;
   cout << "Dividing " << vec[elindex]->Reference()->Id() << endl;
   cout.flush();
//   if(vec[elindex]->Reference()->Id() == 4) {
//	   cmesh.Reference()->Print(locfile);
//	   cmesh.Print(locfile);
//   }

   cmesh.Divide(elindex,subindex);
   cout << " creating ";
   int nsub = subindex.NElements();
   for(int is=0; is<nsub; is++) {
	   cout << vec[subindex[is]]->Reference()->Id() << ' ';
   }
   cout << endl;
   //cmesh.Reference()->Print(locfile);
//   locfile.flush();
//   cmesh.CleanUpUnconnectedNodes();
}

void CoarsenAny(TPZCompMesh &cmesh) {
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
   int nelem = vec.NElements();
	int elindex = rand()%nelem;
   TPZCompEl *cel = vec[elindex];
   while(!cel || !Coarsen(cel,cmesh)) {
		elindex = rand()%nelem;
   	cel = vec[elindex];
   }
//   cmesh.ExpandSolution();
//   cmesh.CleanUpUnconnectedNodes();
}

//ofstream logfile("log.dat");

void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel) {
   int cycle;
   int nelem = cmesh.NElements();
   int i,totalel;
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
	for(cycle = 0; cycle<10; cycle++) {
   	while(nelem < maxel) {
      	DivideAny(cmesh);
         nelem=0;
         totalel = vec.NElements();
         for(i=0; i<totalel; i++) if(vec[i]) nelem++;
      }

      while(nelem > minel) {
      	CoarsenAny(cmesh);
         nelem=0;
         totalel = vec.NElements();
         for(i=0; i<totalel; i++) if(vec[i]) nelem++;
      }
      cmesh.ComputeNodElCon();
//      cmesh.Print(logfile);
//      logfile.flush();
      cmesh.ExpandSolution();
      cmesh.CleanUpUnconnectedNodes();
      cmesh.ExpandSolution();
//      cmesh.Print(logfile);
//      logfile.flush();
//      cmesh.ComputeConnectSequence();
//      if(cycle == 0) cmesh.Print(logfile);
//      logfile.flush();
   }
//   cmesh.Print(logfile);
//   logfile.flush();
}
//____________________PLOT TRIÂNGULO/QUADRILATERO__(Cedric)_____________________
/*void PrintPlot(char *pltname,TPZGeoMesh &gm) {//Gera o arquivo para o VIEW3D
	ofstream out(pltname);
	out<<"Arquivo gerado por Cedric Marcelo Augusto Ayala Bravo\n";
	out<<"2 Dimensional grid\n";
	out<<"dim 3\n";
   int nel = gm.;

   Pix last = MG.NodeMap().last(); //A função atual não pertence a TGeoGrid. O objeto deve ser especificado
   TGeoNod *nglast = (TGeoNod *) MG.NodeMap().contents(last);
   long NodeLast = nglast->Id();
	out<<"coor "<<(NodeLast+1)<< endl;
   int j=0;
   Pix i=MG.NodeMap().first();
   while (i) 			{
    	TGeoNod *node = (TGeoNod *)  MG.NodeMap().contents(i);
		out<<(j+1)<<' '<<node->Coord(0)<<' '<<node->Coord(1)<<" 0"<<endl;
      MG.NodeMap().next(i);
      j++;
	}
	out<<"nomax "<< 4 <<"\n";
   int numel = MG.NumElem();
   out<<"elem "<<numel<<endl;
   j=0;
	i=MG.ElementMap().first();
   while (i) {
    	TGeoEl *element = (TGeoEl *)  MG.ElementMap().contents(i);
      int numsides = element->NumSides();
      out<<(j+1)<<' '<<  numsides <<' '<<(element->NodePtr(0)->Id()+1)<<' '<<(element->NodePtr(1)->Id()+1);
		out<<' '<<(element->NodePtr(2)->Id()+1);
      if ( numsides == 4) out<<' '<<(element->NodePtr(3)->Id()+1);
      out<<endl;
      MG.ElementMap().next(i);
      j++;
	}
}         */

/*   TPZStack<int> elgraph(0);
   TPZStack<int> elgraphindex(0);
   int nnod;
   secondmesh->ComputeElGraph(elgraph,elgraphindex,nnod);
   int nel = elgraphindex.NElements()-1;
   TPZMetis metis(nel,nnod);
   metis.SetElementGraph(elgraph,elgraphindex);
   TPZVec<int> perm(0), iperm(0);
   metis.Resequence(perm,iperm);
   secondmesh->Permute(perm);
*/

/*   TPZVec<char *> scalnames(3),vecnames(1);
   scalnames[0] = "SigmaX";
   scalnames[1] = "SigmaY";
   scalnames[2] = "TauXY";
   vecnames[0]  = "Displacement";
   //vecnames[1]  = "PrincipalStresses";
   an.DefineGraphMesh(2,scalnames,vecnames,"MvwModel.pos");
*/

/*   int bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel1;//elemento geometrico
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;//BC : AutoBuidl gerara o elemento computacional

   bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel2;
   bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;
*/

   //1o exemplo
//   elg2->Divide(vecsub);//2:4567
//   elg1->Divide(vecsub1);//1:891011
//   vecsub[2]->Divide(vecsub);//6:12131415
//	vecsub[1]->Divide(vecsub);//13:16171819
//   elg0->Divide(vecsub);//0:20212223
//   vecsub[3]->Divide(vecsub);//23:24252627
	//2o exemplo
//   elg2->Divide(vecsub);//2:4567
//   elg1->Divide(vecsub1);//1:891011
//   vecsub[0]->Divide(vecsub);//4:12131415
//	vecsub[1]->Divide(vecsub);//13:16171819
//   elg0->Divide(vecsub);//0:20212223
//   vecsub[3]->Divide(vecsub);//23:24252627
