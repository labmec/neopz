/**
 * @file
 * @brief Contains example to iterative analysis
 */
#include "pzconslaw.h"
#include "TPZConsLawTest.h"
#include "pzeulerconslaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
#include "TPZGMSHReadMesh.h"
#include "TPZIterativeAnalysis.h"
#include "TPZFlowCMesh.h"

#include "pzeltype.h"
#include "TPZGeoElement.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"

#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzskylmat.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
#include "pzreal.h"

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
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzdxmesh.h"

#include "pzflowcmesh.h"
#include "pzcmesh.h"


#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;
using namespace pzgeom;
using namespace pzrefine;

void LeituraDaMalha(char *meshfile,TPZStack<TPZGeoEl *> &elem,TPZStack<TPZGeoElSide> &elembc);
void LeituraDaMalha2(char *meshfile,TPZStack<TPZGeoEl *> &elem,TPZStack<TPZGeoElSide> &elembc);
void SetDeltaTime(TPZMaterial *mat,int nstate);
void CriacaoDeNos(int nnodes,double lista[20][3]);
TPZMaterial *Wing2d(int grau,TPZStack<TPZGeoElSide> &elembc);
TPZMaterial *Wing3d(int grau,TPZStack<TPZGeoElSide> &elembc);
void ContagemDeElementos(TPZMaterial *mat);
void Function(TPZVec<REAL> &x,TPZVec<REAL> &result);
//void PostProcess(TPZGeoMesh &gmesh,ostream &out);
void NivelDivide(TPZCompMesh *cmesh);
void Divisao(TPZCompMesh *cmesh);
void TestShapesDescontinous();
static clock_t start,end;//,begin,ttot=0;
void CoutTime(clock_t &start);
static TPZGeoMesh *gmesh = new TPZGeoMesh;
static TPZCompMesh *cmesh = new TPZFlowCompMesh(gmesh);
static TPZVec<REAL> x0(3,0.);
static int grau = 0;
//static int nivel = 0,tipo;
//static int problem=0;
//static REAL pi = 2.0*asin(1.0);
static REAL CFL=-1.0;
static REAL gama = 1.4;

/**
 *  Velocidade do som no ar a 15,9 graus centigrados = 340,9 mt/seg
 */

int main() {
	
	std::string name = "mesh.out";
	ofstream outgm(name.c_str());
	
	cout << "\nGrau do espaco de interpolacao -> 0,1,2,3,... ";
	cin >> grau;
	//  TPZCompElDisc::gDegree = grau;
	//  TPZCompEl::gOrder = grau;
	cmesh->SetDefaultOrder(grau);
	TPZAutoPointer<TPZMaterial> mat;
	TPZStack<TPZGeoEl *> elem;
	TPZStack<TPZGeoElSide> elembc;
	int dim;
	cout << "\nDimensao do problema 2,3\n"
	<< "\t[2D:Wing2d]\n"
	<< "\t[3D:Wing3d]\n"
	<< "\t\t\t";
	cin >> dim;
	TPZGMSHReadMesh readmesh(gmesh);
	if(dim == 2){
		//LeituraDaMalha2("naca0012_2d.msh",elem,elembc);
		readmesh.ReadMesh2D("naca0012_2d.msh",elem,elembc);
		mat = Wing2d(grau,elembc);
	}
	if(dim == 3){
		//LeituraDaMalha("naca0012_3d.msh",elem,elembc);
		readmesh.ReadMesh3D("naca0012_3d.msh",elem,elembc);
		mat = Wing2d(grau,elembc);
	}
	
	TPZVec<int> accumlist(10);
	accumlist[0] = 0;
	accumlist[1] = 1;
	accumlist[2] = 2;
	accumlist[3] = 2;
	accumlist[4] = 3;
	accumlist[5] = 0;
	accumlist[6] = 1;
	accumlist[7] = 2;
	accumlist[8] = 3;
	accumlist[9] = 0;
	
	//  TPZCompMesh *cmesh2 = cmesh->ComputeMesh(accumlist,2);
	// cmesh2->Print(outgm);
	return 0;
	
	if(1) {
		cout << "\ndescontinuo.c::main verificando a consistencia da malha de interfaces\n";
		if(TPZInterfaceElement::main(*cmesh)){
			cout << "->\t\t\tOK!";
		} else {
			cout << "->\tPROBLEMAS COM INTERFACES\n\n";
			//return ContagemDeElementos();
		}
	}
	
	if(1) {
		cout << "\nmain::Imprime malhas\n";
		gmesh->Print(outgm);
		cmesh->Print(outgm);
		outgm.flush();
	}
	
	int numiter,marcha;
	cout << "\nNumero de iteracoes requerida ? : ";
	cin >> numiter;
	//numiter = 1000;
	cout << "main:: Parametro marcha : \n";
	cin >> marcha;
	//marcha = 99;
	if(1){
		cout << "main::SetDeltaTime entre CFL (si nulo sera calculado) -> ";
		cin >> CFL;
		//CFL = 0.2;
		TPZDiffusionConsLaw::fCFL = CFL;
	}
	
	if(1){
		start = clock();
		int qual;
		cout << "main::[Refinamento nivel:1][Refinamento manual:2]\n";
		cin >> qual;
		if(qual == 2) Divisao(cmesh);
		if(qual == 1) NivelDivide(cmesh);
		CoutTime(start);
		int nstate = 4;
		SetDeltaTime(mat.operator->(),nstate);
		if(0){
			gmesh->Print(outgm);
			cmesh->Print(outgm);
			outgm.flush();
		}
	}
	
	if(1){
		cout << "\nmain::Ajuste no contorno e imprime malhas\n";
		cmesh->AdjustBoundaryElements();
		if(0){
			gmesh->Print(outgm);
			cmesh->Print(outgm);
			outgm.flush();
		}
	}
	
	if(1){
		TPZIterativeAnalysis an(cmesh,outgm);
		if(1){//Analysis
			cout << "\nmain::Resolve o sistema\n";
			//TPZStructMatrix *stiff;
			TPZSkylineStructMatrix stiff(cmesh);
			//TPZSpStructMatrix stiff(cmesh);
			an.SetStructuralMatrix(stiff);
			an.Solution().Zero();
			TPZStepSolver solver;
			solver.SetDirect(ELDLt);//ELU, ECholesky
			an.SetSolver(solver);
			if(1){
				REAL tol;
				tol = 1.0e15;// = norma da solu��o inicial + epsilon
				cout << "\nTolerancia ? : " << tol << "\n";
				//cin >> tol;
				//an.SetExact(Solution);
				int resolution=0;
				cout << "main:: Parametro resolution : \n";
				//cin >> resolution;
				resolution = 0; cout << resolution << "\n";
				an.IterativeProcess(name,tol,numiter,mat,marcha,resolution);
				//if(0) PostProcess(*gmesh,outgm);
			}
		}
		ContagemDeElementos(mat.operator->());
	}//if(0/1)
	
	if(0){
		gmesh->Print(outgm);
		cmesh->Print(outgm);
		outgm.flush();
	}
	
	outgm.close();
	if(cmesh) delete cmesh;
	if(gmesh) delete gmesh;
	//AvisoAudioVisual();
	return 0;
}

void LeituraDaMalha(char *meshfile,TPZStack<TPZGeoEl *> &elem,TPZStack<TPZGeoElSide> &elembc){
	
	ifstream mesh(meshfile);
	char title[256];
	int nnodes,number;
	mesh >> title;//$NOD
	mesh >> nnodes;
	gmesh->NodeVec().Resize(nnodes);
	TPZVec<REAL> coord(3);
	int i;
	for(i=0;i<nnodes;i++){
		mesh >> number;
		mesh >> coord[0] >> coord[1] >> coord[2];
		gmesh->NodeVec()[i].Initialize(coord,*gmesh);
	}
	mesh >> title;//$ENDNOD
	mesh >> title;//$ELM
	int numelem;
	mesh >> numelem;
	TPZVec<int> nodes;
	nodes.Resize(3);
	int index,numb,fgt1,nmat,fgt2,nvert;
	TPZVec<int> nos(8);//m�ximo do cubo
	for(i=0;i<numelem;i++){
		mesh >> numb >> fgt1 >> nmat >> fgt2 >> nvert;
		for(i=0;i<nvert;i++) mesh >> nos[i];
		if(nvert == 3){//tri�ngulos
			for(i=0;i<nvert;i++) nodes[i] = nos[i]-1;
			nodes.Resize(3);
			elembc.Push(TPZGeoElSide(gmesh->CreateGeoElement(ETriangle,nodes,1,index),-nmat));
			continue;
		} else if(nvert == 4 && fgt1 != 4){//quadrilateros
			nodes.Resize(4);
			for(i=0;i<nvert;i++) nodes[i] = nos[i]-1;
			elembc.Push(TPZGeoElSide(gmesh->CreateGeoElement(EQuadrilateral,nodes,1,index),-nmat));
			continue;
		} else if(nvert == 4 && fgt1 == 4){//tetraedros
			nodes.Resize(4);
			for(i=0;i<nvert;i++) nodes[i] = nos[i]-1;
			elem.Push(gmesh->CreateGeoElement(ETetraedro,nodes,1,index));
			continue;
		} else if(nvert == 8){//hexahedros
			nodes.Resize(8);
			for(i=0;i<nvert;i++) nodes[i] = nos[i]-1;
			elem.Push(gmesh->CreateGeoElement(ECube,nodes,1,index));
			continue;
		}
	}
	mesh >> title;//$ENDELM
	mesh.close();
}

void LeituraDaMalha2(char *meshfile,TPZStack<TPZGeoEl *> &elem,TPZStack<TPZGeoElSide> &elembc){
	
	ifstream mesh(meshfile);
	char title[256];
	int nnodes,number;
	mesh >> title;//$NOD
	mesh >> nnodes;
	gmesh->NodeVec().Resize(nnodes);
	TPZVec<REAL> coord(3);
	int i;
	for(i=0;i<nnodes;i++){
		mesh >> number;
		mesh >> coord[0] >> coord[1] >> coord[2];
		gmesh->NodeVec()[i].Initialize(coord,*gmesh);
	}
	mesh >> title;//$ENDNOD
	mesh >> title;//$ELM
	int numelem;
	mesh >> numelem;
	TPZVec<int> nodes;
	nodes.Resize(3);
	int index,numb,fgt1,nmat,fgt2,nvert,no1,no2,no3;
	for(i=0;i<numelem;i++){
		mesh >> numb >> fgt1 >> nmat >> fgt2 >> nvert;
		mesh >> no1;
		mesh >> no2;
		nodes[0] = no1-1;
		nodes[1] = no2-1;
		if(nvert == 2){//s� elementos de contorno
			nodes.Resize(2);
			elembc.Push(TPZGeoElSide(gmesh->CreateGeoElement(EOned,nodes,1,index),-nmat));
			continue;
		}
		if(nvert == 3){//s� elementos de volume
			nodes.Resize(3);
			mesh >> no3;
			nodes[2] = no3-1;
			elem.Push(gmesh->CreateGeoElement(ETriangle,nodes,1,index));
		}
	}
	mesh >> title;//$ENDELM
	mesh.close();
}

void SetDeltaTime(TPZMaterial *mat,int nstate){
	
	TPZVec<REAL> x(3,0.0),sol;
	int i;
	x[0] = 0.5;
	x[1] = 0.5;
	if(nstate==5) x[2] = 0.5;
	Function(x,sol);
	REAL prod = 0.0,maxveloc;
	for(i=1;i<nstate-1;i++) prod += sol[i]*sol[i];//(u�+v�+w�)*ro�
	REAL dens2 = sol[0]*sol[0];
	maxveloc = sqrt(prod/dens2);//velocidade
	TPZEulerConsLaw *law = dynamic_cast<TPZEulerConsLaw *>(mat);
	REAL press = law->Pressure(sol);
	if(press < 0) cout << "main::SetDeltaTime pressao negativa, toma valor absoluto para calculo do som\n";
	REAL sound = sqrt(law->Gamma()*press/sol[0]);
	maxveloc += sound;
	//REAL deltax = cmesh->DeltaX();
	REAL deltax = cmesh->LesserEdgeOfMesh();
	//REAL deltax = cmesh->MaximumRadiusOfMesh();
	REAL deltaT = CFL*deltax/maxveloc;
	cout << "main::SetDeltaTime : " << deltaT << endl;
	law->SetDelta(deltaT);
	
}

void Divisao (TPZCompMesh *cmesh,int key){
	
	if(key < 0) return;
	TPZVec<int> csub(0);
	int n1=1;
	while(n1) {
		cout << "\nId do elemento geometrico a dividir ? : ";
		cin >> n1;
		if(n1 < 0) break;
		int nelc = cmesh->ElementVec().NElements();
		int el=0;
		TPZCompEl *cpel=0;
		for(el=0;el<nelc;el++) {
			cpel = cmesh->ElementVec()[el];
			if(cpel && cpel->Reference()->Id() == n1) break;
		}
		if(cpel && el < nelc && cpel->Type() == 16){
			PZError << "main::Divisao elemento interface (nao foi dividido!)\n\n";
			cout << "Elementos divissiveis:\n";
			for(el=0;el<nelc;el++) {
				cpel = cmesh->ElementVec()[el];
				if(cpel && cpel->Type() != 16){
					TPZGeoEl *gel = cpel->Reference();
					if(gel) cout << gel->Id() << ",";
				}
			}
		} else {
			if(!el || el < nelc) cmesh->Divide(el,csub,0);
			else {
				cout << "main::Divisao elemento sem referencia\n";
				ContagemDeElementos(0);
			}
			n1 = 1;
		}
	}
}

void ContagemDeElementos(TPZMaterial *mat){
	
	int poin=0,line=0,tria=0,quad=0,tetr=0,pira=0,pris=0,hexa=0,disc=0,inte=0;
	int nelem = cmesh->ElementVec().NElements();
	int k,totel=0,bcel=0,niv = 0,nivmax=0;
	for(k=0;k<nelem;k++){
		TPZCompEl *comp = cmesh->ElementVec()[k];
		if(!comp) continue;
		totel++;
		if(comp->Reference()->MaterialId() < 0) bcel++;
		niv = comp->Reference()->Level();
		if(nivmax < niv) nivmax = niv;
		if(comp->Type() == 15) disc++;
		if(comp->Type() == 16) inte++;
	}
	nelem = gmesh->ElementVec().NElements();
	int total=0,nivmax2=0;
	for(k=0;k<nelem;k++){
		TPZGeoEl *geo = gmesh->ElementVec()[k];
		if(!geo) continue;
		total++;
		niv = geo->Level();
		if(nivmax2 < niv) nivmax2 = niv;
		if(geo->Reference()){
			int nsides = geo->NSides();
			if(nsides ==  1) poin++;
			if(nsides ==  3) line++;
			if(nsides ==  7) tria++;
			if(nsides ==  9) quad++;
			if(nsides == 15) tetr++;
			if(nsides == 19) pira++;
			if(nsides == 21) pris++;
			if(nsides == 27) hexa++;
		}
	}
	cout << "\nTotal de elementos computacionais  : " << totel;
	cout << "\nTotal de elementos de dominio      : " << abs(totel-bcel);
	cout << "\nTotal de elementos de contorno     : " << bcel;
	cout << "\nTotal de elementos ponto           : " << poin;
	cout << "\nTotal de elementos linha           : " << line;
	cout << "\nTotal de elementos triangulo       : " << tria;
	cout << "\nTotal de elementos quadrilatero    : " << quad;
	cout << "\nTotal de elementos tetraedro       : " << tetr;
	cout << "\nTotal de elementos piramide        : " << pira;
	cout << "\nTotal de elementos prisma          : " << pris;
	cout << "\nTotal de elementos hexaedro        : " << hexa;
	cout << "\nTotal do tipo discontinuo          : " << disc;
	cout << "\nTotal de tipo interface            : " << inte;
	cout << "\nTotal de nos                       : " << gmesh->NodeVec().NElements();
	cout << "\nTotal de elementos geometricos     : " << total;
	cout << "\nTamanho do vetor de connects       : " << cmesh->NConnects();
	cout << "\nTamanho do vetor el. comput.       : " << cmesh->NElements();
	cout << "\nGrau do espa�o de interpola��o     : " << grau;
	cout << "\nNivel maximo comput. atingido      : " << nivmax;
	cout << "\nNivel maximo geomet. atingido      : " << nivmax << endl << endl;
	if(mat){
		cout << "\nPropriedades materiais             : ";
		mat->Print();
	}
	cout << "\nDeltaX                             : " <<  cmesh->DeltaX() << endl;
	cout << "\nLesserEdgeOfEl                     : " <<  cmesh->LesserEdgeOfMesh() << endl;
	cout << "\nMaximumRadiusOfEl                  : " <<  cmesh->MaximumRadiusOfMesh() << endl;
}

int Nivel(TPZGeoEl *gel);
void NivelDivide(TPZCompMesh *cmesh){
	
	TPZVec<int> csub(0);
	int nivel;
	cout << "\nmain::Divisao todos os elementos da malha serao divididos!\n";
	cout << "\nmain::Divisao Nivel da malha final ? : ";
	cin >> nivel;
	//cout << "\nNivel da malha a ser atingido = " << nivel << endl;
	int nelc = cmesh->ElementVec().NElements();
	int el,actual;
	TPZCompEl *cpel;
	TPZGeoEl *gel;
	el = -1;
	while(++el<nelc) {
		cpel = cmesh->ElementVec()[el];
		if(!cpel) continue;
		if(cpel->Type() == 16) continue;
		if(cpel->Material()->Id() < 0) continue;
		gel = cpel->Reference();
		actual = Nivel(gel);
		if(actual < nivel){
			cmesh->Divide(el,csub,0);
			nelc = cmesh->ElementVec().NElements();
			el = -1;
			continue;
		}
	}
}

int Nivel(TPZGeoEl *gel){
	//retorna o n�vel do elemento gel
	if(!gel) return -1;
	TPZGeoEl *fat = gel->Father();
	if(!fat) return 0;
	int niv = 0;
	while(fat){
		fat = fat->Father();
		niv++;
	}
	return niv;
}

void Divisao(TPZCompMesh *cmesh){
	
	//int k=0;
	TPZVec<int> csub(0);
	int n1=1;
	while(n1) {
		cout << "Id do elemento geometrico a dividir ? : ";
		cin >> n1;
		cout << n1 << endl;
		if(n1 < 0) break;
		int nelc = cmesh->ElementVec().NElements();
		int el;
		TPZCompEl *cpel;
		for(el=0;el<nelc;el++) {
			cpel = cmesh->ElementVec()[el];
			if(!cpel) continue;
			if(cpel->Type() == 16) continue;
			if(cpel->Material()->Id() < 0) continue;
			if(cpel && cpel->Reference()->Id() == n1) break;
		}
		if(el < nelc) cmesh->Divide(el,csub,0);
		else cout << "Divisao::Elemento n�o existe\n";
		n1 = 1;
	}
}

void CoutTime(clock_t &start){
    end = clock();
    cout << "\nFim da etapa : "  <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *Wing2d(int grau,TPZStack<TPZGeoElSide> &elembc){
	
	//elemento de volume
	// TPZGeoElement</*TPZShapeTriang,*/TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	// TPZGeoElement</*TPZShapeLinear,*/TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//int interfdim = 1;
	int i;
	// TPZCompElDisc::gInterfaceDimension = interfdim;
	gmesh->BuildConnectivity();
	int nummat = 1;
	int nivel;
	TPZArtDiffType artdiff = LeastSquares_AD;
	cout << "\nmain::Divisao Nivel final da malha ? : ";
	cin >> nivel;
	REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );///0.5;
	REAL delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );//0.5;
	REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
	//calculando novos valores
	delta_t = delta_x*cfl;
	REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
	gama = 1.4;
	cout << "\nDominio [0,1]x[0,1]"
	<< "\nMax df/dx (desconhecido) = 1.0"
	<< "\nCFL = " << cfl
	<< "\ndelta otimo = " << delta
	<< "\nDelta x = " << delta_x
	<< "\ndelta t = " << delta_t
	<< "\ndiffusao = " << artdiff
	<< "\ndelta aproximado = " << delta << endl;
	
	int dim = 2;
	TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
	DebugStop(); //mat->SetForcingFunction(Function);
	cmesh->InsertMaterialObject(mat);
	
	// boundary conditions
	TPZBndCond *bc;
	TPZFMatrix val1(4,4),val2(4,1);
	
	//CC : a vizinhan�a geometrica foi preenchida
	val1.Zero();
	val2.Zero();
	int numbc = elembc.NElements();
	cout << "main::Wing2d criando CC\n";
	for(i=0;i<numbc;i++){
		TPZGeoEl *elgbound = elembc[i].Element();//elemento de contorno
		if(!elgbound){
			cout << "main::Wing2d erro, elemento geometrico nulo\n";
			continue;
		}
		int sideinner = elgbound->NSides()-1;
		TPZGeoElSide elgvolside = elgbound->Neighbour(sideinner);//elemento de volume
		TPZGeoEl *elgvol = elgvolside.Element();
		int volside = elgvolside.Side();
		if(!elgvol){
			cout << "main::Wing2d erro, vizinho geometrico de volume nulo\n";
			continue;
		}
		int typecc = elembc[i].Side();
		if(typecc == -1){
			TPZGeoElBC(elgvol,volside,typecc);//CC da asa: wall
		} else {
			TPZGeoElBC(elgvol,volside,typecc);//CC far field: nonreflecting
		}
		elgbound->RemoveConnectivities();
		int index = gmesh->ElementIndex(elgbound);// identifica o index do elemento
		gmesh->ElementVec()[index] = NULL;
		delete elgbound;
		gmesh->ElementVec().SetFree(index);
	}
	TPZAutoPointer<TPZMaterial> aximat(mat);
	bc = mat->CreateBC(aximat,-1,5,val1,val2);//parede
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-2,2,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-3,2,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-4,2,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-5,2,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	cout << "main::Wing2d fim CC\n";
	
	cout << endl;
	cmesh->AutoBuild();
	
	return mat;
}
////////////////////////////////////
void Function(TPZVec<REAL> &x,TPZVec<REAL> &result){
	
    result.Resize(4);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 0.5;
    REAL v = 0.0;
    REAL p = 1.0;
    REAL vel2 = u*u + v*v;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = p/(gama-1.0) + 0.5 * ro * vel2;
    return;
	
}

//----------------------------------------------------------------------------------------------
TPZMaterial *Wing3d(int grau,TPZStack<TPZGeoElSide> &elembc){
	
	//elemento de volume
	//TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
	//int interfdim = 2;
	int i;
	int nivel;
	// TPZCompElDisc::gInterfaceDimension = interfdim;
	gmesh->BuildConnectivity();
	int nummat = 1;
	TPZArtDiffType artdiff = LeastSquares_AD;
	cout << "\nmain::Divisao Nivel final da malha ? : ";
	cin >> nivel;
	REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );///0.5;
	REAL delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );//0.5;
	REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
	REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
	gama = 1.4;//do ar
	cout << "\nDominio [0,1]x[0,1]"
	<< "\nMax df/dx (desconhecido) = 1.0"
	<< "\nCFL = " << cfl
	<< "\ndelta otimo = " << delta
	<< "\nDelta x = " << delta_x
	<< "\ndelta t = " << delta_t
	<< "\ndiffusao = " << artdiff
	<< "\ndelta aproximado = " << delta << endl;
	
	int dim = 3;
	TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
	DebugStop(); //mat->SetForcingFunction(Function);
	cmesh->InsertMaterialObject(mat);
	
	// boundary conditions
	TPZBndCond *bc;
	TPZFMatrix val1(5,4),val2(5,1);
	
	//CC : a geometric neighboard was filled
	val1.Zero();
	val2.Zero();
	int numbc = elembc.NElements();
	cout << "main::Wing3d criando CC\n";
	for(i=0;i<numbc;i++){
		TPZGeoEl *elgbound = elembc[i].Element();//elemento de contorno
		if(!elgbound){
			cout << "main::Wing3d erro, elemento geometrico nulo\n";
			continue;
		}
		int sideinner = elgbound->NSides()-1;
		TPZGeoElSide elgvolside = elgbound->Neighbour(sideinner);//elemento de volume
		TPZGeoEl *elgvol = elgvolside.Element();
		int volside = elgvolside.Side();
		if(!elgvol){
			cout << "main::Wing3d erro, vizinho geometrico de volume nulo\n";
			continue;
		}
		int typecc = elembc[i].Side();
		if(typecc == -1){
			TPZGeoElBC(elgvol,volside,typecc);//CC da asa: wall
		} else {
			TPZGeoElBC(elgvol,volside,typecc);//CC far field: nonreflecting
		}
		elgbound->RemoveConnectivities();
		int index = gmesh->ElementIndex(elgbound);// identifica o index do elemento
		gmesh->ElementVec()[index] = NULL;
		delete elgbound;
		gmesh->ElementVec().SetFree(index);
	}
	//the domain is a hexahedral 
	TPZAutoPointer<TPZMaterial> aximat(mat);
	bc = mat->CreateBC(aximat,-1,5,val1,val2);//parede na asa
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-2,6,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-3,6,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-4,6,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-5,6,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	bc = mat->CreateBC(aximat,-6,6,val1,val2);//no refletivas
	cmesh->InsertMaterialObject(bc);
	cout << "main::Wing3d fim CC\n";
	cout << "main::Wing3D criando elementos computacionais\n";
	cmesh->AutoBuild();
	cout << "main::Wing3D elementos computacionais gerados\n";
	return mat;
}



