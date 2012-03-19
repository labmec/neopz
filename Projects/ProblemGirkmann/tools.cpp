/*
 *  tools.cpp
 *  
 *  Created by Agnaldo on 4/29/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
using namespace std;

#include "tools.h"

#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"

#include <iostream>
#include <string>
#include <math.h>

#include "pzelasAXImat.h" 
#include "pzfstrmatrix.h"
#include "pzbstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZInterfaceEl.h"


#include "pzlog.h"
#ifdef LOG4CXX
//static LoggerPtr logger(Logger::getLogger("agnaldo.integrate"));
static LoggerPtr loggertensoes(Logger::getLogger("agnaldo.tensoes"));
#endif


const int mat1Id	=  1; // Material casca (recebe Peso proprio)
const int mat2Id	=  2; // Material anel (recebe Peso proprio)
const int mat3Id	=  8; // Material interface

const int mat1ArcGeoUm = -1; //Contorno inferior da casca
const int mat1ArcGeoDois = -2; //Contorno superior da casca
const int mat1ArcBC	= -3; //superficie media da casca (recebe Carregamento distribuido mais peso proprio)
const int mat1EngBC	= -4; //Restricao de deslocamento em x e giro
const int mat1JuncaoPoint = -5; // Ponto na juncao casca e anel, no lado do anel

const int mat2BaseBC = -6; //Contorno da Base do anel/Reacao de apoio
const int mat2BaseBCPoint = -7; //Contorno da Base do anel/Aplicar mola


const int matRefDir1 = 3; //Juncao entre a casca e o anel
const int matPointRefDir1 = 4; //Ponto na Casca: Juncao casca e anel
const int matPointRefDir2 = 5; //Ponto na Casca: Juncao casca e anel
const int mat1EngPoint1 = 6; 
const int mat1EngPoint2 = 7; 

const REAL Pi = atan(1.)*4.;

const int dirichlet = 0;
const int neumann = 1;
const int mista = 2;
const int hidrost = 3;

///----------------------- Construtor -------------------------------
REAL tools::fRc; REAL tools::fh; REAL tools::falpha;
REAL tools::fa; REAL tools::fb; REAL tools::fRho;
REAL tools::fyoung; REAL tools::fpoisson; bool tools::fAnelComPesoLeve;

tools::tools(REAL Rc, REAL h, REAL alpha, REAL a, REAL b, REAL Rho,  REAL young, REAL poisson, bool AnelComPesoLeve)
{
	fRc=0.; fh=0.; falpha=0.;
	fa=0.; fb=0.; fRho=0.;
	fyoung=0.; fpoisson=0.;
	fAnelComPesoLeve = true;
	
	fRc=Rc;
	fh=h;
	falpha=alpha;
	fa=a;
	fb=b;
	fRho=Rho;
	//fTz=Tz;
	fyoung=young;
	fpoisson=poisson;
	fAnelComPesoLeve=AnelComPesoLeve;
}


tools::~tools()
{	
}

///---------------------- MatrixR() ---------------------------------------------------
TPZFMatrix<REAL> tools::MatrixR(REAL theta)
{
	TPZFMatrix<REAL> r(3,3,0.);
	r(0,0) = cos(theta); r(0,1) = sin(theta);
	r(1,0) = -sin(theta); r(1,1) =  cos(theta);
	
	return r;
}

///---------------------- ShellGen() -----------------------------------------------
TPZGeoMesh * tools::MalhaGeoGen(int ndiv, int ndirectdivR,int ndirectdivL,int ndirectdivp, bool interface, int RefDirId)
{
	//ndiv = 2;
	int Qnodes = 6*(ndiv + 1);
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	
	gmesh->NodeVec().Resize(Qnodes);
	
	TPZFMatrix<REAL> v(3,1,0.), vl(3,1,0.), vltemp(3,1,0.);
	v(0,0) = 0.;
	v(1,0) = 1.;
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolArc(3);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	REAL theta = falpha/(ndiv*2);
	REAL Rm = fRc/sin(falpha);//23.3359;(raio da superficie media da casca esferica)
	
	//Indice dos nos
	int id = 0;
	double val1, val2;
	//int i;
	for(int cn = 0; cn <= 2*ndiv; cn++)
	{
		MatrixR(cn*theta).Multiply(v,vl);
		
		vltemp = vl;
		val1 = vltemp(0,0) * (Rm - fh/2.);
		val2 = vltemp(1,0) * (Rm - fh/2.);
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 , val1);//coord X
		Node[id].SetCoord(1 , val2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
		
		vltemp = vl;
		val1 = vltemp(0,0)*Rm;
		val2 = vltemp(1,0)*Rm;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 , val1);//coord X
		Node[id].SetCoord(1 , val2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
		
		vltemp = vl;
		val1 = vltemp(0,0) * (Rm +fh/2.);
		val2 = vltemp(1,0) * (Rm + fh/2.);
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 , val1);//coord X
		Node[id].SetCoord(1 , val2);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//(0,0)
	Node[id].SetNodeId(id);
	Node[id].SetCoord(0 , fRc - fh*sin(falpha)/2.);//coord X
	Node[id].SetCoord(1 , (Rm+fh/2.)*cos(falpha)-fb);//coord Y
	gmesh->NodeVec()[id] = Node[id];
	id++;
	
	//(a,0)
	Node[id].SetNodeId(id);
	Node[id].SetCoord(0 , fRc + fa - fh*sin(falpha)/2.);//coord X
	Node[id].SetCoord(1 , (Rm+fh/2.)*cos(falpha)-fb);//coord Y
	gmesh->NodeVec()[id] = Node[id];
	id++;
	
	//(a,b)
	Node[id].SetNodeId(id);
	Node[id].SetCoord(0 , fRc + fa - fh*sin(falpha)/2.);//coord X
	Node[id].SetCoord(1 , (Rm+fh/2.)*cos(falpha));//coord Y
	gmesh->NodeVec()[id] = Node[id];
	
	//Indice dos Elementos
	//Casca
	id = 0;
//	static int meshcount = 0;
	
	for(int ce = 0; ce < ndiv; ce++)
	{
		TopolQuad[0] = 6*ce;
		TopolQuad[1] = TopolQuad[0] + 6;
		TopolQuad[2] = TopolQuad[1] + 1;
		TopolQuad[3] = TopolQuad[0] + 1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (id,TopolQuad,mat1Id,*gmesh);
		id++;
		
		TopolQuad[0] += 1;	
		TopolQuad[1] += 1;
		TopolQuad[2] += 1;
		TopolQuad[3] += 1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (id,TopolQuad,mat1Id,*gmesh);
		id++;
		
		TopolArc[0] = 6*ce;
		TopolArc[1] = TopolArc[0] + 6;
		TopolArc[2] = TopolArc[0] + 3;
		new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,mat1ArcGeoUm,*gmesh);
		id++;
		
		TopolArc[0] += 1;
		TopolArc[1] += 1;
		TopolArc[2] += 1;
		new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,mat1ArcBC,*gmesh);
		id++;
			
		TopolArc[0] += 1;
		TopolArc[1] += 1;
		TopolArc[2] += 1;
		new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,mat1ArcGeoDois,*gmesh);
		id++;			 
	}
		
		//---------------------
		TopolLine[0] = 0;
		TopolLine[1] = 1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1EngBC,*gmesh);
		id++;
		
		TopolLine[0] = 1;
		TopolLine[1] = 2;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat1EngBC,*gmesh);
		id++;
		
		TopolPoint[0] = 0;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,mat1EngPoint1, *gmesh);
		id++;
		
		TopolPoint[0] = 2;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,mat1EngPoint2, *gmesh);
		id++;
		
		//Indice Elementos do anel
		TopolQuad[0] = Qnodes - 3;
		TopolQuad[1] = Qnodes - 2;
		TopolQuad[2] = Qnodes - 5;
		TopolQuad[3] = Qnodes - 6;
		new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,mat2Id,*gmesh);
		id++;
		
		TopolQuad[0] = Qnodes - 1;
		TopolQuad[1] = Qnodes - 4;
		TopolQuad[2] = Qnodes - 5;
		TopolQuad[3] = Qnodes - 2;
		new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,mat2Id,*gmesh);
		id++;
		
		TopolLine[0] = Qnodes - 3;
		TopolLine[1] = Qnodes - 2;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,mat2BaseBC, *gmesh);
		id++;
		
		TopolPoint[0] = Qnodes-3;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,mat2BaseBCPoint, *gmesh);
		id++;
		
		//--------------------------- Refinamento Direcional ---------------------------------
		TopolPoint[0] = Qnodes - 6;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,matPointRefDir1,*gmesh);
		id++;
		
		TopolPoint[0] = Qnodes - 4;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,matPointRefDir2,*gmesh);
		id++;
		
		TopolLine[0] = Qnodes - 6;
		TopolLine[1] = Qnodes - 5;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matRefDir1,*gmesh);
		id++;
		
		TopolLine[0] = Qnodes - 5;
		TopolLine[1] = Qnodes - 4;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matRefDir1,*gmesh);
		
		//---------  Elemento de interface  --------------------------------------------
		if (interface==true) 
		{
			gmesh->AddInterfaceMaterial(mat2Id,mat1Id, mat3Id); //Adicionar um material de interface associados aos elementos mat2 e mat1 do material.
			gmesh->AddInterfaceMaterial(mat1Id,mat2Id, mat3Id);
		}
		
		gmesh->BuildConnectivity();
		
		//------------------ fazer refinamento direcional ------------------------------------
		set<int> SETmatRefDir11;
		for(int j = 0; j < ndirectdivR; j++)
		{
			int nel = gmesh->NElements();
			for (int iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP11 = gmesh->ElementVec()[iref];
				if(!gelP11) continue;
				SETmatRefDir11.insert(matRefDir1);
				if (RefDirId != 0) {
					int matid= gelP11->MaterialId();
					if(matid==1) TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11);
				}
				else TPZRefPatternTools::RefineDirectional(gelP11, SETmatRefDir11);
			}		
		}
		
		set<int> SETmatRefDir12;
		for(int j = 0; j < ndirectdivL; j++)
		{
			int nel = gmesh->NElements();
			for (int iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP12 = gmesh->ElementVec()[iref];
				if(!gelP12) continue;
				SETmatRefDir12.insert(matRefDir1);
				if (RefDirId != 0) {
					int matid= gelP12->MaterialId();
					if(matid==2) TPZRefPatternTools::RefineDirectional(gelP12, SETmatRefDir12);
				}
				else TPZRefPatternTools::RefineDirectional(gelP12, SETmatRefDir12);
			}		
		}
		
		set<int> SETmatPointRefDir2;
		for(int j = 0; j < ndirectdivp; j++)
		{
			int nel = gmesh->NElements();
			for (int iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP2 = gmesh->ElementVec()[iref];
				if(!gelP2 || gelP2->HasSubElement()) continue;
				SETmatPointRefDir2.insert(matPointRefDir1);
				TPZRefPatternTools::RefineDirectional(gelP2, SETmatPointRefDir2);
			}		
		}
		
		set<int> SETmatPointRefDir3;
		for(int j = 0; j < ndirectdivp; j++)
		{
			int nel = gmesh->NElements();
			for (int iref = 0; iref < nel; iref++)
			{
				TPZVec<TPZGeoEl*> filhos;
				TPZGeoEl * gelP3 = gmesh->ElementVec()[iref];
				if(!gelP3) continue;
				SETmatPointRefDir3.insert(matPointRefDir2);
				TPZRefPatternTools::RefineDirectional(gelP3, SETmatPointRefDir3);
			}		
		}
	
	return gmesh;
}

//--------------------------- RefinamentoUniforme() ------------------------------------
void tools::RefinamentoUniforme(TPZGeoMesh & gMesh, int &nh)
{
	//int h =nh;
	for ( int ref = 0; ref < nh; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh.NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gMesh.ElementVec() [i];
			//int ind =  gel->MaterialId();
			if ( gel->Dimension() == 2  || gel->MaterialId()== matRefDir1 ) gel->Divide ( filhos );
		}//for i
	}//ref
}

///--------------------------------------- AreaBaseAnel() ------------------------------- 
REAL tools::AreaBaseAnel()
{ 
	REAL R1 = fRc - fh*sin(falpha)/2. ;
	REAL R2 = fRc + (fa - fh*sin(falpha)/2. );
	REAL area = Pi*(R2*R2 - R1*R1);
	return area;
}

///---------------------------------------  VMaterial() -----------------------------------
REAL tools::VMaterial()
{
	REAL Rm2 = fRc*tan(falpha/2.);
	REAL Rm = fRc/sin(falpha); 
	REAL Vcalota = 2.*Pi*fh*Rm*Rm2 ;
	
	REAL aa = fRc - fh*sin(falpha)/2. ;
	REAL bb = fRc + fa - fh*sin(falpha)/2.;
	REAL cc = fRc + fh*sin(falpha)/2.;
	
	REAL CG1 = aa + 0.5*(cc - aa);
	REAL CG2 = bb - 0.5*(bb - cc);
	REAL CG3 = CG2;
	REAL CG4 = aa + 2.*(cc - aa)/3.;
	
	REAL A1 = (fh*sin(falpha))*(fb - fh*cos(falpha));
	REAL A2 = (fa - fh*sin(falpha))*(fb - fh*cos(falpha));
	REAL A3 = (fa - fh*sin(falpha))*(fh*cos(falpha));
	REAL A4 = 0.5*(fh*cos(falpha))*(fh*sin(falpha));
	
	REAL Vbase = 2.*Pi*(CG1*A1 + CG2*A2 + CG3*A3 + CG4*A4);
	
	REAL VTotal; 
	REAL VCalculado = 0.;
	bool temp = fAnelComPesoLeve;
	if (temp==false) {
		VTotal = Vcalota+Vbase;
		Vcalota = 48.029962432109642; //resultado calculado numericamente, mais preciso
		Vbase = 28.719973111966411; //resultado calculado numericamente
		VCalculado =  Vcalota+Vbase;	
	}else{
		VTotal = Vcalota;
		Vcalota = 48.029962432109642; //resultado calculado numericamente, mais preciso
		VCalculado = Vcalota;
	}
	
	return VCalculado;
}

///------------------------------------------- MalhaCompGen() ---------------------- 
TPZCompMesh * tools::MalhaCompGen(TPZGeoMesh * gmesh, int p)
{
	
	REAL pesoProprio = fRho;
	REAL Vmat = VMaterial();
	REAL Abase = AreaBaseAnel();
	REAL FDistribuidaBase = (fRho*Vmat)/Abase;
//	REAL BigNum =TPZMaterial::gBigNumber;
	
	REAL fr = 0.;
	REAL fz_casca, fz_anel;
	fz_casca = -pesoProprio;  //(sobrecarga[ KN/m2 ])
	fz_anel = fz_casca;
	
	bool temp = fAnelComPesoLeve;
	if (temp==true) fz_anel=0.;
	//REAL fz = 0.; //so no caso de cond contorno hidrostatica
	
	TPZAutoPointer<TPZMaterial> matCasca = new TPZElasticityAxiMaterial(mat1Id, fyoung, fpoisson, fr, fz_casca);
	TPZAutoPointer<TPZMaterial> matAnel = new TPZElasticityAxiMaterial(mat2Id, fyoung, fpoisson, fr, fz_anel);
	TPZAutoPointer<TPZMaterial> matInterface = new TPZElasticityAxiMaterial(mat3Id, fyoung, fpoisson, fr, fz_casca);
	
	TPZManVector<REAL> Orig(3);		Orig[0]  = 0.;		Orig[1]  = 0.;		Orig[2]  = 0.;
	TPZManVector<REAL> AxisR(3);		AxisR[0] = 1.;	AxisR[1] = 0.;	AxisR[2] = 0.;
	TPZManVector<REAL> AxisZ(3);		AxisZ[0] = 0.;	AxisZ[1] = 1.;	AxisZ[2] = 0.;
	
	TPZElasticityAxiMaterial * aximat1 = dynamic_cast<TPZElasticityAxiMaterial*>(matCasca.operator->());
	TPZElasticityAxiMaterial * aximat2 = dynamic_cast<TPZElasticityAxiMaterial*>(matAnel.operator->());
	TPZElasticityAxiMaterial * aximat3 = dynamic_cast<TPZElasticityAxiMaterial*>(matInterface.operator->());
	
	aximat1->SetOrigin(Orig, AxisZ, AxisR);
	aximat2->SetOrigin(Orig, AxisZ, AxisR);
	aximat3->SetOrigin(Orig, AxisZ, AxisR);
	
	
	///Computational Mesh
	TPZCompEl::SetgOrder(p);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetAllCreateFunctionsContinuous();
	
	cmesh->InsertMaterialObject(matCasca);
	cmesh->InsertMaterialObject(matAnel);
	cmesh->InsertMaterialObject(matInterface);
	
	{
		//Boundary Conditions
		
		//Inserir Mola em um ponto da base do anel 
		TPZFMatrix<REAL> k1(2,2,0.), k2(2,2,0.), f(2,1,0.);
		k1(1,1) = 1.e-3;
		TPZAutoPointer<TPZMaterial> ContBC = matAnel->CreateBC(matAnel, mat2BaseBCPoint, mista, k1, f);
		cmesh->InsertMaterialObject(ContBC);
		
		k2(0,0) = 1.;
		TPZAutoPointer<TPZMaterial> ContBC2 = matCasca->CreateBC(matCasca, mat1EngBC, mista, k2, f);
		//cmesh->InsertMaterialObject(ContBC2);
		
		//Reacao de apoio no anel
		TPZFMatrix<REAL> Reac1(2,2,0.), Reac2(2,1,0.);
		Reac2(1,0) = FDistribuidaBase; 
		TPZAutoPointer<TPZMaterial> Cont = matAnel->CreateBC(matAnel, mat2BaseBC, neumann, Reac1, Reac2);
		cmesh->InsertMaterialObject(Cont);
		
		//Calcular momento e cortante
		TPZFMatrix<REAL> Reac12(2,2,0.), Reac22(2,1,0.);
		TPZAutoPointer<TPZMaterial> ContBC1 = matAnel->CreateBC(matAnel, matRefDir1, neumann, Reac12, Reac22);
		cmesh->InsertMaterialObject(ContBC1);
		
		
		
		//------------------- cond contorno Hidrostatica: so para validar o programa---------------------
		//TPZFMatrix<REAL> val1(2,2,0.),val2(2,1,0.);
		//		REAL peso = 100.;
		//				
		//		val2(0,0) = -peso;
		//		TPZAutoPointer<TPZMaterial> CondBc1 = matCasca->CreateBC(matCasca, mat1ArcGeoDois, hidrost, val1,val2);
		//		cmesh->InsertMaterialObject(CondBc1);
		//		
		//		TPZAutoPointer<TPZMaterial> CondBc2 = matCasca->CreateBC(matCasca, mat1EngBC, hidrost, val1,val2);
		//		cmesh->InsertMaterialObject(CondBc2);
		//		
		//		val2(0,0) = peso; 
		//		TPZAutoPointer<TPZMaterial> CondBc3 = matCasca->CreateBC(matCasca, mat1ArcGeoUm, hidrost, val1,val2);
		//		cmesh->InsertMaterialObject(CondBc3);
		//		
		//		TPZAutoPointer<TPZMaterial> CondBc4 = matAnel->CreateBC(matAnel, mat2BaseBC, hidrost, val1,val2);
		//		cmesh->InsertMaterialObject(CondBc4);
		//		
		//		TPZAutoPointer<TPZMaterial> CondBc5 = matAnel->CreateBC(matAnel, mat2BaseBCDois, hidrost, val1,val2);
		//		cmesh->InsertMaterialObject(CondBc5);
		//		
		//		TPZFMatrix<REAL> k1(2,2,0.), k2(2,2,0.), f(2,1,0.);
		//		k1(1,1) = 1.;
		//		TPZAutoPointer<TPZMaterial> ContBC6 = matAnel->CreateBC(matAnel, mat2BaseBCPoint, mista, k1, f);
		//		cmesh->InsertMaterialObject(ContBC6);
		//		
		//		k2(0,0) = 1.;
		//		TPZAutoPointer<TPZMaterial> ContBC7 = matCasca->CreateBC(matCasca, mat1EngBC, mista, k2, f);
		//		cmesh->InsertMaterialObject(ContBC7);
		//		
		//		TPZFMatrix<REAL> Reac12(2,2,0.), Reac22(2,1,0.);
		//		TPZAutoPointer<TPZMaterial> ContBC = matAnel->CreateBC(matAnel, matRefDir1, neumann, Reac12, Reac22);
		//		cmesh->InsertMaterialObject(ContBC);
		//-------------------------------------------------------------------------------
	}
		
	cmesh->AutoBuild();
	
	//ofstream arg("cmesh.txt");
	//	cmesh->Print(arg);
	
	return cmesh;
}


///--------------------------------MalhaCompMeshWithInterface()--------------------------------------------
TPZCompMesh * tools::MalhaCompMeshWithInterface(TPZGeoMesh * gmesh, int p, REAL simetric, REAL penalidade)
{
	
	REAL pesoProprio = fRho;
	REAL Vmat = VMaterial();
	REAL Abase = AreaBaseAnel();
	REAL FDistribuidaBase = (fRho*Vmat)/Abase;
	
	REAL fr = 0.;
	REAL fz_casca, fz_anel;
	fz_casca = -pesoProprio;  //(sobrecarga[ KN/m2 ])
	fz_anel = fz_casca;
	
	bool temp = fAnelComPesoLeve;
	if (temp==true) fz_anel=0.;
	
	TPZAutoPointer<TPZMaterial> matCasca = new TPZElasticityAxiMaterial(mat1Id, fyoung, fpoisson, fr, fz_casca, simetric, penalidade);
	TPZAutoPointer<TPZMaterial> matInterface = new TPZElasticityAxiMaterial(mat3Id, fyoung, fpoisson, fr,  fz_casca, simetric, penalidade);
	TPZAutoPointer<TPZMaterial> matAnel = new TPZElasticityAxiMaterial(mat2Id, fyoung, fpoisson, fr, fz_anel, simetric, penalidade);
	
	TPZManVector<REAL> Orig(3);		Orig[0]  = 0.;		Orig[1]  = 0.;		Orig[2]  = 0.;
	TPZManVector<REAL> AxisR(3);		AxisR[0] = 1.;	AxisR[1] = 0.;	AxisR[2] = 0.;
	TPZManVector<REAL> AxisZ(3);		AxisZ[0] = 0.;	AxisZ[1] = 1.;	AxisZ[2] = 0.;
	
	TPZElasticityAxiMaterial * aximat1 = dynamic_cast<TPZElasticityAxiMaterial*>(matCasca.operator->());
	TPZElasticityAxiMaterial * aximat2 = dynamic_cast<TPZElasticityAxiMaterial*>(matAnel.operator->());
	TPZElasticityAxiMaterial * aximat3 = dynamic_cast<TPZElasticityAxiMaterial*>(matInterface.operator->());
	
	aximat1->SetOrigin(Orig, AxisZ, AxisR);
	aximat2->SetOrigin(Orig, AxisZ, AxisR);
	aximat3->SetOrigin(Orig, AxisZ, AxisR);
	
	///Computational Mesh
	TPZCompEl::SetgOrder(p);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->InsertMaterialObject(matCasca);
	cmesh->InsertMaterialObject(matAnel);
	cmesh->InsertMaterialObject(matInterface);
	
	{
		//Boundary Conditions
		
		//Inserir Mola em um ponto da base do anel 
		TPZFMatrix<REAL> k1(2,2,0.), k2(2,2,0.), f(2,1,0.);
		k1(1,1) = 1.e-3;
		
		TPZAutoPointer<TPZMaterial> BCPoint7 = matAnel->CreateBC(matAnel, mat2BaseBCPoint, mista, k1, f);
		cmesh->InsertMaterialObject(BCPoint7);
		
		k2(0,0) = 1.;
		TPZAutoPointer<TPZMaterial> ContBC3 = matCasca->CreateBC(matCasca, mat1EngPoint1, mista, k2, f);
		//cmesh->InsertMaterialObject(ContBC3);
		
		//Reacao de apoio no anel
		TPZFMatrix<REAL> Reac1(2,2,0.), Reac2(2,1,0.);
		Reac2(1,0) = FDistribuidaBase; 
		TPZAutoPointer<TPZMaterial> Cont = matAnel->CreateBC(matAnel, mat2BaseBC, neumann, Reac1, Reac2);
		cmesh->InsertMaterialObject(Cont);
		
		//Calcular momento e cortante
		TPZFMatrix<REAL> Reac12(2,2,0.), Reac22(2,1,0.);
		TPZAutoPointer<TPZMaterial> ContBC2 = matAnel->CreateBC(matAnel, matRefDir1, neumann, Reac12, Reac22);
		cmesh->InsertMaterialObject(ContBC2);
		
	}
	
	//AQUI: SetAllCreateFunctionsContinuous
    cmesh->SetAllCreateFunctionsContinuous();
	
	//AQUI: AutoBuild(mat1)
	set<int> SETmat1;
	SETmat1.insert(mat1Id);
	SETmat1.insert(mat1EngBC);
	cmesh->AutoBuild(SETmat1);
	gmesh->ResetReference();
	
	//AQUI: AutoBuild(mat2)
	set<int> SETmat2;
	SETmat2.insert(mat2Id);
	SETmat2.insert(mat2BaseBCPoint);
	SETmat2.insert(mat2BaseBC);
	SETmat2.insert(matRefDir1);
	cmesh->AutoBuild(SETmat2);
	cmesh->LoadReferences(); 
	
	//AQUI: Criar elemento de interface
	for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh->ElementVec()[el];
		if(!compEl) continue;
		int matId = compEl->Reference()->MaterialId();
		if((matId== matRefDir1)||(matId== matPointRefDir1)||(matId== matPointRefDir2))
		{
			compEl->Reference()->ResetReference();
		}
	}
	
	for(int el = 0; el < cmesh->ElementVec().NElements(); el++)
	{
		TPZCompEl * compEl = cmesh->ElementVec()[el];
		if(!compEl) continue;
		int index = compEl ->Index();
		if((compEl) && (compEl->Dimension() == 2) && (compEl->Reference()->MaterialId() ==mat2Id))
		{
			//cout<<"----------------------------------"<<endl;
			//compEl->Reference()->Print();
			//cout << "--------------------------------"<<endl;
			
			TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
			int nsides = compEl->Reference()->NSides();
			TPZStack<TPZCompElSide> neighbours;
			
			int isini;//inicial 
			if (nsides==7) isini=3;
			else isini=4;
			
			for (int is=isini; is<nsides-1; is++)
			{
				neighbours.Resize(0);
				TPZCompElSide celside(compEl,is);
				celside.EqualLevelElementList(neighbours, 1, 1);
				//celside.LowerLevelElementList(1);
				//celside.RemoveDuplicates(neighbours);
				int nneig = neighbours.NElements();
				if (nneig>1) DebugStop();
				if (nneig==1)
				{
					int idneig = neighbours[0].Reference().Element()->MaterialId();
					if (idneig==mat1Id)
					{
						//TPZInterfaceElement * FaceElem = 
						InterpEl->CreateInterface(is,true);
					}
				}
			} 
		}
	}
	return cmesh;
}


///-----------------------------------------SolveSist()---------------------------------------------
#define VTK
void tools::SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh, int sim)
{
	//TPZFStructMatrix full(fCmesh);
	//TPZFrontStructMatrix<TPZFrontNonSym> full(fCmesh);
	//full.SetNumberOfThreads(3);
	
	TPZStepSolver<REAL> step;
	if(sim==1){
		TPZBandStructMatrix full(fCmesh);
		an.SetStructuralMatrix(full);
		step.SetDirect(ELU);
	} //caso nao simetrico
	else{
		TPZSkylineStructMatrix full(fCmesh);
		an.SetStructuralMatrix(full);
		step.SetDirect(ELDLt);
	}//caso simetrico
			
	an.SetSolver(step);
	an.Run();
	static int count = 0;
	if (count%2) {
		//testing code
		TesteInterface(fCmesh, an.Solution());
	}
	
	count++;
	
	/*	TPZAutoPointer<TPZGuiInterface> gui;
	 TPZFMatrix<REAL> rhs, residual;
	 TPZAutoPointer<TPZMatrix> test = full.CreateAssemble(rhs, gui );
	 test->MultAdd(an.Solution(),rhs,residual,1,-1);
	 REAL resnorm = Norm(residual);
	 std::cout << "residual norm " << resnorm << std::endl;
	 */	
	TPZAutoPointer<TPZMaterial> mat1 = fCmesh->FindMaterial(mat1Id);
	TPZAutoPointer<TPZMaterial> mat2 = fCmesh->FindMaterial(mat2Id);
	//TPZElasticityAxiMaterial * aximat1 = dynamic_cast<TPZElasticityAxiMaterial*>(mat1.operator->());
	//TPZElasticityAxiMaterial * aximat2 = dynamic_cast<TPZElasticityAxiMaterial*>(mat2.operator->());
	
	
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	
#ifdef VTK
	TPZManVector<std::string,10> scalnames(3), vecnames(0);
	scalnames[0] = "Sigmarr";
	scalnames[1] = "Sigmazz";
	scalnames[2] = "Sigmatt";
	
	std::string plotfile("saidaT.vtk");
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
#endif	
}


///---------------------------CalcCortMomento()-------------------------------------------------
std::ofstream outf("TensorTensoes.txt");
TPZVec<REAL> tools::CalcCortMomento(TPZCompMesh  *malha)
{
//	REAL QalphaR = 0.;
//	REAL QalphaL = 0.;
//	REAL T1zR = 0.;
//	REAL T1zL = 0.;
//	REAL T1z = 0.;
	static int count = 0;
	std::string estr1,estr2;
#ifdef LOG4CXX
	{
		std::stringstream sout1, sout2,sout;
		sout1 << "Mat1Ind"<< count;
		sout2 << "Mat2Ind" << count;
		estr1 = sout1.str();
		estr2 = sout2.str();
		count++;
		sout << estr1 << " = {};\n";
		sout << estr2 << " = {};\n";
		LOGPZ_DEBUG(loggertensoes,sout.str())
	}
#endif
	outf<<endl;
	outf << "\n=====================  DADOS DAS TENSOES ===============================\n";
	
	int el;
	const int nelem = malha->NElements();
	TPZVec<REAL> Resultados(9,0.);
	
	///Loop sobre elementos computacionais 
	for(el=0; el < nelem; el++) 
	{
		TPZCompEl * Cel = malha->ElementVec()[el];
		if(!Cel) continue;
//		TPZInterfaceElement *inter = dynamic_cast<TPZInterfaceElement *> (Cel);
		int matId = Cel->Material()->Id();
		
		if (matId == matRefDir1) // if(inter)
		{
			//cout << "--------------------------------------"<<endl;
			//			Cel->Reference()->Print();
			//			cout<<"----------------------------------------"<<endl;
			
			TPZCompElSide cels(Cel,2); //  side 2
			TPZStack<TPZCompElSide> neighbours;
			neighbours.Resize(0);
			cels.EqualLevelElementList(neighbours,1,0);
			int nneig = neighbours.NElements();
			
			//			cout << "Num de Vizinho === " << nneig <<endl;
			for (int jn=0; jn<nneig; jn++) 
			{
				int idneigh =neighbours[jn].Reference().Element()->MaterialId();
				
				TPZFMatrix<REAL> n(3,1,0.);//normal
				TPZFMatrix<REAL> tn(3,1,0.);//tangente
				REAL nr = 0., nz=0., tnr=0.,tnz=0.;
				REAL Tr=0., Tz=0.;
//				REAL w, detJac;
				
				if(idneigh == mat1Id)
				{
					//cout << endl;
					//cout<<"===Id dos vizinhos ==="<< idneigh <<endl;
					//neighbours[jn].Reference().Element()->Print();
					
					
					tn(0,0) = sin(falpha);
					tn(1,0) =  cos(falpha);
					MatrixR(Pi/2.).Multiply(tn,n);
					
					tnr = tn(0,0);	tnz = tn(1,0);
					nr = n(0,0); nz = n(1,0);
					
					TPZGeoEl *gel =Cel->Reference();
					TPZIntPoints *intr = gel->CreateSideIntegrationRule(gel->NSides()-1,14);
					TPZGeoElSide gels(gel,2); // Elemento geom.  e lado
					TPZGeoElSide neigh = neighbours[jn].Reference();
					TPZGeoElSide neighhigh(neigh.Element(),neigh.Element()->NSides() - 1);
					
					//cout << endl;
					//cout<<"===neighhigh ==="<<endl;
					//neighhigh.Element()->Print();												
					
					TPZTransform tr1 = gels.NeighbourSideTransform(neigh);
					TPZTransform tr2 = neigh.SideToSideTransform(neighhigh);
					TPZTransform trf = tr2.Multiply(tr1);
					
					///Regra de integracao comeca aqui
					REAL Tn=0., Ts=0.;
					REAL r=0.,z=0., Zc=0.;
					Zc = fRc/tan(falpha);
					int in;
					int nn = intr->NPoints();
					for(in=0; in<nn; in++)
					{
						REAL w, detJac;
						TPZManVector<REAL,3> loc(1),loc2d(2), loc2(2),loc3(2), xco(3), xco1(3), xco2(3),xco3(3); 
						TPZFNMatrix<10> jac(2,2),invJac(2,2),axes(3,3);
						TPZManVector<REAL,5> solsigr(1), solsigz(1), solsigt(1), soltaurz(1);
						
						intr->Point(in,loc,w);
						gel->Jacobian(loc,jac,axes,detJac,invJac);
						trf.Apply(loc,loc2d);
						TPZCompEl *celneigh = neighbours[jn].Element();
						
						//celneigh->Solution(data, var, sol);
						celneigh->Solution(loc2d,3,solsigr);
						celneigh->Solution(loc2d,4,solsigz);
						celneigh->Solution(loc2d,5,solsigt);
						celneigh->Solution(loc2d,6,soltaurz);
						gel->X(loc,xco);
						celneigh->Reference()->X(loc2d,xco2);
						r = xco2[0];
						z = xco2[1];
						
#ifdef LOG4CXX
						{
							std::stringstream outfile;
							outfile.precision(15);
							outfile<<endl;
							outfile << "(*====== DADOS DAS TENSOES =================\n X \n nr nz \n r w detjac alpha Rc \n Tensor *)\n";
							outfile << "AppendTo[" << estr1 << ",{";
							outfile << "\n{" << xco[0] <<" , "<<xco[1] << "}," << endl;
							outfile << "{ " << nr << " , " << nz << " }," << endl;
							outfile << "{ " << r << " , " << w << " , " << fabs(detJac) << " , " << falpha << " , " << fRc << "}," << endl;
							outfile << "{ {" << solsigr << " , " << soltaurz << " , 0. }," << endl;
							outfile << "{ " << soltaurz << " , " << solsigz << ", 0. }," << endl;
							outfile << "{ 0. , 0. , " << solsigt << "} }" << endl;
							outfile << "}];\n";
							LOGPZ_DEBUG(loggertensoes,outfile.str());
						}
#endif
						
						//-------------arquivo textos do tensor das tensoes
						{	
							//QalphaR +=((solsigr[0]*nr + soltaurz[0]*nz)*tnr + (solsigz[0]*nz + soltaurz[0]*nr)*tnz)*xco2[0]*w*fabs(detJac)/fRc;
							//							T1zR +=(soltaurz[0]*nr + solsigz[0]*nz)*xco2[0]*w*fabs(detJac);
							//									
							//							outf<<endl;
							//							outf <<" Id do material ==> "<< celneigh->Reference()->MaterialId()<<endl;
							//							outf <<"coordX = {"<< xco[0] << ", "<<xco[1] <<"}"<< endl;
							//							outf <<"coordX2 = {"<< xco2[0] << ", "<<xco2[1] <<"}"<< endl;
							//							outf <<"peso = "<< w <<"\t detJac = "<<detJac<<endl;
							//							outf <<"TR = {{"<<solsigr<<", "<<soltaurz<<", "<< 0.<<"}\n"
							//							<<"    {"<<soltaurz<<", "<<solsigz <<", "<< 0.<<"}\n"
							//							<<"    {"<<0. <<", "<< 0.<<", "<< solsigt<<"}}"<<endl;
							//							outf << "normal = {"<<nr<<", "<<nz<<"}"<<endl;
							//							outf << "t = {"<<tnr<<", "<<tnz<<"}"<<endl;
							//							outf << "QR = "<<QalphaR <<endl;
							//							outf << "T1zR = "<<T1zR <<endl;
						}
						//----------------------------------------------------------------------------
						
						
						Tr = solsigr[0]*nr + soltaurz[0]*nz;
						Tz = soltaurz[0]*nr + solsigz[0]*nz;
						
						Tn = Tr*cos(falpha) - Tz*sin(falpha);
						Ts = Tr*sin(falpha) + Tz*cos(falpha);
						
						
						Resultados[0] += (1./fRc)*(Tz*(r-fRc) - Tr*(z-Zc))*r*w*fabs(detJac); //Momento a direita
						Resultados[1] += (1./fRc)*(Tr*sin(falpha) + Tz*cos(falpha))*r*w*fabs(detJac); //Cortante a direita
						
						REAL test = (soltaurz[0]*nr + solsigz[0]*nz)*xco2[0]*w*fabs(detJac)*M_PI;//teste T1z do philippe
						//std::cout << "Left contribution " << test << std::endl;
						test += Resultados[6];
						Resultados[6] = test;
						
					}///fim da regra de integracao
					delete intr;
				}//if: matId1d
				
				else{
					//cout << endl;
					//cout<<"===Id dos vizinhos ==="<< idneigh <<endl;
					//neighbours[jn].Reference().Element()->Print();
					
					tn(0,0) = -sin(falpha);
					tn(1,0) = -cos(falpha);
					MatrixR(Pi/2.).Multiply(tn,n);
					
					tnr = tn(0,0);	tnz = tn(1,0);
					nr = n(0,0); nz = n(1,0);
					
					TPZGeoEl *gel =Cel->Reference();
					TPZIntPoints *intr = gel->CreateSideIntegrationRule(gel->NSides()-1,14);
					TPZGeoElSide gels(gel,2); // Elemento geom.  e lado
					TPZGeoElSide neigh = neighbours[jn].Reference();
					TPZGeoElSide neighhigh(neigh.Element(),neigh.Element()->NSides() - 1);
					
					//cout << endl;
					//cout<<"===neighhigh ==="<<endl;
					//neighhigh.Element()->Print();												
					
					TPZTransform tr1 = gels.NeighbourSideTransform(neigh);
					TPZTransform tr2 = neigh.SideToSideTransform(neighhigh);
					TPZTransform trf = tr2.Multiply(tr1);
					
					///Regra de integracao comeca aqui
					REAL Tn=0., Ts=0.;
					REAL r=0.,z=0., Zc=0.;
					Zc = fRc/tan(falpha);
					int in;
					for(in=0; in<intr->NPoints(); in++)
					{
						REAL w, detJac;
						TPZVec<REAL> loc(1),loc2d(2), loc2(2),loc3(2), xco(3), xco1(3), xco2(3),xco3(3); 
						TPZFMatrix<REAL> jac(2,2),invJac(2,2),axes(3,3);
						TPZVec<REAL> solsigr(1), solsigz(1), solsigt(1), soltaurz(1);
						
						intr->Point(in,loc,w);
						gel->Jacobian(loc,jac,axes,detJac,invJac);
						trf.Apply(loc,loc2d);
						TPZCompEl *celneigh = neighbours[jn].Element();
						
						//celneigh->Solution(data, var, sol);
						celneigh->Solution(loc2d,3,solsigr);
						celneigh->Solution(loc2d,4,solsigz);
						celneigh->Solution(loc2d,5,solsigt);
						celneigh->Solution(loc2d,6,soltaurz);
						gel->X(loc,xco);
						celneigh->Reference()->X(loc2d,xco2);
						r = xco2[0];
						z = xco2[1];
						
#ifdef LOG4CXX
						{
							std::stringstream outfile;
							outfile.precision(15);
							outfile<<endl;
							outfile << "(*====== DADOS DAS TENSOES =================\n X \n r w detjac alpha Rc \n Tensor *)\n";
							outfile << "AppendTo[" << estr2 << ",{";
							outfile << "\n{" << xco[0] <<" , "<<xco[1] << "}," << endl;
							outfile << "{ " << nr << " , " << nz << " }," << endl;
							outfile << "{ " << r << " , " << w << " , " << fabs(detJac) << " , " << falpha << " , " << fRc << "}," << endl;
							outfile << "{ {" << solsigr << " , " << soltaurz << " , 0. }," << endl;
							outfile << "{ " << soltaurz << " , " << solsigz << ", 0. }," << endl;
							outfile << "{ 0. , 0. , " << solsigt << "} }" << endl;
							outfile << "}];\n";
							LOGPZ_DEBUG(loggertensoes,outfile.str());
						}
#endif
						
						
						//-------------arquivo textos do tensor das tensoes
						{	
							//QalphaL +=((solsigr[0]*nr + soltaurz[0]*nz)*tnr + (solsigz[0]*nz + soltaurz[0]*nr)*tnz)*xco2[0]*w*fabs(detJac)/fRc;
							//							T1zL +=(soltaurz[0]*nr + solsigz[0]*nz)*xco2[0]*w*fabs(detJac);
							//							
							//							outf<<endl;
							//							outf <<" Id do material ==> "<< celneigh->Reference()->MaterialId()<<endl;
							//							outf <<"coordX = {"<< xco[0] << ", "<<xco[1] <<"}"<< endl;
							//							outf <<"coordX2 = {"<< xco2[0] << ", "<<xco2[1] <<"}"<< endl;
							//							outf <<"peso = "<< w <<"\t detJac = "<<detJac<<endl;
							//							outf <<"TL = {{"<<solsigr<<", "<<soltaurz<<", "<< 0.<<"}\n"
							//							<<"    {"<<soltaurz<<", "<<solsigz <<", "<< 0.<<"}\n"
							//							<<"    {"<<0. <<", "<< 0.<<", "<< solsigt<<"}}"<<endl;
							//							outf << "normal = {"<<nr<<", "<<nz<<"}"<<endl;
							//							outf << "t = {"<<tnr<<", "<<tnz<<"}"<<endl;
							//							outf << "QL = "<<QalphaL <<endl;
							//							outf << "T1zL = "<<T1zL <<endl;
						}
						//----------------------------------------------------------------------------
						
						Tr = solsigr[0]*nr + soltaurz[0]*nz;
						Tz = soltaurz[0]*nr + solsigz[0]*nz;
						
						Tn = Tr*cos(falpha) - Tz*sin(falpha);
						Ts = Tr*sin(falpha) + Tz*cos(falpha);
						
						
						Resultados[2] += (1./fRc)*(Tz*(r-fRc) - Tr*(z-Zc))*r*w*fabs(detJac); //Momento a esquerda
						Resultados[3] += (1./fRc)*(Tr*sin(falpha) + Tz*cos(falpha))*r*w*fabs(detJac); //Cortante a esquerda
						REAL teste = (soltaurz[0]*nr + solsigz[0]*nz)*xco2[0]*w*fabs(detJac)*M_PI;//teste do philippe
						//std::cout << "Right contribution " << teste << std::endl;
						teste += Resultados[7];
						Resultados[7]  = teste;
						//outras formulas
						//Resultados[0] += (-1./fRc)*Tn*((r-fRc)/sin(falpha))*r*w*fabs(detJac); //Momento
						//Resultados[1] += (1./fRc)*Ts*r*w*fabs(detJac); //Cortante
					}///fim da regra de integracao
					
					delete intr;
				}//else
				
			}//for: jn
		}//if: matId=matRefDir1
		Resultados[4] = (Resultados[0] - Resultados[2])/2.; //momento medio
		Resultados[5] = (Resultados[1] - Resultados[3])/2.; //cortante medio
		Resultados[8] = (Resultados[6] - Resultados[7]);//teste do philippe
	}//for: el
	return Resultados;
}


///----------------------------- PrintInterface()--------------------------------------------------
void tools::PrintInterface(TPZCompMesh  *malha)
{
	int el;
	const int nelem = malha->NElements();
	
	TPZInterfaceElement *face;
	
	///Loop sobre elementos computacionais 
	for(el=0; el < nelem; el++) 
	{
		TPZCompEl * Cel = malha->ElementVec()[el];
		if(!Cel) continue;
		int matId = Cel->Material()->Id();
		
		face = dynamic_cast<TPZInterfaceElement *> (Cel);
		
		if(face){
			
			int index = Cel->Index();
			cout<<endl;
			cout << "Ok ... O elemento de Index "<< index <<" e matId = "<< matId<<" eh de  interface"<<endl;
			Cel->Reference()->Print();
			cout <<endl;
		}
	}if(!face) cout << "nao tem interface"<< endl;
}

///-------------------------TesteInterface()-----------------------------------------
void tools::TesteInterface(TPZCompMesh *cmesh, TPZFMatrix<REAL> &solution)
{
	TPZSkylineStructMatrix full(cmesh);
	std::set<int> materialIds;
	//materialIds.insert(mat1Id);
	materialIds.insert(mat3Id);
	full.SetMaterialIds(materialIds);
	
	TPZAutoPointer<TPZGuiInterface> gui;
	TPZFMatrix<REAL> rhs, residual;
	TPZAutoPointer<TPZMatrix<REAL> > test = full.CreateAssemble(rhs, gui );
	test->MultAdd(solution,rhs,residual,1,-1);
	//	REAL resnorm = Norm(residual);
	//	std::cout << "residual norm " << resnorm << std::endl;
	
	REAL accvertical = 0.;
	std::set<int> corner1;
	CornerConnects(cmesh,corner1,mat1Id);
	std::set<int>::iterator it;
	for (it=corner1.begin(); it!=corner1.end(); it++) {
		int seqnum = cmesh->ConnectVec()[*it].SequenceNumber();
		int eqnum = cmesh->Block().Position(seqnum)+1;
		accvertical += residual(eqnum,0);
	}
	std::cout << "Accumulated vertical residual " << accvertical << std::endl;
}


///---------------------- CornerConnects()-----------------------------------------------
void tools::CornerConnects(TPZCompMesh *cmesh, std::set<int> &indices, int matid)
{
	int nel = cmesh->NElements();
	int iel;
	for (iel=0; iel<nel; iel++) {
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		int matidgel = gel->MaterialId();
		if (matidgel != matid) {
			continue;
		}
		int ncorners = gel->NCornerNodes();
		int ic;
		for (ic=0; ic<ncorners; ic++) {
			indices.insert(cel->ConnectIndex(ic));
		}
	}
}



