/**
 * @file
 * @brief Tests for sub structuration
 * @author Philippe Devloo
 * @since 2006
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzmat1dlin.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzelasmat.h"
#include "TPZVTKGeoMesh.h"

#include "TPZYCMohrCoulombPV.h"


#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

TPZGeoMesh *CreateGMesh(int nel, REAL elsize);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);

int mainEla2D();
int mainSimples();
int mainRascunho();

TPZGeoMesh * GMesh2D();
TPZCompMesh *CMesh2D(TPZGeoMesh *gmesh, const int porder);


using namespace std;

int main(int argc, char *argv[])
{
	return mainEla2D();
}

int mainEla2D()
{
	const int porder = 1;
	TPZGeoMesh *gmesh = GMesh2D();
	ofstream file("gMeshElast.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
	TPZCompMesh *cmesh = CMesh2D(gmesh,porder);
	
	TPZCompEl *cel = cmesh->ElementVec()[0];
	int ncon = cel->NConnects();
	std::cout << "NCon = " << ncon << std::endl;
	int ord = cel->Connect(8).Order();
	std::cout << "ord = " << ord << std::endl;
	int nstate = cel->Connect(8).NState();
	std::cout << "nstate = " << nstate << std::endl;
	
	TPZAnalysis an(cmesh);
	TPZStepSolver <REAL> step;
	step.SetDirect(ECholesky);
	an.SetSolver(step);
	TPZSkylineStructMatrix skylstr(cmesh);
	an.SetStructuralMatrix(skylstr);
	an.Run();
	std::cout << cmesh->Solution() << std::endl;
	
	TPZStack <std::string> vecnames,scalnames;
	
	vecnames.Push("displacement");
	vecnames.Push("Strain");
	scalnames.Push("SigmaY");
  an.DefineGraphMesh(2, scalnames, vecnames, "outElast.vtk");
	int resolution = 2;
	an.PostProcess(resolution);
	
	return 0;
}

int mainRascunho()
{
	// TESTE DO MOHRCOULOMBPV
	
	TPZYCMohrCoulombPV *MohrCoulombPV = new TPZYCMohrCoulombPV;
	
	TPZManVector<REAL, 3> sigtrialvec(3,0.), sigprojectvec(3,0.);
	TPZYCMohrCoulombPV::TComputeSequence toto;
	toto.fGamma.resize(1);
	sigtrialvec[0] = 500.;
	TPZManVector<TFad<3,REAL>,3 > sigtrialfad(3), sigprojfad(3);
	for (int i = 0; i < 3; i++) 
	{
		sigtrialfad[i].val() = sigtrialvec[i];
		sigtrialfad[i].fastAccessDx(i) = 1;
	}
	
	//MohrCoulombPV->ReturnMapPlane<REAL>(sigtrialvec, sigprojectvec, toto);
	
	toto.fGamma[0] = 0.;
	//MohrCoulombPV->ReturnMapPlane<TFad<3,REAL> >(sigtrialfad, sigprojfad, toto);
	std::cout << "sigtrialfad:\n" << sigtrialfad << "\nsigprojfad:\n" << sigprojfad << std::endl;
	
	//REAL critplat = MohrCoulombPV->PhiPlane<REAL>(sigprojectvec);
	//std::cout << "CritPlat = " << critplat << std::endl;
	
	//  critplat = MohrCoulombPV->PhiPlane<REAL>(sigmaproject);
	//	std::cout << "CritPlat = " << critplat << std::endl;
	
	
	// Testando Left Edge
	toto.fGamma.Resize(2,0.);
	toto.fGamma.Fill(0.);
	std::cout << "toto.fGamma = " << toto.fGamma << std::endl;
	sigtrialvec[0] = 500.;
	sigtrialvec[1] = 0.;
	sigtrialvec[2] = 0.;	
	sigprojectvec.Fill(0.);
	std::cout << "sigprojectvec = " << sigprojectvec << std::endl;
	//MohrCoulombPV->ReturnMapLeftEdge(sigtrialvec, sigprojectvec, toto);
	//MohrCoulombPV->ReturnMapRightEdge(sigtrialvec, sigprojectvec, toto);
	
	//	template<class T>
	//	bool ReturnMapLeftEdge(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected,
	//												 TComputeSequence &memory)
	
	// Será que está na superfície?
	//    T PhiPlane(typename TPZTensor<T>::TPZDecomposed &sigma) const
	
	
	
	//	template<class T>
	//	bool ReturnMapPlane(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected, 
	//											TComputeSequence &memory)
	
	
	
	// FIM TESTE DO MOHRCOULOMBPV
	
}

int mainSimples(int argc, char *argv[])
{
	
	// Gerar Problema Simples
	
	REAL dom = 1.;
	int nel = 3, pOrder = 1;
	REAL elsize = dom/nel;
	TPZGeoMesh *gmesh = CreateGMesh(nel, elsize);
	
	std::ofstream out("gmesh.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
	
	TPZCompMesh *cmesh = CMesh(gmesh, pOrder);
	
	// Resolvendo o Sistema
	TPZAnalysis an(cmesh);
	TPZSkylineStructMatrix skyl(cmesh);
	an.SetStructuralMatrix(skyl);
	TPZStepSolver<REAL> step;
	step.SetDirect(ECholesky);
	an.SetSolver(step);
	an.Run();
	TPZFMatrix<REAL> solucao;
	solucao = cmesh->Solution();
	std::ofstream outomar("omar.nb");
	solucao.Print("Sol",outomar,EMathematicaInput);
	solucao.Print("Sol");
	
	//Print(const char *name, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	// PosProcessamentoiag
	
	 
	std::cout << "FINISHED!" << std::endl;
}


	 
TPZGeoMesh *CreateGMesh(int nel, REAL elsize)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	
	int nnodes = nel + 1;
	gmesh->NodeVec().Resize(nnodes);
	int mat1d = 1, bc0 = -1, bc1 = -2;
	
	// Colocando nos na malha
	for (int i = 0 ; i < nnodes; i++) 
	{
		const REAL pos = i * elsize;
		TPZVec <REAL> coord(3,0.);
		coord[0] = pos;
		gmesh->NodeVec()[i].SetCoord(coord);
		gmesh->NodeVec()[i].SetNodeId(i);
	}
	
	// Criando Elementos
	TPZVec <int> topol(2), TopolPoint(1);
	int id = 0;
	
	for (int iel = 0; iel < nel; iel++) 
	{
		const int ino1 = iel;
		const int ino2 = iel + 1;
		topol[0] = ino1;
		topol[1] = ino2;
		
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,topol,mat1d,*gmesh);
		id++;
	}
	
	// Cond Contorno esquerda
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc0,*gmesh);
	id++;

	// Cond Contorno Direita
	TopolPoint[0] = nnodes-1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc1,*gmesh);
	
	gmesh->BuildConnectivity();
	
	return gmesh;
}

TPZGeoMesh * GMesh2D()
{
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	const int nnodes = 4;
	gmesh->NodeVec().Resize(nnodes);
	const int matid = 1, bc0 = -1, bc1 = -2;
	
	int ino = 0;
	TPZManVector <REAL,3> coord(3,0.);
	// no 0
	coord[0] = 0.;
	coord[1] = 0.1;
	gmesh->NodeVec()[ino].SetCoord(coord);
	gmesh->NodeVec()[ino].SetNodeId(ino);
	ino++;

	// no 1
	coord[0] = 1.;
	coord[1] = 0.;
	gmesh->NodeVec()[ino].SetCoord(coord);
	gmesh->NodeVec()[ino].SetNodeId(ino);
	ino++;
	
	// no 2
	coord[0] = 1.;
	coord[1] = 1.;
	gmesh->NodeVec()[ino].SetCoord(coord);
	gmesh->NodeVec()[ino].SetNodeId(ino);
	ino++;
	
	// no 3
	coord[0] = 0.;
	coord[1] = 1.;
	gmesh->NodeVec()[ino].SetCoord(coord);
	gmesh->NodeVec()[ino].SetNodeId(ino);
	ino++;
	
	TPZManVector <int,4> TopolQuad(4);
	TPZManVector <int,1> TopolPoint(1);
	TPZManVector <int,2> TopolLine(2);
	
	int id = 0;
	const int nel = 1;

	for (int ip = 0; ip < 4; ip++) {
		TopolQuad[ip] = ip;
	}
	new TPZGeoElRefPattern < pzgeom::TPZGeoQuad > (id, TopolQuad, matid, *gmesh);
	id++;
	
	// Cond Contorno esquerda baixo
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc0,*gmesh);
	id++;

	// Cond Contorno esquerda baixo
	TopolPoint[0] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc0,*gmesh);
	id++;

	
	// Cond Contorno esquerda baixo
	TopolPoint[0] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc1,*gmesh);
	id++;
	
	// Cond Contorno esquerda baixo
	TopolPoint[0] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint> (id,TopolPoint,bc1,*gmesh);
	id++;
	 
	/*
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (id,TopolLine,bc1,*gmesh);
	*/
	gmesh->BuildConnectivity();

	return gmesh;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
	const int dim = 1, k = 1, c = 0, b = 1, f = 1;
	const int matId = 1, bc0 = -1, bc1 = -2;
	const int dirichlet = 0, neumann = 1, mixed = 2;
	TPZFMatrix <REAL> xk(1,1,k), xc(1,1,c), xb(1,1,b), xf(1,1,f);
	
	// Criando material
	TPZMat1dLin *material;
	material = new TPZMat1dLin(matId);
	material->SetMaterial(xk,xc,xb,xf);
	TPZMaterial * mat(material);
	material->NStateVariables();
	
	///criar malha computacional
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
	
	// Inserindo material na malha
	cmesh->InsertMaterialObject(mat);
		
	///Inserir condicao de contorno esquerda
	TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0, dirichlet, val1, val2);
	
	// Condicao de contorno da direita
	TPZMaterial * BCond1 = material->CreateBC(mat, bc1, dirichlet, val1, val2);
	
	//cmesh->SetAllCreateFunctionsContinuous();
	//cmesh->InsertMaterialObject(mat);
	cmesh->InsertMaterialObject(BCond0);
	cmesh->InsertMaterialObject(BCond1);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	//cmesh->AdjustBoundaryElements();
	//cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;    
	
	
}


TPZCompMesh *CMesh2D(TPZGeoMesh *gmesh, const int porder)
{
	const int dim = 2;
	const int matid = 1, bc0 = -1, bc1 = -2;
	const int dirichlet = 0, neumann = 1, mixed = 2;
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(porder);
	cmesh->SetDimModel(dim);

	REAL Elast = 1000.;
	REAL nu = 0.;
	REAL fx(0.),fy(0.);
	
	TPZElasticityMaterial *elast2d = new TPZElasticityMaterial(matid,Elast,nu,fx,fy);
	
	cmesh->InsertMaterialObject(elast2d);
	
	TPZFNMatrix <4> val1(2,2,0.), val2(2,1,0.);
	TPZMaterial *BCcond0 = elast2d->CreateBC(elast2d,bc0,dirichlet,val1,val2);
	cmesh->InsertMaterialObject(BCcond0);
	
	val2(1,0) = -10.;
	TPZMaterial *BCcond1 = elast2d->CreateBC(elast2d,bc1,neumann,val1,val2);
	cmesh->InsertMaterialObject(BCcond1);
	
	cmesh->AutoBuild();
	
	return cmesh;
}