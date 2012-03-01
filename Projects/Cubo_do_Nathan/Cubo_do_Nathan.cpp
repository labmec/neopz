/*
 *  Cubo_do_Nathan.cpp
 *  PZ
 *
 *  Created by Nathan Shauer on 6/14/11.
 *  Copyright 2011 Unicamp. All rights reserved.
 *
 */
#include <iostream>
#include "pzelast3d.h"
#include "tpzcube.h"
#include "pzgmesh.h"
#include "TPZGeoCube.h"
#include "pzviscoelastic.h"
#include "pzstepsolver.h"
#include "pzseqsolver.h"
#include "tpzautopointer.h"
#include "pzgmesh.h"
#include "tpzarc3d.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzskylstrmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
//#include "pzgeopoint.h"
//#include "pzvec.h"
#include "pzgeotetrahedra.h"
#include "pzpostprocanalysis.h"
#include "pzinterpolationspace.h"



#include <iostream>
#include <cstdlib>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("main"));
#endif

#define MAIN3



void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh);

TPZGeoMesh *MalhaCubo();

void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

#ifdef MAINELA
int main()
{
	// Cria matriz rigidez de um cubo de -1 a 1 para material elástico
	int nummat = 1;
	REAL Ela = 1000, poisson = 0.2;
	TPZVec<REAL> force(3,0.);
	TPZElasticity3D elast3d(nummat, Ela, poisson, force);
	TPZMaterialData cubodata;
	TPZVec <REAL> pt(3,0.);
	TPZFMatrix phi(3,1,0.);
	TPZFMatrix dphi(3,3,0.);
	dphi(0,0) = 1;
	dphi(1,1) = 1;
	dphi(2,2) = 1;
	cubodata.phi = phi;
	cubodata.dphix = dphi;
	cubodata.x = pt;
	TPZFMatrix ek(9,9,0), ef(9,1,0);
	REAL weight = 8.;
	elast3d.Contribute(cubodata, weight, ek, ef);
	ek.Print("ek: ");	
}
#endif

#ifdef MAIN1
int main()
{
	int nummat = 1;
	REAL Ela = 1000, poisson = 0.2;
	TPZVec<REAL> force(3,0.);
	TPZMaterialData cubodata;
	TPZVec <REAL> pt(3,0.);
	TPZFMatrix phi(3,1,0.);
	TPZFMatrix dphi(3,3,0.);
	dphi(0,0) = 1;
	dphi(1,1) = 1;
	dphi(2,2) = 1;
	cubodata.phi = phi;
	cubodata.dphix = dphi;
	cubodata.x = pt;
	TPZFMatrix ek(9,9,0), ef(9,1,0);
	REAL weight = 8.;
	
	// Cria matriz rigidez de um cubo de -1 a 1 para material viscoelástico
	int viscoid = 2;
	int numsteps = 1000;
	//TPZMatWithMem <TPZFMatrix, TPZElasticity3D> viscomem;
	REAL lambdaE = 277.778, muE = 416.667, lambdaV = 0, muV = 0, alphaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alphaT = 0.01;
	ek.Zero();
	ef.Zero();
	TPZViscoelastic viscomat(viscoid, Ela, poisson, lambdaV, muV, alphaT, force);
	int index;
	TPZFNMatrix<6> qsi(6,1,0.);
	viscomat.SetDefaultMem(qsi);
	index = viscomat.PushMemItem();
	cubodata.intPtIndex = index;
	ef(0,0) = 1.*4.+1.*4.;
	viscomat.Contribute(cubodata, weight, ek, ef);
	TPZAutoPointer <TPZMatrix> matrix = new TPZFMatrix(ek);	
	
	//Solve the system ek.result=ef
	TPZStepSolver solver;
	solver.SetMatrix(matrix);
	solver.SetDirect(ECholesky);
	TPZFMatrix result(9,1,0.);
	solver.Solve(ef, result);
	// system solved
	
	viscomat.SetUpdateMem();
	int i,j,k;	
	TPZFMatrix dsol(3,3,0.);
	for (k = 0 ; k < 3 ; k++)
	{
		for (i = 0 ; i < 3 ; i++)
		{
			for (j = 0 ; j < 3 ; j++)
			{
				dsol(j,k) += result(3*i+k,0)*dphi(j,i);
			}
		}	
	}
	cubodata.dsol = dsol;
	TPZFMatrix Stress;
	viscomat.ComputeStressTensor(Stress, cubodata); // Nao pode dar update no qsi para calcular a tensao!
	TPZFMatrix Strain;
	viscomat.ComputeStrainTensor(Strain, dsol);
	TPZVec <REAL> Strainxx(numsteps+1), Stressxx(numsteps+1);
	Strainxx[0] = Strain(0,0);
	Stressxx[0] = Stress(0,0);
	viscomat.UpdateQsi(cubodata);
	
	
	int iv;
	for(iv = 0 ; iv < numsteps ; iv++)
	{
		ek.Zero();
		ef.Zero();
		ef(0,0) = 1.*4.+1.*4.;
		viscomat.Contribute(cubodata, weight, ek, ef);

		matrix = new TPZFMatrix(ek);
		//ek.Print("ek");
		//matrix->Print("matrix");
		
		//Solve the system ek.result=ef
		solver.ResetSolver();
		solver.SetMatrix(matrix);
		solver.SetDirect(ECholesky);
		result.Zero();
		solver.Solve(ef, result);
		// system solved
		//result.Print("result");

		int i,j,k;	
		TPZFMatrix dsol(3,3,0.);
		for (k = 0 ; k < 3 ; k++)
		{
			for (i = 0 ; i < 3 ; i++)
			{
				for (j = 0 ; j < 3 ; j++)
				{
					dsol(j,k) += result(3*i+k,0)*dphi(j,i);
				}
			}	
		}
		cubodata.dsol = dsol;
		Stress.Zero();
		viscomat.ComputeStressTensor(Stress, cubodata); // Nao pode dar update no qsi para calcular tensao!
		Stressxx[iv+1] = Stress(0,0);
		Strain.Zero();
		viscomat.ComputeStrainTensor(Strain, dsol);
		Strainxx[iv+1] = Strain(0,0);
		
		viscomat.UpdateQsi(cubodata);
	}
	cout << "Strain" << endl;
	cout << Strainxx << endl;
	cout << "Stress" << endl;
	cout << Stressxx << endl;

	ofstream out("strain.txt");
	for (int is = 0 ; is < Strainxx.NElements() ; is++)
	{
		out << Strainxx[is] << endl;
	}
	ofstream out2("stresscte.txt");
	for (int iss = 0 ; iss < Strainxx.NElements() ; iss++)
	{
		out2 << Stressxx[iss] << endl;
	}
	
	return 0;
}
#endif


#ifdef MAIN2
int main()
{
	int nummat = 1;
	REAL Ela = 1000, poisson = 0.2;
	TPZVec<REAL> force(3,0.);
	TPZMaterialData cubodata;
	TPZVec <REAL> pt(3,0.);
	TPZFMatrix phi(3,1,0.);
	TPZFMatrix dphi(3,3,0.);
	dphi(0,0) = 1;
	dphi(1,1) = 1;
	dphi(2,2) = 1;
	cubodata.phi = phi;
	cubodata.dphix = dphi;
	cubodata.x = pt;
	TPZFMatrix ek(9,9,0), ef(9,1,0);
	REAL weight = 8.;
	
	// Cria matriz rigidez de um cubo de -1 a 1 para material viscoelástico
	int viscoid = 2;
	int numsteps = 1000;
	//TPZMatWithMem <TPZFMatrix, TPZElasticity3D> viscomem;
	REAL lambdaE = 277.778, muE = 416.667, lambdaV = 0, muV = 0, alphaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alphaT = 0.01;
	ek.Zero();
	ef.Zero();
	TPZViscoelastic viscomat(viscoid, Ela, poisson, lambdaV, muV, alphaT, force);
	int index;
	TPZFNMatrix<6> qsi(6,1,0.);
	viscomat.SetDefaultMem(qsi);
	index = viscomat.PushMemItem();
	cubodata.intPtIndex = index;
	TPZFMatrix Dsol(3,3,0.);
	Dsol(0,0) = 0.001;
	cubodata.dsol = Dsol;
	TPZFMatrix Stress;
	viscomat.ComputeStressTensor(Stress,cubodata);	
	TPZVec <REAL> Stressxx(numsteps+1);
	Stressxx[0] = Stress(0,0);

	viscomat.Contribute(cubodata, weight, ek, ef);
	TPZAutoPointer <TPZMatrix> matrix = new TPZFMatrix(ek);	
	
	viscomat.SetUpdateMem();
	viscomat.UpdateQsi(cubodata);
		
	int iv;
	for(iv = 0 ; iv < numsteps ; iv++)
	{
		ek.Zero();
		ef.Zero();
		viscomat.Contribute(cubodata, weight, ek, ef);
		Stress.Zero();
		viscomat.ComputeStressTensor(Stress, cubodata);
		Stressxx[iv+1] = Stress(0,0);
		
		viscomat.UpdateQsi(cubodata);
	}
	
	cout << Stressxx << endl;
	
	ofstream out("stress.txt");
	for (int is = 0 ; is < Stressxx.NElements() ; is++)
	{
		out << Stressxx[is] << endl;
	}
	
}
#endif

#ifdef MAIN3
int main()
{
	InitializePZLOG();
	TPZGeoMesh *gmesh = 0;
	gmesh = MalhaCubo();
	int porder = 1;
	//TPZCompEl::SetgOrder(porder);
	
	// Rodando com p=1
	TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
	
	InsertViscoElasticity(cmesh);

	cmesh->SetDefaultOrder(porder);

	//cmesh->SetAllCreateFunctionsContinuousWithMem();
	cmesh->SetAllCreateFunctionsContinuous(); //elastic

	cmesh->AutoBuild();

	TPZSkylineStructMatrix skylstruct(cmesh);
	TPZStepSolver step;
	step.SetDirect(ECholesky);
	TPZAnalysis an(cmesh);
	an.SetStructuralMatrix(skylstruct);
	an.SetSolver(step);
	an.Run();	
	
	// Rodando com p=2
	TPZAutoPointer<TPZCompMesh> cmesh2 = new TPZCompMesh(gmesh);
	
	InsertViscoElasticity(cmesh2);
	
	cmesh2->SetDefaultOrder(2);
	//cmesh->SetAllCreateFunctionsContinuousWithMem();
	cmesh2->SetAllCreateFunctionsContinuous(); //elastic
	cmesh2->AutoBuild();
	
	TPZSkylineStructMatrix skylstruct2(cmesh2);
	TPZStepSolver step2;
	step2.SetDirect(ECholesky);
	TPZAnalysis an2(cmesh2);
	an2.SetStructuralMatrix(skylstruct2);
	an2.SetSolver(step2);
	an2.Assemble();	
	//an2.Run();	
	
	//an.Solution().Print("Solution");

#ifdef LOG4CXX
	{
		std::stringstream str;
		an.Solution().Print("Solution",str);
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
	
#ifdef LOG4CXX
	{
		std::stringstream str;
		an2.Solution().Print("Solution",str);
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
	
	/*
	std::map<int ,TPZAutoPointer<TPZMaterial> > materialmap(cmesh->MaterialVec());
	std::map<int ,TPZAutoPointer<TPZMaterial> >::iterator it;
	for (it = materialmap.begin(); it != materialmap.end() ; it++) 
	{
		TPZAutoPointer<TPZMaterial> mat = it->second;
		TPZViscoelastic *vmat = dynamic_cast< TPZViscoelastic *> (mat.operator->());
		if(vmat)
		{
			vmat->SetUpdateMem();
		}
	}
	*/ 
	
	
	TPZAdmChunkVector<TPZCompEl *> ElementVec = cmesh->ElementVec();
	int nel = ElementVec.NElements();
	int iel;
	for (iel = 0; iel < nel; iel++) 
	{
		TPZCompEl *el = ElementVec[iel];
		TPZInterpolationSpace *insp = dynamic_cast<TPZInterpolationSpace *> (el);
		insp->PRefine(2);
	}
	cmesh->ExpandSolution();
	/*
	// Resolvendo p=2 com p=1 interpolada
  TPZAutoPointer<TPZMatrix> mat = an2.Solver().Matrix();
	TPZFMatrix carga = an2.Rhs();
	mat->Print("rigidez");
	carga.Print("carga");
	mat->SolveDirect(carga,ECholesky);
	carga.Print("Solucao");
	 */
	
	//Olhando residuo
	TPZAutoPointer<TPZMatrix> mat = an2.Solver().Matrix();
	TPZFMatrix carga = an2.Rhs();
	const TPZFMatrix sol = cmesh->Solution();
	TPZFMatrix result;
	mat->Multiply(sol,result);
	
	if (cmesh->Solution().Rows() != result.Rows()) 
	{
		std::cout << "Solucoes de tamanho diferente" << std::endl;
		DebugStop();
	}
	TPZFMatrix sub(cmesh->Solution().Rows(),cmesh->Solution().Cols());
	for (int is = 0 ; is < cmesh->Solution().Rows() ; is++) 
	{
		sub(is,0) = result(is,0) - carga(is,0);
	}
	
#ifdef LOG4CXX
	{
		std::stringstream str;
		sub.Print("Subtracao",str);
		LOGPZ_DEBUG(logger,str.str());
	}
#endif
	

  int dimension = 3;
	int resolution = 1;
  TPZVec <std::string> vecnames(4), scalnames(0);
	cmesh->Solution() = sub;
	
//  scalnames[0] = "StressX";
	vecnames[0] = "Displacement";
	vecnames[1] = "DisplacementX";
	vecnames[2] = "DisplacementY";
	vecnames[3] = "DisplacementZ";
	//scalnames[0] = "ViscoStressX";
	std::string plotfile("cubinho.vtk");
	
	an.DefineGraphMesh(dimension, scalnames, vecnames, plotfile);
	an.PostProcess(resolution);

	std::string plotfile2("cubinho2.vtk");
	an2.DefineGraphMesh(dimension, scalnames, vecnames, plotfile2);
	an2.PostProcess(resolution);
	
	return 0;
	
	TPZPostProcAnalysis postan(&an);
	TPZVec <int> matids(3);
	matids[0] = 1;
	TPZVec <std::string> varName(3);
	varName[0] = "Displacement";
	varName[1] = "PrincipalStrain";
	varName[2] = "ViscoStressX";
	postan.SetPostProcessVariables(matids, varName);
	TPZFStructMatrix bobo(postan.Mesh());
	postan.SetStructuralMatrix(bobo);
 	postan.TransferSolution();
	std::string postplot("discont.vtk");
	vecnames.resize(2);
	vecnames[0] = "Displacement";
	vecnames[1] = "PrincipalStrain";
	scalnames.resize(1);
	scalnames[0] = "ViscoStressX";

	postan.DefineGraphMesh(dimension, scalnames, vecnames, postplot);
	postan.PostProcess(resolution);

    for (int istep = 0 ; istep < 5 ; istep++)
    {
      an.AssembleResidual();
      an.Solve();
			postan.TransferSolution();
		  postan.PostProcess(resolution);
    }
        
	//an.Solution().Print("Solution",std::cout);
    //an.AssembleResidual();
    //an.Solve();
	//an.Solution().Print("Solution",std::cout);
	
	//void TPZAnalysis::DefineGraphMesh(int dim, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile) 
//0.00207065
		
	//postan.Run();

	//postan.Solution().Print("Deus do ceu:", std::cout);
		
	//TPZAutoPointer <TPZGuiInterface> toto;
	//TPZFMatrix rhs;
	//skylstruct.Assemble(rhs,toto);

	return 0;
}

#endif


void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1, neumann = 1, mixed = 2;
//	int dirichlet = 0;
	int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5, dirp2 = -6;
	TPZManVector<REAL> force(3,0.);
	//force[1] = 0.;
	REAL Ela = 1000, poisson = 0.; 
	REAL lambdaV = 0, muV = 0, alphaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alphaT = 0.01;	
    
    
	//TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, Ela, poisson, lambdaV, muV, alphaT, force);
	TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, Ela, poisson, force);

	//TPZFNMatrix<6> qsi(6,1,0.);
	//viscoelast->SetDefaultMem(qsi); //elast
	//int index = viscoelast->PushMemItem(); //elast
	TPZAutoPointer<TPZMaterial> viscoelastauto(viscoelast);
	mesh->InsertMaterialObject(viscoelastauto);
		
	// Neumann em x = 1;
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
	val2(0,0) = 1.;
	TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto4(bc4);
	mesh->InsertMaterialObject(bcauto4);
	
	// Neumann em x = -1;
	val2(0,0) = -1.;
	TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto5(bc5);
	mesh->InsertMaterialObject(bcauto5);
	
	val2.Zero();
	// Dirichlet em -1 -1 -1 xyz;
	val1(0,0) = 1.;
	val1(1,1) = 1.;
	val1(2,2) = 1.;
	TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto1(bc1);
	mesh->InsertMaterialObject(bcauto1);
	
	// Dirichlet em 1 -1 -1 yz;
	val1(0,0) = 0.;
	val1(1,1) = 1.;
	val1(2,2) = 1.;
	TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto2(bc2);
	mesh->InsertMaterialObject(bcauto2);
	
	// Dirichlet em 1 1 -1 z;
	val1(0,0) = 0.;
	val1(1,1) = 0.;
	val1(2,2) = 1.;
	TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto3(bc3);
	mesh->InsertMaterialObject(bcauto3);
	
}



using namespace std;

TPZGeoMesh *MalhaCubo()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "../cube1.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZManVector <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
		
		//(Node, idbcnode, *gMesh);	
		
	}
	
	{
		read.close();
		read.open(FileName.c_str());
				
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int el;
		int neumann1 = -4, neumann2 = -5, dirp2 = -6;
        int index = 0;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index = el;

			//TPZGeoEl * tetra = gMesh->CreateGeoElement(ETetraedro, TopolTetra, matElId, index);
			//TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
//<<<<<<< .mine

		}
		
		//gMesh->BuildConnectivity();
		
		// Colocando as condicoes de contorno
		for(el=0; el<numelements; el++)
		{

 
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			
			// na face x = 1
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == 1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann1);		
			}
			
			// Na face x = -1
			ncoordzVec.Resize(0);
			sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];

				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == -1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann2);		
			}

		}
	
		//gMesh->BuildConnectivity();
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
		gMesh->BuildConnectivity();
	}
	
	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
	//	for (int el = 0; el < numnodes-1; el++) 
	//	{
	//		Nodefind[el] = gMesh->NodeVec()[el];
	//
	//	}
	//	Nodefind.Print(std::cout);
	//	std::cout.flush();
	
	//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
	//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
	
	ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	
	std::ofstream out("Cube.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
	
	return gMesh;
	
}


/*
TPZGeoMesh *MalhaCubo()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName;
	FileName = "../cube1.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZVec <int> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	int idbcnode = -2;
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
		
		//(Node, idbcnode, *gMesh);     
		
	}
	
	{
		read.close();
		read.open(FileName.c_str());
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int el;
		int neumann1 = -4, neumann2 = -5;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index;
			//TPZGeoEl * tetra = gMesh->CreateGeoElement(ETetraedro, TopolTetra, matElId, index);
			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
			
			// Colocando as condicoes de contorno
			TPZVec <TPZGeoNode> Nodefinder(4);
			TPZVec <REAL> nodecoord(3);
			// na face x = 1
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				cout << TopolTetra[i] << endl;
				Nodefinder[i] = gMesh->NodeVec()[TopolTetra[i]];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == 1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = TopolTetra[i];
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				cout << ncoordzVec << endl;
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann1);         
			}
			
			// Na face x = -1
			ncoordzVec.Resize(0);
			sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				Nodefinder[i] = gMesh->NodeVec()[TopolTetra[i]];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[0] == -1.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = TopolTetra[i];
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,neumann2);         
			}
			
			
		}
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
		gMesh->BuildConnectivity();
	}
	
	// identificando as superficies que terao cond de contorno. Coord z dos 3 nos = 0
	//      for (int el = 0; el < numnodes-1; el++) 
	//      {
	//              Nodefind[el] = gMesh->NodeVec()[el];
	//
	//      }
	//      Nodefind.Print(std::cout);
	//      std::cout.flush();
	
	//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
	//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);
	
	//ofstream arg("malhaPZ.txt");
	//gMesh->Print(arg);
	
	std::ofstream out("Cube.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
	
	return gMesh;
	
}
*/
 
//TESTE Step Solver
/*
 TPZAutoPointer <TPZMatrix> matrixteste = new TPZFMatrix(2,2,0.);
 TPZFMatrix b(2,1,0.);
 matrixteste->operator()(0,0) = 3.;
 matrixteste->operator()(0,1) = 2.;
 matrixteste->operator()(1,0) = 5.;
 matrixteste->operator()(1,1) = -2.;
 b(0,0) = 7.;
 b(1,0) = 1.;
 matrixteste->Print("matrixteste: ");
 b.Print("b:");
 TPZStepSolver solver;
 solver.SetMatrix(matrixteste);
 solver.SetDirect(ELU);
 TPZFMatrix result(2,1,0.);
 result.Print("desloc:");
 solver.Solve(b, result);
 */
//TESTE Step Solver

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
	// look for an element/corner node whose distance is close to start
	TPZGeoNode *gn1 = gr->FindNode(x);
	int iel;
	int nelem = gr->ElementVec().NElements();
	TPZGeoEl *gel;
	for (iel = 0; iel<nelem; iel++) {
		gel = gr->ElementVec()[iel];
		if(!gel) continue;
		int nc = gel->NCornerNodes();
		int c;
		for (c=0; c<nc; c++) {
			TPZGeoNode *gn = gel->NodePtr(c);
			if (gn == gn1) {
				break;
			}
		}
		if (c<nc) {
			TPZGeoElBC(gel, c, bc);
			return;
		}
	}
}

