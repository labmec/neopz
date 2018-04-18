/**
 * @file
 * @brief Tests for sub structuration
 * @author Philippe Devloo
 * @since 2006
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "tpzdohrassembly.h"

#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZTimer.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"

#include "pzvtkmesh.h"

#include "pzlog.h"

#include <fstream>
#include <string>

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
static LoggerPtr loggernathan(Logger::getLogger("pz.nathan"));
#endif

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaPredio();
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
REAL Height(TPZGeoMesh *gmesh);
int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height);

using namespace std;
int main1(int argc, char *argv[]);


int main(int argc, char *argv[])
{
    return main1(argc, argv);
}

int main1(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	TPZTimer total;
	total.start();
	std::cout << "COMECA O TEMPO"<< std::endl;
	
	int dimension = 3;
	int dim = 2;
	int maxlevel = 5;
	int sublevel = 3;
	int plevel = 1;
	TPZPairStructMatrix::gNumThreads = 20;
	int numthreads = 20;
	//	tempo.fNumthreads = numthreads;	// alimenta timeTemp com o numero de threads
	TPZGeoMesh *gmesh = 0;
	{
		TPZCompEl::SetgOrder(plevel);
		
		TPZAutoPointer<TPZCompMesh> cmesh;
		
		if(0)
		{
			TPZGenSubStruct sub(dim,maxlevel,sublevel);
			cmesh = sub.GenerateMesh();
			cmesh->SetDimModel(dim);
			gmesh = cmesh->Reference();
		}
		else 
		{
			dim = 3;
			if (1) // Predio Viscoso
			{
				int dimension = 3;
				gmesh = MalhaPredio();
				cmesh = new TPZCompMesh(gmesh);
				cmesh->SetDimModel(3);
				cmesh->SetDefaultOrder(plevel);
				cmesh->SetAllCreateFunctionsContinuousWithMem(dimension);
				InsertViscoElasticity(cmesh);
				cmesh->AutoBuild();
				
			}
			else // Cubo Viscoso
			{
				int dimension = 3;
				gmesh = MalhaCubo();
				cmesh = new TPZCompMesh(gmesh);
				cmesh->SetDimModel(3);
				cmesh->SetDefaultOrder(plevel);
				cmesh->SetAllCreateFunctionsContinuousWithMem(dimension);
				InsertViscoElasticityCubo(cmesh);
				cmesh->AutoBuild();
			}
		}
		
		std::cout << "Numero de equacoes " << cmesh->NEquations() << std::endl;
		
		int numthread_assemble = 20;
		int numthread_decompose = 20;
		TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
		TPZDohrStructMatrix dohrstruct(cmeshauto);
		
		dohrstruct.IdentifyExternalConnectIndexes();
		
		std::cout << "Substructuring the mesh\n";
		//	TPZfTime timetosub; // init of timer
		//REAL height = Height(gmesh);
		//int nsubstruct = SubStructure(cmesh, height/2);
		
		dohrstruct.SubStructure(16);

		//	tempo.ft0sub = timetosub.ReturnTimeDouble();  // end of timer
		//	std::cout << tempo.ft0sub << std::endl;
		
		//	sub.SubStructure();
		
		
		 //Teste Skyline
        /*
		TPZSkylineStructMatrix skyl(cmesh);
		TPZFMatrix<REAL> rhsfake(cmesh->NEquations(),1,0);
		int numsubmesh = cmesh->NElements();
		TPZAutoPointer<TPZGuiInterface> fakegui = new TPZGuiInterface;
		int nel = cmesh->NElements();
		for (int iel = 0 ; iel < nel ; iel++)
		{
			TPZSubCompMesh *subcompmesh = dynamic_cast<TPZSubCompMesh*>(cmesh->ElementVec()[iel]);
			if(subcompmesh)
			{
				subcompmesh->SetAnalysisSkyline(0,0,fakegui);
			}
		}
		TPZMatrix<REAL> *stiff2 = skyl.CreateAssemble(rhsfake, fakegui,numthread_assemble,numthread_decompose);
		*/
		
#ifdef LOG4CXX
		{
			std::stringstream str;
			cmesh->Print(str);
			LOGPZ_DEBUG(logger,str.str());
		}
#endif

		
		dohrstruct.SetNumThreads(numthreads);
		
		TPZAutoPointer<TPZGuiInterface> gui;
		TPZFMatrix<STATE> rhs(cmesh->NEquations(),1,0.);
        
		TPZMatrix<STATE> *matptr = dohrstruct.Create();
		
		 
		dohrstruct.Assemble(*matptr,rhs,gui,numthread_assemble,numthread_decompose);

	
		TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
		TPZAutoPointer<TPZMatrix<STATE> > precond = dohrstruct.Preconditioner();
		
		{
			std::ofstream out("DohrCerta2.txt");
			TPZFMatrix<REAL> Subtract(dohr->Rows(),dohr->Rows()), unitary(dohr->Rows(),dohr->Rows());
			unitary.Identity();
			TPZFMatrix<REAL> result;
			dohr->Multiply(unitary, result);
			result.Print("DohrCerta2", out);
			
		}

/*	
#ifdef LOG4CXX
		{  
			std::ofstream out("DohrErrada.txt"), outRhsCerto("RhsSkyl.txt"), outRhsErrado("RhsDohrmann.txt");
			TPZFMatrix<REAL> Subtract(dohr->Rows(),dohr->Rows()), unitary(dohr->Rows(),dohr->Rows());
			unitary.Identity();
			TPZFMatrix<REAL> result;
			dohr->Multiply(unitary, result);
			std::ofstream out2("Dohr_Certa.txt");
			result.Print("DohrCerta",out2);
			for (int i = 0 ; i < dohr->Rows(); i++) 
			{
				for (int j = 0 ; j < dohr->Rows(); j++) 
				{
					double temp = result(i,j) - stiff2->Get(i,j); 
					if (temp < 1e-10) 
					{
							temp = 0;
					}
					Subtract(i,j) = temp; 
				}
			}
			std::stringstream str;
			result.Print("DohrmannErrada", out);
			stiff2->Print("Skyl",out);
			Subtract.Print("Subtract", out);
			rhsfake.Print("RhsSkyl", outRhsCerto);
			rhs.Print("RhsDohrmann", outRhsErrado);
			LOGPZ_DEBUG(logger,str.str());
		}
#endif
 */
		
        
		int neq = dohr->Rows();
        
		TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
		std::cout << "Numero de equacoes " << neq << std::endl;
        
		TPZStepSolver<STATE> pre(precond);
		pre.SetMultiply();
		TPZStepSolver<STATE> cg(dohr);
		//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
		
		cg.SetCG(500,pre,5.e-6,0);
		cg.Solve(rhs,diag);

		diag.Print("diag");

        
		TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
		if (!dohrptr) {
			DebugStop(); 
		}
         
		
		dohrptr->AddInternalSolution(diag);
        
        TPZMaterial * mat = cmeshauto->FindMaterial(1);
		int nstate = mat->NStateVariables();
		int nscal = 0, nvec = 0;
		if(nstate ==1) 
		{
			nscal = 1;
		}
		else
		{
			nvec = 1;
		}
		TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
		if(nscal == 1)
		{
			scalnames[0]="state";            
		}
		else
		{
			vecnames[0] = "state";
		}
		
		//cmeshauto->Solution().Print();
		
		std::string postprocessname("ugabuga.vtk");
		TPZVTKGraphMesh vtkmesh(cmesh.operator->(),dim,mat,scalnames,vecnames);
		vtkmesh.SetFileName(postprocessname);
		vtkmesh.SetResolution(0);
		int numcases = 1;
		
		// Iteracoes de tempo
		int istep = 0, nsteps = 80;
		vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(istep, 1.);

		
		typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > subtype;
		const subtype &sublist = dohrptr->SubStructures(); 
		subtype::const_iterator it = sublist.begin();
		int subcount=0;
		while (it != sublist.end()) 
		{
			TPZFMatrix<STATE> subext,subu;
			dohrptr->fAssembly->Extract(subcount,diag,subext);
			(*it)->UGlobal(subext,subu);
			TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
			submesh->LoadSolution(subu);
			subu.Print();
			
		
			std::map<int ,TPZMaterial * > materialmap(submesh->MaterialVec());
			std::map<int ,TPZMaterial * >::iterator itmat;
			for (itmat = materialmap.begin(); itmat != materialmap.end() ; itmat++) 
			{
				TPZMaterial * mat = itmat->second;
				TPZViscoelastic *vmat = dynamic_cast< TPZViscoelastic *> (mat);
				if(vmat)
				{
					vmat->SetUpdateMem();
				}
			}	
			         
			subcount++;
			it++;
		}
		
        /*
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo iterativo",sout);
			LOGPZ_INFO(loggernathan,sout.str())
		}
#endif	
	*/	
		
		//ViscoElastico
		
        vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(istep+1, 1.);
        
		std::cout << "To seguindo!!!" << std::endl;
		for (istep = 2 ; istep < nsteps ; istep++)
		{
			TPZAutoPointer<TPZGuiInterface> guifake;
			dohrstruct.Assemble(rhs, guifake);
			cg.Solve(rhs,diag);	
			
			dohrptr->AddInternalSolution(diag);
			
			// Colocando a solucao na malha
			typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > subtype;
			const subtype &sublist = dohrptr->SubStructures(); 
			subtype::const_iterator it = sublist.begin();
			int subcount=0;
			while (it != sublist.end()) 
			{
				TPZFMatrix<STATE> subext,subu;
				dohrptr->fAssembly->Extract(subcount,diag,subext);
				(*it)->UGlobal(subext,subu);
				TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
				submesh->LoadSolution(subu);
				subcount++;
				it++;
			}

			vtkmesh.DrawMesh(numcases);
			vtkmesh.DrawSolution(istep, 1.);	
		}
	}
	
	total.stop();
	std::cout << "TEMPO = " << total.seconds() << std::endl;
	
	delete gmesh;

	return EXIT_SUCCESS;
}

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	STATE E = 1.e6;
	STATE poisson = 0.3;
	TPZManVector<STATE> force(3,0.);
	force[1] = 20.;
	TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
	TPZMaterial * elastauto(elast);
	TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
	TPZMaterial * bcauto(bc);
	mesh->InsertMaterialObject(elastauto);
	mesh->InsertMaterialObject(bcauto);
}

void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	STATE Ela = 1.e6;
	STATE poisson = 0.2;
	TPZManVector<STATE> force(3,0.);
	force[2] = -20.;
	STATE ElaE = 1000000., poissonE = 0.2, ElaV = 950000., poissonV = 0.14; 
	
	STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alpha = 1.;	
	deltaT = 0.1;
	
	TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
	viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
	
	TPZMaterial * viscoelastauto(viscoelast);
	TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = viscoelast->CreateBC(viscoelastauto, -1, 0, val1, val2);
	TPZFNMatrix<6> qsi(6,1,0.);
	viscoelast->SetDefaultMem(qsi); //elast
	int index = viscoelast->PushMemItem(); //elast
	TPZMaterial * bcauto(bc);
	mesh->InsertMaterialObject(viscoelastauto);
	mesh->InsertMaterialObject(bcauto);	
}

void InsertViscoElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1, neumann = 1, mixed = 2;
	//	int dirichlet = 0;
	int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5, dirp2 = -6;
	TPZManVector<STATE> force(3,0.);
	//force[1] = 0.;

	STATE ElaE = 1000., poissonE = 0.2, ElaV = 970., poissonV = 0.14; 

	STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alpha = 1.;	
	deltaT = 0.1;
	
	TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
	viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
	//TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
	//TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
	
	TPZFNMatrix<6> qsi(6,1,0.);
	viscoelast->SetDefaultMem(qsi); //elast
	int index = viscoelast->PushMemItem(); //elast
	TPZMaterial * viscoelastauto(viscoelast);
	mesh->InsertMaterialObject(viscoelastauto);
	
	// Neumann em x = 1;
	TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
	val2(0,0) = 1.;
	TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
	TPZMaterial * bcauto4(bc4);
	mesh->InsertMaterialObject(bcauto4);
	
	// Neumann em x = -1;
	val2(0,0) = -1.;
	TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
	TPZMaterial * bcauto5(bc5);
	mesh->InsertMaterialObject(bcauto5);
	
	val2.Zero();
	// Dirichlet em -1 -1 -1 xyz;
	val1(0,0) = 1e4;
	val1(1,1) = 1e4;
	val1(2,2) = 1e4;
	TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
	TPZMaterial * bcauto1(bc1);
	mesh->InsertMaterialObject(bcauto1);
	
	// Dirichlet em 1 -1 -1 yz;
	val1(0,0) = 0.;
	val1(1,1) = 1e4;
	val1(2,2) = 1e4;
	TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
	TPZMaterial * bcauto2(bc2);
	mesh->InsertMaterialObject(bcauto2);
	
	// Dirichlet em 1 1 -1 z;
	val1(0,0) = 0.;
	val1(1,1) = 0.;
	val1(2,2) = 1e4;
	TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
	TPZMaterial * bcauto3(bc3);
	mesh->InsertMaterialObject(bcauto3);
	
}

TPZGeoMesh *MalhaPredio()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/Substruct_Viscoelastic/";
	FileName += "8andares02.txt";
	
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
		int matBCid = -1;
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
			
			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		// Colocando as condicoes de contorno
		
		for(el=0; el<numelements; el++)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			// na face z = 0
			TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
				if (nodecoord[2] == 0.)
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
				TPZGeoElBC(tetraSide,matBCid);	
			}
		}
	}
	
	ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	ofstream predio("GeoPredio.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true); 
	
	return gMesh;
}

TPZGeoMesh *MalhaCubo()
{
	int numnodes=-1;
	int numelements=-1;
	
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/Substruct_Viscoelastic/";
	FileName += "cube1.txt";
	
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
			
			int index = el;
			
			TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		
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
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
	}
	
	ofstream arg("malhaPZ1BC.txt");
	gMesh->Print(arg);
	
	std::ofstream out("Cube.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh, out, true);
	
	return gMesh;
}

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

REAL Height(TPZGeoMesh *gmesh)
{
	TPZAdmChunkVector<TPZGeoNode> &nodevec = gmesh->NodeVec();
	int nnodes = nodevec.NElements();
	int in;
	REAL maxz = 0.;
	for (in=0; in<nnodes; in++) {
		REAL z = nodevec[in].Coord(2);
		maxz = (maxz < z) ? z : maxz;
	}
	return maxz;
}

int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height)
{
	int nelem = cmesh->NElements();
	TPZManVector<int> subindex(nelem,-1);
	int iel;
	int nsub = 0;
	for (iel=0; iel<nelem; iel++) 
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		if (!gel) {
			continue;
		}
		int nsides = gel->NSides();
		TPZManVector<REAL> center(gel->Dimension(),0.), xco(3,0.);
		gel->CenterPoint(nsides-1,center);
		gel->X(center,xco);
		REAL z = xco[2];
		int floor = (int) z/height;
		nsub = (floor+1) > nsub ? (floor+1) : nsub;
		subindex[iel] = floor;
	}
	
#ifdef DEBUG 
	{
		TPZGeoMesh *gmesh = cmesh->Reference();
		int nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int cel;
		int nel = cmesh->NElements();
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = cmesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = subindex[cel];
		}
		ofstream vtkfile("partition.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	
	int isub;
	TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
	for (isub=0; isub<nsub; isub++) 
	{
		int index;
		std::cout << '^'; std::cout.flush();
		submeshes[isub] = new TPZSubCompMesh(cmesh,index);
		
		if (index < subindex.NElements()) 
		{
			subindex[index] = -1;
		}
	}
	for (iel=0; iel<nelem; iel++) 
	{
		int domindex = subindex[iel];
		if (domindex >= 0) 
		{
			TPZCompEl *cel = cmesh->ElementVec()[iel];
			if (!cel) 
			{
				continue;
			}
			submeshes[domindex]->TransferElement(cmesh.operator->(),iel);
		}
	}
	cmesh->ComputeNodElCon();
	for (isub=0; isub<nsub; isub++) 
	{
		submeshes[isub]->MakeAllInternal();
		std::cout << '*'; std::cout.flush();
	}
	
	cmesh->ComputeNodElCon();
	cmesh->CleanUpUnconnectedNodes();
	return nsub;
}
