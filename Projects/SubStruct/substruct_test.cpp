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
#endif


void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);

using namespace std;

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	TPZTimer total;
	total.start();
	std::cout << "COMECA O TEMPO"<< std::endl;
	
	int dim = 3;
	int plevel = 1;
	
	TPZPairStructMatrix::gNumThreads = 0; //RAUL
	int numthreads = 2;
	//	tempo.fNumthreads = numthreads;	// alimenta timeTemp com o numero de threads
	TPZGeoMesh *gmesh = 0;
	{
		TPZCompEl::SetgOrder(plevel);
		
		TPZAutoPointer<TPZCompMesh> cmesh;

		gmesh = MalhaCubo();
		cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDimModel(3);
		cmesh->SetDefaultOrder(plevel);
		InsertElasticityCubo(cmesh);
		cmesh->AutoBuild();		
		
		std::cout << "Numero de equacoes " << cmesh->NEquations() << std::endl;
		
        {
            std::ofstream sout("cmesh_ref.txt");
            cmesh->Print(sout);
        }
        
		int numthread_assemble = 0; //RAUL
		int numthread_decompose = 0; //RAUL
		TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
		TPZDohrStructMatrix dohrstruct(cmeshauto);
		
		dohrstruct.IdentifyExternalConnectIndexes();
		
		std::cout << "Substructuring the mesh\n";
		//	TPZfTime timetosub; // init of timer
		//REAL height = Height(gmesh);
		//int nsubstruct = SubStructure(cmesh, height/2);
		
		int nsubstruct = 4; //RAUL
		dohrstruct.SubStructure(nsubstruct);
		//	tempo.ft0sub = timetosub.ReturnTimeDouble();  // end of timer
		//	std::cout << tempo.ft0sub << std::endl;
		
		//	sub.SubStructure();
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
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

		dohrstruct.Assemble(*matptr,rhs, gui, numthread_assemble, numthread_decompose);
		
		TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
		TPZAutoPointer<TPZMatrix<STATE> > precond = dohrstruct.Preconditioner();
	
		int neq = dohr->Rows();
	
		TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
		
		std::cout << "Numero de equacoes " << neq << std::endl;
		
		TPZStepSolver<STATE> pre(precond);
		pre.SetMultiply();
		TPZStepSolver<STATE> cg(dohr);
		//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
		
		cg.SetCG(500,pre,1.e-8,0);
		cg.Solve(rhs,diag);
		
		
		TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
		if (!dohrptr) {
			DebugStop();
		}
		
		dohrptr->AddInternalSolution(diag);
		
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
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			diag.Print("Resultado do processo iterativo",sout);
			LOGPZ_INFO(loggerconverge,sout.str())
		}
#endif	
		
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
		TPZManVector<std::string> scalnames(nscal),vecnames(nvec),tennames(0);
		if(nscal == 1)
		{
			scalnames[0]="state";            
		}
		else
		{
			vecnames[0] = "state";
		}
		std::string postprocessname("dohrmann_elastic.vtk");
        if(!mat){
            DebugStop();
        }
        std::set<int> matids;
        int matid = mat->Id();
        matids.insert(matid);
		TPZVTKGraphMesh vtkmesh(cmesh.operator->(),dim,matids,scalnames,vecnames,tennames);
		vtkmesh.SetFileName(postprocessname);
		vtkmesh.SetResolution(1);
		int numcases = 1;
	
		vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(0, 1.);
	}
	
	total.stop();
	std::cout << "TEMPO = " << total.seconds() << std::endl;
	
	delete gmesh;
	
	return EXIT_SUCCESS; 
}


void InsertElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1, neumann = 1, mixed = 2;
	//	int dirichlet = 0;
	int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
	TPZManVector<STATE> force(3,0.);
	//force[1] = 0.;
    
	STATE ElaE = 1000., poissonE = 0.2;   //, poissonV = 0.1, ElaV = 100.; 
    
	STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alpha = 1.;	
	deltaT = 0.01;
	
	//TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
	//viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
	//TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
	TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
	
	TPZFNMatrix<6> qsi(6,1,0.);
	//viscoelast->SetDefaultMem(qsi); //elast
	//int index = viscoelast->PushMemItem(); //elast
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

TPZGeoMesh *MalhaCubo()
{
	int64_t numnodes=-1;
	int64_t numelements=-1;
	
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/SubStruct/";
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
	
	TPZManVector <int64_t> TopolTetra(4);
	
	const int64_t Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int64_t nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int64_t in;
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
		
		int64_t l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		
		int64_t el;
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
			
			int64_t index = el;
			
			new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		
		// Colocando as condicoes de contorno
		for(el=0; el<numelements; el++)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			
			// na face x = 1
			TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
			for (int i = 0; i < 4; i++) 
			{
				int64_t pos = tetra->NodeIndex(i);
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
				int64_t pos = tetra->NodeIndex(i);
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
	int64_t iel;
	int64_t nelem = gr->ElementVec().NElements();
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
