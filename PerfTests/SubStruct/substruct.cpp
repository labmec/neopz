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
#include <chrono>
#include <thread>
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
#include "TPZParallelUtils.h"
#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"

#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"

#include "pzvtkmesh.h"

#include "pzlog.h"

#include <string>
#include "TPZSimpleTimer.h"


#ifdef PZ_LOG
static TPZLogger loggerconverge("pz.converge");
static TPZLogger logger("main");
#endif

//#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef HAS_GETRUSAGE
#include <sys/resource.h> // getrusage
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh);
void InsertViscoElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh);
TPZGeoMesh *MalhaPredio();
TPZGeoMesh *MalhaCubo();
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
REAL Height(TPZGeoMesh *gmesh);
int SubStructure(TPZAutoPointer<TPZCompMesh> cmesh, REAL height);

#define PERF_START(obj)\
auto start_##obj = std::chrono::system_clock::now();
#define PERF_STOP(obj)\
auto end_##obj = std::chrono::system_clock::now();\
std::chrono::duration<double> elapsed_##obj = end_##obj - start_##obj;
#define PERF_PRINT(obj){\
std::cout << std::string("time_")+std::string(#obj)+std::string(": ") << elapsed_##obj.count() << "s\n";\
}

enum meshChoice {EPredio, ECubo};
constexpr static int nHRef{3};
int main(int argc, char *argv[])
{
    /*OPTIONS*/
    meshChoice msh = EPredio;
    constexpr int dim{3};
    constexpr int initPOrder{2};
    constexpr int nSubMeshes{3};
    unsigned int nThreads{std::thread::hardware_concurrency()};
    constexpr bool useTBB{false};
    constexpr int numIt{500};//CG solver
    /*CODE*/
    TPZSimpleTimer serial;
    
#ifdef USING_LIKWID
    likwid_manager_t likwid_manager;
#endif

#ifdef USING_TBB
    task_scheduler_init init;
#endif
    
    int main_ret = EXIT_SUCCESS;
    
    PERF_START(total_rst);
    
    TPZPairStructMatrix::gNumThreads = nThreads;
    TPZGeoMesh  *gmesh = 0;
    TPZAutoPointer<TPZCompMesh> cmeshauto = 0;
    TPZDohrStructMatrix* dohrstruct = 0;
    TPZFMatrix<STATE> *rhs = NULL;
    TPZMatrix<STATE> *matptr = 0;

    switch(msh)
      {
      case EPredio:
        {
          gmesh = MalhaPredio();
          cmeshauto = new TPZCompMesh(gmesh);
          cmeshauto->SetDefaultOrder(initPOrder);
          cmeshauto->SetDimModel(dim);
          InsertElasticity(cmeshauto);
          cmeshauto->AutoBuild();
          break;
        }
      case ECubo:
        {
          gmesh = MalhaCubo();
          cmeshauto = new TPZCompMesh(gmesh);
          cmeshauto->SetDimModel(dim);
          cmeshauto->SetDefaultOrder(initPOrder);
          //cmeshauto->SetAllCreateFunctionsContinuousWithMem();
          //cmeshauto->SetAllCreateFunctionsContinuous();
          InsertViscoElasticityCubo(cmeshauto);
          cmeshauto->AutoBuild();
          break;
        }
      }
   
    dohrstruct = new TPZDohrStructMatrix(cmeshauto);
    dohrstruct->IdentifyExternalConnectIndexes();
    dohrstruct->SubStructure(nSubMeshes);
        
#ifdef PZ_LOG
    {
      std::stringstream str;
      cmeshauto->Print(str);
      if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,str.str());
    }
#endif
      
    PERF_START(create_rst);
    matptr = dohrstruct->Create();
    PERF_STOP(create_rst);
    PERF_PRINT(create_rst);
    
    TPZAutoPointer<TPZMatrix<STATE> > precond = NULL;
    
    dohrstruct->SetNumThreads(nThreads);
        
    PERF_START(assemble_rst);
    TPZAutoPointer<TPZGuiInterface> gui;
    rhs = new TPZFMatrix<STATE>(cmeshauto->NEquations(),1,0.);
    if(useTBB)
      dohrstruct->AssembleTBB(*matptr,*rhs, gui);
    else
      dohrstruct->Assemble(*matptr,*rhs, gui, nThreads, nThreads);
    PERF_STOP(assemble_rst);    
    PERF_PRINT(assemble_rst);
    PERF_START(precond_rst);
    precond = dohrstruct->Preconditioner();
    PERF_STOP(precond_rst);
    PERF_PRINT(precond_rst);
    /* Work after checkpoint 3 */
    TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
        
    int neq = dohr->Rows();
    TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
               
    TPZStepSolver<STATE> pre(precond);
    pre.SetMultiply();
    TPZStepSolver<STATE> cg(dohr);
        
    /* Configure the CG solver to iterate:
       - until it converges (residual <= 1.e-8), or
       - until it reaches 500 itearations.
    */
    cg.SetCG(numIt,pre,1.e-8,0);
        
    PERF_START(solve_rst);
    cg.Solve(*rhs,diag);
    PERF_STOP(solve_rst);
    PERF_PRINT(solve_rst);
    /* checking if the solver converged */
    if (cg.GetTolerance() > 1.e-8)
      {
        std::cerr << "ERROR: solver do not converged with the limit of iterations."  << std::endl;
        exit(1);
      }
        
    TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr =
      dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
        
    if (!dohrptr) {
      DebugStop();
    }
        
    dohrptr->AddInternalSolution(diag);

        
        
       
    typedef std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > subtype;
    const subtype &sublist = dohrptr->SubStructures();
    subtype::const_iterator it = sublist.begin();
    int subcount=0;
    while (it != sublist.end()) {
            
      TPZFMatrix<STATE> subext,subu;
      dohrptr->fAssembly->Extract(subcount,diag,subext);
      (*it)->UGlobal(subext,subu);
      TPZCompMesh *submesh = SubMesh(cmeshauto, subcount);
      submesh->LoadSolution(subu);
            
      //ViscoElastico
      //Atualizando memoria do material
      std::map<int ,TPZMaterial * > materialmap(submesh->MaterialVec());
      std::map<int ,TPZMaterial * >::iterator itmat;
      for (itmat = materialmap.begin(); itmat != materialmap.end() ; itmat++)
        {
          TPZMaterial * mat = itmat->second;
          TPZViscoelastic *vmat = dynamic_cast< TPZViscoelastic *> (mat);
          if(vmat)
            {
              vmat->SetUpdateMem(true);
            }
        }
      subcount++;
      it++;
    }
        
#ifdef PZ_LOG
    {
      std::stringstream sout;
      diag.Print("Resultado do processo iterativo",sout);
      LOGPZ_INFO(loggerconverge,sout.str())
        }
#endif
        
    TPZMaterial * mat = cmeshauto->FindMaterial(1);
    std::set<int> matIds;
    matIds.insert(1);
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
    TPZManVector<std::string> scalnames(nscal),vecnames(nvec), tensnames(0);
    if(nscal == 1)
      {
        scalnames[0]="state";
      }
    else
      {
        vecnames[0] = "state";
      }
    std::string postprocessname("dohrmann_visco.vtk");
    TPZVTKGraphMesh vtkmesh(cmeshauto.operator->(),dim,matIds,scalnames,vecnames, tensnames);
    vtkmesh.SetFileName(postprocessname);
    vtkmesh.SetResolution(1);
    int numcases = 1;
        
        
    // Iteracoes de tempo
    int istep = 0;
    vtkmesh.DrawMesh(numcases);
    vtkmesh.DrawSolution(istep, 1.);
    
    PERF_STOP(total_rst);
    

    PERF_PRINT(create_rst);
    PERF_PRINT(assemble_rst);
    PERF_PRINT(precond_rst);    
    PERF_PRINT(solve_rst);    
    PERF_PRINT(total_rst);
    
    if (gmesh != NULL) delete gmesh;
    
    return main_ret;
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
    TPZManVector<STATE> force(3,0.);
    force[1] = 20.;
    STATE ElaE = 1000., poissonE = 0.2, ElaV = 100., poissonV = 0.1;
    
    STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alpha = 1.;
    deltaT = 0.01;
    
    TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
    viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
    
    TPZMaterial * viscoelastauto(viscoelast);
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    TPZBndCond *bc = viscoelast->CreateBC(viscoelastauto, -1, 0, val1, val2);
    TPZFNMatrix<6,STATE> qsi(6,1,0.);
    viscoelast->SetDefaultMem(qsi); //elast
    viscoelast->PushMemItem(); //elast
    TPZMaterial * bcauto(bc);
    mesh->InsertMaterialObject(viscoelastauto);
    mesh->InsertMaterialObject(bcauto);
}

void InsertViscoElasticityCubo(TPZAutoPointer<TPZCompMesh> mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //      int dirichlet = 0;
    int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;
    TPZManVector<STATE> force(3,0.);
    //force[1] = 0.;
    REAL Ela = 1000, poisson = 0.;
    REAL lambdaV = 0, muV = 0, alphaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alphaT = 0.01;
    
    
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, Ela, poisson, lambdaV, muV, alphaT, force);
    TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, Ela, poisson, force);
    
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

TPZGeoMesh *MalhaPredio()
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
#ifndef MACOSX
	std::string FileName = "8andares02.txt";
#else
  std::string FileName = "../8andares02.txt";
#endif
	{
		bool countnodes = false;
		bool countelements = false;
		
		std::ifstream read (FileName.c_str());
		if (!read.is_open()) {
      std::cerr << "Could not open file: " << FileName << std::endl;
      exit(1);
		}
        
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
	
	TPZVec <int64_t> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int64_t nodeId = 0, elementId = 0, matElId = 1;
	
	std::ifstream read;
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
			
			int64_t index = el;
			
			new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		// Colocando as condicoes de contorno
		
		pzutils::ParallelFor(0,numelements,[&](int el)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			// na face z = 0
			TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
			for (int i = 0; i < 4; i++)
			{
				int64_t pos = tetra->NodeIndex(i);
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
		});
	}
	
	std::ofstream arg("malhaPZ.txt");
	gMesh->Print(arg);
	std::ofstream predio("GeoPredio.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true);
	
	return gMesh;
}

TPZGeoMesh *MalhaCubo()
{
	int numnodes=-1;
	int numelements=-1;
	
	std::string FileName = "cube1.txt";
	{
		bool countnodes = false;
		bool countelements = false;

		std::ifstream read (FileName.c_str());
		if (!read.is_open()) {
            std::cerr << "Could not open file: " << FileName << std::endl;
            exit(1);
		}
		
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
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int nodeId = 0, elementId = 0, matElId = 1;
	
	std::ifstream read;
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
			
			int64_t index = el;
			
			new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		
		// Colocando as condicoes de contorno
		pzutils::ParallelFor(0,numelements,[&](int el)
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
			
		});
		
		TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
		yz[0] = 1.;
		z[2] = -1;
		int bcidxyz = -1, bcidyz = -2, bcidz = -3;
		SetPointBC(gMesh, xyz, bcidxyz);
		SetPointBC(gMesh, yz, bcidyz);
		SetPointBC(gMesh, z, bcidz);
		
	}
	
	/* refine mesh */
  for ( int ref = 0; ref < nHRef; ref++ ){
    TPZVec<TPZGeoEl *> filhos;
    int64_t n = gMesh->NElements();
    for ( int64_t i = 0; i < n; i++ ){
      TPZGeoEl * gel = gMesh->ElementVec() [i];
		        
      if(!gel) continue;
      if(gel->Dimension() < 1) continue;
      if(gel->HasSubElement()) continue;
                
      gel->Divide (filhos);
    }
  }
	
	std::ofstream arg("malhaPZ1BC.txt");
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
	int64_t nelem = cmesh->NElements();
	TPZManVector<int64_t> subindex(nelem,-1);
	int64_t iel;
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
	
#ifdef PZDEBUG
	{
		TPZGeoMesh *gmesh = cmesh->Reference();
		int64_t nelgeo = gmesh->NElements();
		TPZVec<int> domaincolor(nelgeo,-999);
		int64_t cel;
		int64_t nel = cmesh->NElements();
		for (cel=0; cel<nel; cel++) {
			TPZCompEl *compel = cmesh->ElementVec()[cel];
			if(!compel) continue;
			TPZGeoEl *gel = compel->Reference();
			if (!gel) {
				continue;
			}
			domaincolor[gel->Index()] = subindex[cel];
		}
		std::ofstream vtkfile("partition.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, domaincolor);
	}
#endif
	
	int isub;
	TPZManVector<TPZSubCompMesh *> submeshes(nsub,0);
	for (isub=0; isub<nsub; isub++)
	{
		int64_t index;
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
