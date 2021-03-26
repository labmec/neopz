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

#include "pzfstrmatrix.h"

#ifdef PZ_LOG
static TPZLogger loggerconverge("pz.converge");
static TPZLogger logger("main");
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
// nao funcione
int main2(int argc, char *argv[]);

// read starting from checkpoint1
int main3(int argc, char *argv[]);

// read starting from checkpoint2
int main4(int argc, char *argv[]);

// read starting from checkpoint3
int main5(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    return main1(argc, argv);
}

int main1(int argc, char *argv[])
{
	
	TPZTimer total;
	total.start();
	std::cout << "COMECA O TEMPO"<< std::endl;
	
	int dim = 2;
	int maxlevel = 5;
	int sublevel = 3;
	int plevel = 1;
	TPZPairStructMatrix::gNumThreads = 4;
	int numthreads = 4;
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
			if (1) // Predio Elastisco
			{
				gmesh = MalhaPredio();
				cmesh = new TPZCompMesh(gmesh);
				cmesh->SetDimModel(3);
				InsertElasticity(cmesh);
				cmesh->AutoBuild();
				
			}
			else // Cubo Viscoso
			{
				gmesh = MalhaCubo();
				cmesh = new TPZCompMesh(gmesh);
				cmesh->SetDimModel(3);
				cmesh->SetDefaultOrder(plevel);
				cmesh->SetAllCreateFunctionsContinuous();
                InsertElasticity(cmesh);
//				cmesh->SetAllCreateFunctionsContinuousWithMem();
//				InsertViscoElasticityCubo(cmesh);
				cmesh->AutoBuild();
			}
		}
/*		
        {
            std::ofstream sout("cmesh_nw.txt");
            cmesh->Print(sout);
        }
*/
		std::cout << "Numero de equacoes " << cmesh->NEquations() << std::endl;
		
		int numthread_assemble = 4;
		int numthread_decompose = 4;
		TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
        //
        TPZFStructMatrix fullstruct(cmeshauto);
        fullstruct.SetNumThreads(numthread_assemble);
        int64_t sz = cmeshauto->NEquations();
        TPZFMatrix<STATE> rhs_t(sz, 1);
        fullstruct.Assemble(rhs_t, 0);
        return 0;
        
		TPZDohrStructMatrix dohrstruct(cmeshauto);
		
		dohrstruct.IdentifyExternalConnectIndexes();
		
		std::cout << "Substructuring the mesh\n";
		//	TPZfTime timetosub; // init of timer
		//REAL height = Height(gmesh);
		//int nsubstruct = SubStructure(cmesh, height/2);
		
		dohrstruct.SubStructure(8);
		//	tempo.ft0sub = timetosub.ReturnTimeDouble();  // end of timer
		//	std::cout << tempo.ft0sub << std::endl;
		
		//	sub.SubStructure();
#ifdef PZ_LOG
		if(logger.isDebugEnabled())
		{
			std::stringstream str;
			cmesh->Print(str);
			LOGPZ_DEBUG(logger,str.str());
		}
#endif
/*
        {
            TPZFileStream CheckPoint1;
            CheckPoint1.OpenWrite("CheckPoint1.txt");
            cmesh->Reference()->Write(CheckPoint1, 0);
            cmesh->Write(CheckPoint1, 0);
            dohrstruct.Write(CheckPoint1);
        }
        {
            TPZFileStream CheckPoint1;
            CheckPoint1.OpenRead("CheckPoint1.txt");
            TPZGeoMesh locgmesh;
            locgmesh.Read(CheckPoint1, 0);
            TPZCompMesh *loccmesh = new TPZCompMesh(&locgmesh);
            TPZAutoPointer<TPZCompMesh> loccmeshauto(loccmesh);
            loccmesh->Read(CheckPoint1, &locgmesh);
            TPZDohrStructMatrix locdohrstruct(loccmeshauto);
            locdohrstruct.Read(CheckPoint1);
            
        }
*/		
		dohrstruct.SetNumThreads(numthreads);
		
		TPZAutoPointer<TPZGuiInterface> gui;
		TPZFMatrix<STATE> rhs(cmesh->NEquations(),1,0.);
        
        TPZMatrix<STATE> *matptr = dohrstruct.Create();
/*				
        {
            TPZFileStream CheckPoint2;
            CheckPoint2.OpenWrite("CheckPoint2.txt");
            cmesh->Reference()->Write(CheckPoint2, 0);
            cmesh->Write(CheckPoint2, 0);
            matptr->Write(CheckPoint2, 1);
            dohrstruct.Write(CheckPoint2);
        }
        {
            TPZFileStream CheckPoint2;
            CheckPoint2.OpenRead("CheckPoint2.txt");
            TPZGeoMesh gmesh;
            gmesh.Read(CheckPoint2,0);
            TPZCompMesh *cmesh = new TPZCompMesh;
            cmesh->Read(CheckPoint2, &gmesh);
            TPZAutoPointer<TPZCompMesh> loccmeshauto(cmesh);
            TPZMatrix<STATE> *mat;
            mat = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint2, 0));
            delete mat;
            TPZDohrStructMatrix locdohrstruct(loccmeshauto);
            locdohrstruct.Read(CheckPoint2);
            
        }
*/        
        std::cout << "Assembling the substructures\n";
        dohrstruct.Assemble(*matptr,rhs, gui, numthread_assemble, numthread_decompose);

		TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
		TPZAutoPointer<TPZMatrix<STATE> > precond = dohrstruct.Preconditioner();
/*		
        {
            TPZFileStream CheckPoint3;
            CheckPoint3.OpenWrite("CheckPoint3.txt");
            cmesh->Reference()->Write(CheckPoint3, 0);
            cmesh->Write(CheckPoint3, 0);
            dohr->Write(CheckPoint3, 1);
            precond->Write(CheckPoint3, 1);
            rhs.Write(CheckPoint3, 0);
        }
        {
            TPZFileStream CheckPoint3;
            CheckPoint3.OpenRead("CheckPoint3.txt");
            TPZGeoMesh gmesh;
            gmesh.Read(CheckPoint3, 0);
            TPZCompMesh cmesh;
            cmesh.Read(CheckPoint3, &gmesh);
            TPZMatrix<STATE> *matdohr;
            matdohr = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint3, 0));
            TPZMatrix<STATE> *matprecond;
            matprecond = dynamic_cast<TPZMatrix<STATE> *>(TPZSavable::Restore(CheckPoint3, matdohr));
            delete matprecond;
            delete matdohr;
            TPZFMatrix<STATE> rhsloc;
            rhsloc.Read(CheckPoint3, 0);
        }
*/        
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
					vmat->SetUpdateMem();
				}
			}	
			subcount++;
			it++;
		}
		
#ifdef PZ_LOG
		if(loggerconverge.isInfoEnabled())
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
		TPZManVector<std::string> scalnames(nscal),vecnames(nvec), tennames;
		if(nscal == 1)
		{
			scalnames[0]="state";            
		}
		else
		{
			vecnames[0] = "state";
		}
		std::string postprocessname("dohrmann_visco.vtk");
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
		
		
		// Iteracoes de tempo
		int istep = 0;
		vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(istep, 1.);
	}
	
	total.stop();
	std::cout << "TEMPO = " << total.seconds() << std::endl;
	
	delete gmesh;

	return EXIT_SUCCESS;
}

int main3(int argc, char *argv[])
{
    int numthreads = 0;
    int dim = 2;

    TPZGeoMesh *gmesh = NULL;
    {
        int numthread_assemble = 0;
        int numthread_decompose = 0;

        TPZCompMesh *loccmesh = new TPZCompMesh(gmesh);
        TPZAutoPointer<TPZCompMesh> cmesh(loccmesh);
        TPZAutoPointer<TPZDohrStructMatrix> dohrstruct = new TPZDohrStructMatrix(cmesh.operator ->());
        {
            TPZPersistenceManager::OpenRead("CheckPoint1.txt");
            gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
            cmesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
            dohrstruct = dynamic_cast<TPZDohrStructMatrix *>(TPZPersistenceManager::ReadFromFile());
            TPZPersistenceManager::CloseRead();
        }
        {
            TPZPersistenceManager::OpenWrite("CheckPoint6.txt");
            TPZPersistenceManager::WriteToFile(gmesh);
            TPZPersistenceManager::WriteToFile(cmesh.operator ->());
            TPZPersistenceManager::WriteToFile(dohrstruct.operator ->());
            TPZPersistenceManager::CloseWrite();
            
        }
        dim = cmesh->Dimension();
        
        dohrstruct->SetNumThreads(numthreads);
        
        TPZAutoPointer<TPZGuiInterface> gui;
        TPZFMatrix<STATE> rhs(cmesh->NEquations(),1,0.);
        
        TPZMatrix<STATE> *matptr = dohrstruct->Create();
        {
            TPZPersistenceManager::OpenWrite("CheckPoint4.txt");
            TPZPersistenceManager::WriteToFile(cmesh->Reference());
            TPZPersistenceManager::WriteToFile(cmesh.operator ->());
            TPZPersistenceManager::WriteToFile(matptr);
            TPZPersistenceManager::WriteToFile(dohrstruct.operator ->());
            TPZPersistenceManager::CloseWrite();
        }
        {
            TPZPersistenceManager::OpenRead("CheckPoint4.txt");
            TPZGeoMesh *gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
            TPZCompMesh *cmesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
            TPZAutoPointer<TPZCompMesh> loccmeshauto(cmesh);
            TPZMatrix<STATE> *mat = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            delete mat;
            TPZDohrStructMatrix *locdohrstruct = dynamic_cast<TPZDohrStructMatrix*>(TPZPersistenceManager::ReadFromFile());
            locdohrstruct->SetMesh(loccmeshauto);
            delete gmesh;
        }
        dohrstruct->Assemble(*matptr,rhs, gui, numthread_assemble, numthread_decompose);
        
        TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
        TPZAutoPointer<TPZMatrix<STATE> > precond = dohrstruct->Preconditioner();
        
        {
            TPZPersistenceManager::OpenWrite("CheckPoint5.txt");
            TPZPersistenceManager::WriteToFile(cmesh->Reference());
            TPZPersistenceManager::WriteToFile(cmesh.operator ->());
            TPZPersistenceManager::WriteToFile(dohr.operator ->());
            TPZPersistenceManager::WriteToFile(precond.operator ->());
            TPZPersistenceManager::WriteToFile(&rhs);
            TPZPersistenceManager::CloseWrite();
        }
        {
            TPZPersistenceManager::OpenRead("CheckPoint5.txt");
            TPZGeoMesh *gmesh = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::ReadFromFile());
            TPZCompMesh *cmesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matdohr = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matprecond = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            TPZFMatrix<STATE> *rhsloc = dynamic_cast<TPZFMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            delete rhsloc;
            delete matprecond;
            delete matdohr;
            delete cmesh;
            delete gmesh;
        }
        
        int neq = dohr->Rows();
        
        
        TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
        TPZStepSolver<STATE> pre(precond);
        pre.SetMultiply();
        TPZStepSolver<STATE> cg(dohr);
        //  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
        
        cg.SetCG(500,pre,1.e-8,0);
        cg.Solve(rhs,diag);

        
        std::cout << "Numero de equacoes " << neq << std::endl;
        
        TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
        if (!dohrptr) {
            DebugStop();
        }
        
    #ifdef PZ_LOG
        {
            std::stringstream sout;
            diag.Print("Resultado do processo antes do ajuste",sout);
            LOGPZ_INFO(loggerconverge,sout.str())
        }
    #endif
        
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
            TPZCompMesh *submesh = SubMesh(cmesh, subcount);
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
                    vmat->SetUpdateMem();
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
        
        TPZMaterial * mat = cmesh->FindMaterial(1);
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
        std::string postprocessname("dohrmann_visco.vtk");
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
        
        
        // Iteracoes de tempo
        int istep = 0;
        vtkmesh.DrawMesh(numcases);
        vtkmesh.DrawSolution(istep, 1.);

    }
    delete gmesh;

    return EXIT_SUCCESS;
    
}

int main4(int argc, char *argv[])
{
    int numthreads = 0;
    int dim = 2;
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    {
        TPZCompMesh *loccmesh = new TPZCompMesh(gmesh);
        TPZAutoPointer<TPZCompMesh> cmesh(loccmesh);
        int numthread_assemble = 0;
        int numthread_decompose = 0;
        TPZDohrStructMatrix *dohrstruct;
        
        dim = cmesh->Dimension();
        
        
        TPZAutoPointer<TPZGuiInterface> gui;
        
        TPZMatrix<STATE> *matptr;
        {
            
            TPZPersistenceManager::OpenRead("CheckPoint2.txt");
            gmesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::ReadFromFile());
            cmesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
            matptr = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            dohrstruct = dynamic_cast<TPZDohrStructMatrix *>(TPZPersistenceManager::ReadFromFile());
            dohrstruct->SetMesh(cmesh);
            dohrstruct->SetNumThreads(numthreads);

            
        }
        TPZFMatrix<STATE> rhs(cmesh->NEquations(),1,0.);
        dohrstruct->Assemble(*matptr,rhs, gui, numthread_assemble, numthread_decompose);
        
        TPZAutoPointer<TPZMatrix<STATE> > dohr = matptr;
        TPZAutoPointer<TPZMatrix<STATE> > precond = dohrstruct->Preconditioner();
//        {
//            TPZFileStream CheckPoint3;
//            CheckPoint3.OpenWrite("CheckPoint5.txt");
//            cmesh->Reference()->Write(CheckPoint3, 0);
//            cmesh->Write(CheckPoint3, 0);
//            dohr->Write(CheckPoint3, 1);
//            precond->Write(CheckPoint3, 1);
//            rhs.Write(CheckPoint3, 0);
//        }

        {
            TPZPersistenceManager::OpenWrite("CheckPoint6.txt");
            TPZPersistenceManager::WriteToFile(cmesh->Reference());
            TPZPersistenceManager::WriteToFile(cmesh.operator ->());
            TPZPersistenceManager::WriteToFile(dohr.operator ->());
            TPZPersistenceManager::WriteToFile(precond.operator ->());
            TPZPersistenceManager::WriteToFile(&rhs);
            TPZPersistenceManager::CloseWrite();
        }
        {
            TPZPersistenceManager::OpenRead("CheckPoint6.txt");
            TPZAutoPointer<TPZGeoMesh> gmesh = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::ReadFromFile());
            TPZAutoPointer<TPZCompMesh> cmesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matdohr = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matprecond = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            TPZAutoPointer<TPZFMatrix<STATE>> rhsloc = dynamic_cast<TPZFMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            delete matprecond;
            delete matdohr;
        }
        
        int neq = dohr->Rows();
        
        
        
        TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
        TPZStepSolver<STATE> pre(precond);
        pre.SetMultiply();
        TPZStepSolver<STATE> cg(dohr);
        //  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
        
        cg.SetCG(500,pre,1.e-8,0);
        cg.Solve(rhs,diag);
        
        std::cout << "Numero de equacoes " << neq << std::endl;
        
        TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
        if (!dohrptr) {
            DebugStop();
        }
        
#ifdef PZ_LOG
        {
            std::stringstream sout;
            diag.Print("Resultado do processo antes do ajuste",sout);
            LOGPZ_INFO(loggerconverge,sout.str())
        }
#endif
        
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
            TPZCompMesh *submesh = SubMesh(cmesh, subcount);
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
                    vmat->SetUpdateMem();
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
        
        TPZMaterial * mat = cmesh->FindMaterial(1);
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
        std::string postprocessname("dohrmann_visco.vtk");
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
        
        
        // Iteracoes de tempo
        int istep = 0;
        vtkmesh.DrawMesh(numcases);
        vtkmesh.DrawSolution(istep, 1.);
        
        delete dohrstruct;
    }
    delete gmesh;
    
    return EXIT_SUCCESS;
    
}

int main5(int argc, char *argv[])
{
    int numthreads = 0;
    int dim = 2;
    
    TPZAutoPointer<TPZGeoMesh> gmesh;
    {
        TPZAutoPointer<TPZGuiInterface> gui;
        TPZAutoPointer<TPZFMatrix<STATE>> rhs;
        
        TPZAutoPointer<TPZMatrix<STATE> > dohr;
        TPZAutoPointer<TPZMatrix<STATE> > precond;
        
        TPZAutoPointer<TPZCompMesh> cmesh;
        {
            TPZPersistenceManager::OpenRead("CheckPoint3.txt");
            gmesh = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::ReadFromFile());
            cmesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matdohr = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            TPZMatrix<STATE> *matprecond = dynamic_cast<TPZMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
            dohr = matdohr;
            precond = matprecond;
            rhs = dynamic_cast<TPZFMatrix<STATE> *>(TPZPersistenceManager::ReadFromFile());
        }
        
        cmesh->SetReference(gmesh.operator ->());
//        int numthread_assemble = 0;
  //      int numthread_decompose = 0;
        TPZDohrStructMatrix dohrstruct(cmesh);
        dohrstruct.SetNumThreads(numthreads);
        
        dim = cmesh->Dimension();
        
        {
            TPZPersistenceManager::OpenWrite("CheckPoint7.txt");
            TPZPersistenceManager::WriteToFile(cmesh->Reference());
            TPZPersistenceManager::WriteToFile(cmesh.operator ->());
            TPZPersistenceManager::WriteToFile(dohr.operator ->());
            TPZPersistenceManager::WriteToFile(precond.operator ->());
            TPZPersistenceManager::WriteToFile(rhs.operator ->());
            TPZPersistenceManager::CloseWrite();
        }
        
        int neq = dohr->Rows();
        
        
        
        TPZFMatrix<STATE> diag(neq,1,0.), produto(neq,1);
        
        TPZStepSolver<STATE> pre(precond);
        pre.SetMultiply();
        TPZStepSolver<STATE> cg(dohr);
        //  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
        
        cg.SetCG(500,pre,1.e-8,0);
        cg.Solve(*(rhs.operator ->()),diag);
        
        
        std::cout << "Numero de equacoes " << neq << std::endl;
        
        TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr = dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohr.operator->());
        if (!dohrptr) {
            DebugStop();
        }
        
#ifdef PZ_LOG
        {
            std::stringstream sout;
            diag.Print("Resultado do processo antes do ajuste",sout);
            LOGPZ_INFO(loggerconverge,sout.str())
        }
#endif
        
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
            TPZCompMesh *submesh = SubMesh(cmesh, subcount);
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
                    vmat->SetUpdateMem();
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
        
        TPZMaterial * mat = cmesh->FindMaterial(1);
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
        std::string postprocessname("dohrmann_visco.vtk");
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
        
        
        // Iteracoes de tempo
        int istep = 0;
        vtkmesh.DrawMesh(numcases);
        vtkmesh.DrawSolution(istep, 1.);
        
    }
    
    return EXIT_SUCCESS;
    
}

int main2(int argc, char *argv[])
{
	
	/*
	 TPZFMatrix<STATE> teste(2,2);
	 TPZFMatrix<STATE> parte;
	 teste(0,0)=1;
	 teste(0,1)=2;
	 teste(1,0)=3;
	 teste(1,1)=4;
	 teste.GetSub(0,0,2,1,parte);
	 cout << parte << endl;
	 cout << "Hello, world!" << endl;
	 */
	
	
	/**
	 TPZDohrSubstruct meuobjeto;
	 TPZDohrMatrix *matriz = new TPZDohrMatrix();
	 TPZDohrPrecond *precond = new TPZDohrPrecond();
	 TPZStepSolver dohrprecond(precond);
	 dohrprecond.SetMultiply();
	 TPZStepSolver cg(matriz);
	 cg.SetCG(10,dohrprecond,1.e-7,1);
	 TPZFMatrix<STATE> rhs,result;
	 cg.Solve(rhs,result);*/
	//meuobjeto.
	/*  int dim = 2;
	 TPZGenSubStruct sub(dim,6,3);*/
	//	int dim = 2;
	//	int maxlevel = 6;
	//	int sublevel = 3;
	//	int plevel = 3;
	int dim = 2;
	int maxlevel = 4;
	int sublevel = 1;
	int plevel = 2;
	TPZGenSubStruct sub(dim,maxlevel,sublevel);
//	int nk = 8;
	//	int ik;
	//	for(ik=1; ik<nk; ik++)
	//	{
	//	sub.fK[ik] = 1.;//50.*ik;
	//	}
	//sub.fMatDist = TPZGenSubStruct::RandomMat;
	
	TPZCompEl::SetgOrder(plevel);
	
	sub.GenerateMesh();
	
	/*
	 TPZAutoPointer<TPZDohrAssembly> dohrassembly = new TPZDohrAssembly;
	 TPZDohrMatrix<TPZDohrSubstruct> *dohrptr = new TPZDohrMatrix<TPZDohrSubstruct>(dohrassembly);
	 TPZAutoPointer<TPZMatrix> dohr(dohrptr);
	 sub.InitializeDohr(dohr,dohrassembly);
	 
	 
	 // loop over the substructures
	 // This is a lengthy process which should run on the remote processor
	 //	void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct> substruct,  TPZDohrAssembly &dohrassembly);
	 
	 
	 dohrptr->Initialize();
	 #ifdef PZ_LOG
	 {
	 std::stringstream sout;
	 dohrptr->Print("DohrMatrix without condensation", sout);
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif
	 TPZDohrPrecond<TPZDohrSubstruct> *precondptr = new TPZDohrPrecond<TPZDohrSubstruct>(*dohrptr,dohrassembly);
	 precondptr->Initialize();
	 TPZAutoPointer<TPZMatrix> precond(precondptr);
	 
	 */
	
	TPZAutoPointer<TPZDohrAssembly<STATE> > dohrassembly2 = new TPZDohrAssembly<STATE>;
	TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohrptr2 = new TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> >(dohrassembly2);
	dohrptr2->SetNumThreads(4);
	TPZAutoPointer<TPZMatrix<STATE> > dohr2(dohrptr2);
	sub.InitializeDohrCondense(dohr2,dohrassembly2);
	dohrptr2->Initialize();
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		dohr2->Print("The dohr matrix condensed",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
#ifdef PZ_LOG
	if(loggerconverge.isDebugEnabled())
	{
			std::stringstream sout;
			sout << "Three dimensional substructures, maxlevel " << maxlevel << " level of substructures " << sublevel << std::endl;
			sout << "Number of substructures " << dohrptr2->SubStructures().size() << std::endl;
			sout << "Interpolation order " << plevel;
			LOGPZ_DEBUG(loggerconverge,sout.str());
	}
#endif
	
	
	TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> > *precondptr2 = new TPZDohrPrecond<STATE,TPZDohrSubstructCondense<STATE> >(*dohrptr2,dohrassembly2);
	precondptr2->Initialize();
	
	
#ifdef PZ_LOG
	if(loggerconverge.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Printing after creating the preconditioner\n";
		dohrptr2->Print("After creating the preconditioner",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str());
	}
#endif
	
	
	
	TPZAutoPointer<TPZMatrix<STATE> > precond2(precondptr2);
	
	
	TPZFMatrix<STATE> diag(dohr2->Rows(),1,5.), produto(dohr2->Rows(),1), produto2(dohr2->Rows(),1);
	precondptr2->Multiply(diag,produto);
#ifdef PZ_LOG
	if(loggerconverge.isDebugEnabled())
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	
	precondptr2->Multiply(diag,produto2);
#ifdef PZ_LOG
	if(loggerconverge.isDebugEnabled())
	{
		std::stringstream sout;
		produto2.Print("O valor do produto2", sout );
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	//#define TOTAL
#ifdef TOTAL
	{
		int dim=dohr->Rows();
		
		TPZFMatrix<STATE> teste1(dim,dim);
		teste1.Zero();
		TPZFMatrix<STATE> teste2(dim,dim);
		teste2.Zero();
		
		int i,j;
		TPZFMatrix<STATE> col(dim,1); //column of the identity matrix
		TPZFMatrix<STATE> resul(dohr->Rows(),1);
		for (i=0;i<dim;i++) {
			col.Zero();
			col(i,0)=1;
			precondptr->MultAdd(col,col,resul,1,0,0,1);
			for (j=0;j<dim;j++) {
				teste1(i,j) = resul(j,0);
				teste2(i,j) = resul(j,0);
			}
		}
		teste1.Transpose();
		teste1 -= teste2;
		std::cout << "Norma da diferenca das matrizes " << Norm(teste1) << std::endl;
	}
#endif
	std::cout << "Numero de equacoes " << dohr2->Rows() << std::endl;
	//  produto.Print("The value of the product is");
#ifndef MAKEINTERNAL
	diag(0,0) = 0.;
#endif
	dohr2->Multiply(diag,produto);
	dohrptr2->AdjustResidual(produto);
	
#ifdef PZ_LOG
	if(loggerconverge.isDebugEnabled())
	{
		std::stringstream sout;
		produto.Print("O valor do produto", sout );
		diag.Print("O valor da diagonal",sout);
		LOGPZ_DEBUG(loggerconverge,sout.str())
	}
#endif
	diag.Zero();
	TPZStepSolver<STATE> pre(precond2);
	pre.SetMultiply();
	TPZStepSolver<STATE> cg(dohr2);
	//  void SetCG(const int numiterations,const TPZMatrixSolver &pre,const STATE tol,const int FromCurrent);
	
	cg.SetCG(500,pre,1.e-8,0);
	cg.Solve(produto,diag);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		diag.Print("Resultado do processo antes do ajuste",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	
	dohrptr2->AddInternalSolution(diag);
#ifdef PZ_LOG
	{
		std::stringstream sout;
		diag.Print("Resultado do processo iterativo",sout);
		LOGPZ_INFO(loggerconverge,sout.str())
	}
#endif
	//diag.Print("Resultado do solve");
	/* Solve
	 TPZFMatrix<STATE> *teste = new TPZFMatrix(2,2);
	 (*teste)(0,0)=1;
	 (*teste)(0,1)=2;
	 (*teste)(1,0)=3;
	 (*teste)(1,1)=4;
	 TPZStepSolver coef;
	 coef.SetMatrix(teste);
	 coef.SetDirect(ELU);
	 TPZFMatrix<STATE> resul(2,2);
	 resul(0,0)=2;
	 resul(0,1)=3;
	 resul(1,0)=4;
	 resul(1,1)=5;
	 TPZFMatrix<STATE> res(2,2);
	 coef.Solve(resul,res);
	 cout << res << endl;*/
	
	
	return EXIT_SUCCESS;
}

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh)
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

void InsertViscoElasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
//	STATE Ela = 1.e6;
//	STATE poisson = 0.2;
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
	//	int dirichlet = 0;
	int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
	TPZManVector<STATE> force(3,0.);
	//force[1] = 0.;
    
	STATE ElaE = 1000., poissonE = 0.2, ElaV = 100., poissonV = 0.1; 
    
	STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
	lambdaV = 11.3636;
	muV = 45.4545;
	alpha = 1.;	
	deltaT = 0.01;
	
	TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
	viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
	//TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
	//TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
	
	TPZFNMatrix<6,STATE> qsi(6,1,0.);
	viscoelast->SetDefaultMem(qsi); //elast
	viscoelast->PushMemItem(); //elast
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
	int64_t numnodes=-1;
	int64_t numelements=-1;
	
	string FileName, dirname = PZSOURCEDIR;
	FileName = dirname + "/Projects/SubStruct/";
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
	
	TPZVec <int64_t> TopolTetra(4);
	
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
		
		for(el=0; el<numelements; el++)
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
	
	ofstream arg("malhaPZ1BC_NW.txt");
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
	int64_t nnodes = nodevec.NElements();
	int64_t in;
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
		ofstream vtkfile("partition.vtk");
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
