#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "pzl2projection.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "tpzlinearwave.h"
#include <time.h>
#include <stdio.h>

// Using Log4cXX as logging tool
//
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poisson3d"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poisson3d.data"));
#endif
//
// End Using Log4cXX as logging tool

TPZCompMesh * ComputationalMesh(TPZGeoMesh * gmesh, int p);
void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZCompMesh* CMesh);
void PosProcess(int Dimension, TPZAnalysis &an, std::string plotfile, int div);
void RefinamentoUniforme(TPZGeoMesh  *gmesh, int nh);
void Run(int PolynomialOrder, int  Href, std::string GeoGridFile, int div);
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh);
TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);



int main()
{
	#ifdef LOG4CXX
		InitializePZLOG();
	#endif		
/*
    TPZReadGIDGrid readgid;
    std::string name("oneD.dump");
    TPZGeoMesh *gmesh = readgid.GeometricGIDMesh(name);
    gmesh->Print();
    delete gmesh;
 */
	//	Files to read
	std::string GeoGridFile;	
//	GeoGridFile = "SQDomain.dump";
	GeoGridFile = "2Dwellbore.dump";

	int Href = 0;
	int PolynomialOrder = 2;
	int div = 0;	

	// run the problem
	Run(PolynomialOrder, Href, GeoGridFile, div);

	
	return 0;
}

void Run(int PolynomialOrder, int Href, std::string GeoGridFile, int div)

{
    int pOrder = PolynomialOrder;
	TPZReadGIDGrid myreader;
	TPZGeoMesh * gmesh = myreader.GeometricGIDMesh(GeoGridFile);	
	{
		//	Print Geometrical Base Mesh
		std::ofstream arg1("BaseGeoMesh.txt");
		gmesh->Print(arg1);
		std::ofstream file1("BaseGeoMesh.vtk");	
		//	In this option true -> let you use shrink paraview filter
		//	PrintGMeshVTK(gmesh,file1);	
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
	}	
	
	//	Modifying Geometric Mesh

	RefinamentoUniforme(gmesh, Href);
	{
		//	Print Geometrical refined Base Mesh
		std::ofstream arg1("RefineGeoMesh.txt");
		gmesh->Print(arg1);
		std::ofstream file1("RefineGeoMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
	}	
	
	TPZCompMesh * cmesh = ComputationalMesh(gmesh,PolynomialOrder);

	TPZAnalysis *MyAnalysis = new TPZAnalysis(cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	MyAnalysis->SetStructuralMatrix(strskyl);
	TPZStepSolver<STATE> direct;
	direct.SetDirect(ELDLt);
	MyAnalysis->SetSolver(direct);
    
    REAL deltaT = 0.5;
    REAL maxTime = 10.0;
    SolveSystemTransient(deltaT, maxTime, MyAnalysis,cmesh);
    
// 	MyAnalysis->Run();
// 	std::string plotfile("MyProblemSolution.vtk");		
// 	PosProcess(myreader.fProblemDimension,*MyAnalysis, plotfile, div);	
	
	
}

void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh){
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gMesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			int GelDimension = gel->Dimension();
			if (GelDimension == 2 || GelDimension == 1) 
			{
				gel->Divide (filhos);
			}
		}//for i
	}//ref
}

void PosProcess(int Dimension, TPZAnalysis &an, std::string plotfile, 	int div){
	TPZManVector<std::string,10> scalnames(2), vecnames(1);
	scalnames[0] = "Pressure";	
    scalnames[1] = "WaveSpeed";
	vecnames[0]= "Velocity";
	
	
	const int dim = Dimension;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("Mesh.txt");
	an.Print("nothing",out);
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * gmesh, int p)
{
    int matid = 1;
    int dim = 2;
    REAL wavespeed = 1.0;

    ///Computational Mesh
    TPZCompEl::SetgOrder(p);
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);    
    cmesh->SetAllCreateFunctionsContinuous();
    
    TPZMaterial * Air = new TPZLinearWave(matid,dim);
    cmesh->InsertMaterialObject(Air);
    
    {
        //Boundary Conditions
        TPZFMatrix<STATE> k1(dim,dim,0.), k2(dim,dim,0.);
        TPZMaterial * BCD = Air->CreateBC(Air, 2, 0, k1, k2);
        cmesh->InsertMaterialObject(BCD);
        
        TPZMaterial * BCN = Air->CreateBC(Air, 3, 1, k1, k2);
        cmesh->InsertMaterialObject(BCN);
    }   
        
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();    
    
    return cmesh;
}

void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZCompMesh* CMesh)
{
    
    TPZFMatrix<STATE> Patn;
    TPZFMatrix<STATE> PatnMinusOne;    
//  {
//      TPZBFileStream load;
//      load.OpenRead("MultiphaseSaturationSol.bin");
//      SolutiontoLoad.Read(load,0);
//      meshvec[2]->LoadSolution(SolutiontoLoad);
//      TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);        
//  }

  
    std::string OutPutFile = "WaveSolution";
    TPZMaterial *mat1 = CMesh->FindMaterial(1);    
    
    TPZLinearWave * material1 = dynamic_cast<TPZLinearWave *>(mat1);  
    //    TPZMultiphase * material2 = dynamic_cast<TPZMultiphase *>(mat2);      
    material1->SetTimeStep(deltaT);
    
    //  Starting Newton Iterations
    TPZFMatrix<STATE> DeltaX = CMesh->Solution();
    TPZFMatrix<STATE> Uatn = CMesh->Solution();
    TPZFMatrix<STATE> Uatk = CMesh->Solution();      
    
    
    REAL TimeValue = 0.0;
    REAL Tolerance = 1.0e-7;
    int cent = 0;
    int MaxIterations = 50;
    TimeValue = cent*deltaT;
    REAL NormValue =1.0;
    bool StopCriteria = false;
    TPZFMatrix<STATE> RhsAtnMinusOne, RhsAtn, RhsAtnPlusOne, Residual;

    
    std::string outputfile;
    outputfile = OutPutFile;
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcess(material1->Dimension(),*NonLinearAn,outputfile,2);      

    std::cout << " Starting the time computations. " << std::endl;  
    while (TimeValue < maxTime)
    {
        
        material1->SetMinusOneState();
        CMesh->LoadSolution(PatnMinusOne);
        NonLinearAn->AssembleResidual();
        RhsAtnMinusOne = NonLinearAn->Rhs();

        material1->SetNState();
        CMesh->LoadSolution(Patn);        
        NonLinearAn->AssembleResidual();
        RhsAtn = NonLinearAn->Rhs();        
        
        material1->SetPlusOneState();
        CMesh->LoadSolution(Patn);        
        NonLinearAn->Assemble();
        RhsAtnPlusOne = NonLinearAn->Rhs();
        
        Residual= RhsAtnMinusOne + RhsAtn + RhsAtnPlusOne;       
        NormValue = Norm(Residual);
        


        
        int iterations= 0;      
        while (NormValue > Tolerance)
        {       
            
            Residual*=-1.0;
            NonLinearAn->Rhs()=Residual;
            NonLinearAn->Solve();           
            DeltaX = NonLinearAn->Solution();
            Uatk = (Uatn + DeltaX);
            
            
            CMesh->LoadSolution(Uatn + DeltaX);          
            
#ifdef LOG4CXX
            if(logdata->isDebugEnabled())
            {
                std::stringstream sout;
                sout.precision(20);
                Residual.Print(sout);
                Uatk.Print(sout);       
                LOGPZ_DEBUG(logdata,sout.str());
            }
#endif          


            material1->SetPlusOneState();
            NonLinearAn->Assemble();
            RhsAtnPlusOne = NonLinearAn->Rhs();
            Residual= RhsAtnMinusOne + RhsAtn + RhsAtnPlusOne;
            NormValue = Norm(Residual); 
                
            
#ifdef LOG4CXX
            if(logdata->isDebugEnabled())
            {
                std::stringstream sout;
                sout.precision(15);             
                Uatk.Print(sout);
                Residual.Print("Res = ",sout,EMathematicaInput);
                LOGPZ_DEBUG(logdata,sout.str());
            }
#endif      

            
        
            
            iterations++;
            std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;
            if (iterations == MaxIterations) 
            {
                StopCriteria = true;
                std::cout << " Time Step number = : " << iterations  << "\n Exceed max iterations numbers = : " << MaxIterations <<  std::endl;                 
                break;
            }

                
            Uatn = Uatk;
            
        }   

        outputfile = OutPutFile;
        std::stringstream outputfiletemp;
        outputfiletemp << outputfile << ".vtk";
        std::string plotfile = outputfiletemp.str();
        PosProcess(material1->Dimension(),*NonLinearAn,outputfile,2);      
        
        if (StopCriteria) {
            std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;       
            break;
        }
        
        cent++;
        TimeValue = cent*deltaT;
        
        std::cout << " Time Step :  " << cent  << "  Time :  " << TimeValue <<  std::endl; 
        
        PatnMinusOne = Patn;
        Patn = Uatk;
        
    }
   
}


// Setting up initial conditions

void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
    
    TPZSkylineStructMatrix full(Cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();

}

TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure);
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    
    return cmesh;
    
}

void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    disp[0] = 1.0;
    
}
