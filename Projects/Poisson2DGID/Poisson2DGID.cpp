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

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "pzskylstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpoisson3d.h"
#include <time.h>
#include <stdio.h>
#include "pzpostprocanalysis.h"
#include "TPZSSpStructMatrix.h"

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


void PosProcess(int Dimension, TPZAnalysis &an, std::string plotfile, int div);
void RefinamentoUniforme(TPZGeoMesh  *gmesh, int nh);
void Run(int PolynomialOrder, int  Href, std::string GeoGridFile, int div);



int main()
{
	#ifdef LOG4CXX
		InitializePZLOG();
	#endif
    
	//	Files to read
	std::string GeoGridFile;
	GeoGridFile = "PoissonGIDMeshVar.dump";

	int Href = 2;
	int PolynomialOrder = 1;
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
//    {
//        //    Print Geometrical Base Mesh
//        std::ofstream arg1("BaseGeoMesh.txt");
//        gmesh->Print(arg1);
//        std::ofstream file1("BaseGeoMesh.vtk");
//        //    In this option true -> let you use shrink paraview filter
//        //    PrintGMeshVTK(gmesh,file1);
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);
//    }
//
//    //    Modifying Geometric Mesh
//
    
    RefinamentoUniforme(gmesh, Href);

//    {
//        //    Print Geometrical refined Base Mesh
//        std::ofstream arg1("RefineGeoMesh.txt");
//        gmesh->Print(arg1);
//        std::ofstream file1("RefineGeoMesh.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);
//    }
	
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(myreader.fProblemDimension);
	cmesh->SetAllCreateFunctionsContinuous();
	TPZVec <TPZMatPoisson3d *> materialist(myreader.MatNumber,0);
	TPZVec <TPZMaterial *> MatList(myreader.MatNumber,0);
	
	

//	This part depends on Material parameters	
	for (int imat = 0 ; imat < myreader.MatNumber; imat++)
	{
		
		REAL conv1 = myreader.fMaterialDataVec[imat].fProperties[1];
		TPZVec<REAL> convdir1(3,0);
		REAL flux1 = myreader.fMaterialDataVec[imat].fProperties[5];	
		
		for (int i=0; i < 3; i++) 
		{
			convdir1[i]=myreader.fMaterialDataVec[imat].fProperties[i+2];
		}	
		
		materialist[imat] = new TPZMatPoisson3d(myreader.fMaterialDataVec[imat].fMatID,myreader.fProblemDimension); 
		materialist[imat]->SetParameters(myreader.fMaterialDataVec[imat].fProperties[0], conv1, convdir1);// K
		materialist[imat]->SetInternalFlux(flux1);
		MatList[imat] = materialist[imat];
		cmesh->InsertMaterialObject(MatList[imat]);	
		
		for (int ibc = 0; ibc < myreader.BCNumber ; ibc++) 
		{
			if(imat+1 == int(myreader.fBCMaterialDataVec[ibc].fProperties[0])) 
			{	
				TPZFMatrix<STATE> val11(1,1,myreader.fBCMaterialDataVec[ibc].fProperties[2]), val21(1,1,myreader.fBCMaterialDataVec[ibc].fProperties[3]);	
				TPZMaterial * BC = materialist[imat]->CreateBC(MatList[imat], myreader.fBCMaterialDataVec[ibc].fMatID,int(myreader.fBCMaterialDataVec[ibc].fProperties[1]), val11, val21);
				cmesh->InsertMaterialObject(BC);	
			}
		}
	}
//	End: This part depends on Material parameters	
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	{	
				
        cmesh->SetDimModel(2);
		cmesh->AutoBuild();
	}
	

	TPZAnalysis MyAnalysis (cmesh);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		cmesh->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
	
#endif

	TPZSymetricSpStructMatrix strskyl(cmesh);
	MyAnalysis.SetStructuralMatrix(strskyl);
    strskyl.SetNumThreads(0);
	TPZStepSolver<STATE> direct;
	direct.SetDirect(ELDLt);
	MyAnalysis.SetSolver(direct);
	MyAnalysis.Run();
    std::string plotfile("MyProblemSolution.vtk");
    PosProcess(myreader.fProblemDimension,MyAnalysis, plotfile, div);
    
    TPZPostProcAnalysis post_an;
    post_an.SetCompMesh(MyAnalysis.Mesh());
    
    TPZManVector<int,10> post_pro(3);
    post_pro[0] = 1;
    post_pro[1] = 2;
    post_pro[2] = 3;
    TPZManVector<std::string,10> var_names(1);
    var_names[0] = "Pressure";
    post_an.SetPostProcessVariables(post_pro, var_names);
    
    TPZFStructMatrix structmatrix(post_an.Mesh());
    structmatrix.SetNumThreads(0);
    post_an.SetStructuralMatrix(structmatrix);
    post_an.TransferSolution();
    
    std::string plotfile_post("MyProblemSolution_post.vtk");
    PosProcess(myreader.fProblemDimension,post_an, plotfile_post, div);
    
    return;
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
	TPZManVector<std::string,10> scalnames(1), vecnames(0);
	scalnames[0] = "Pressure";	
//    vecnames[0]= "MinusKGradU";
	
	
	const int dim = Dimension;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("Mesh.txt");
	an.Print("nothing",out);
}

