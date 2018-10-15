#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include <pzgengrid.h>
#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzelasmat.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"


#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

static bool oldmat = true;

// Defintions of Implemented Methods
TPZCompMesh *ComputationalElasticityMesh(TPZGeoMesh * gmesh,int pOrder);
void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter = 20);

//	This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

//	These are tools for spatial and polynomial refinement and Postprocess of solutions
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
void UniformRefinement(TPZGeoMesh  *gMesh, int nh);
void UniformRefinement(TPZGeoMesh *gMesh, int nh, int MatId);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void ReservoirPressure(const TPZVec<REAL> &x, TPZVec<STATE> &p,  TPZFMatrix<STATE> &gradp);

int main(int argc, char *argv[])
{

    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/Elasticity2D/";
    FileName += "Elasticity2DLog.cfg";
    InitializePZLOG(FileName);
#endif

    std::string GridFileName;
    GridFileName = dirname + "/Projects/Elasticity2D/";
    GridFileName += "Unit.dump";
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    
    TPZVec<int64_t> PointTopology(1);
    PointTopology[0]=0;
    int matPoint = 6;
    int64_t index;
    gmesh->CreateGeoElement(EPoint, PointTopology, matPoint, index);
    
    PointTopology[0]=1;
    int matPoint2 = 7;
    gmesh->CreateGeoElement(EPoint, PointTopology, matPoint2, index);
    
    gmesh->BuildConnectivity();
    

    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMesh.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("GeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }

    int Href = 1;
    int PElasticity = 2;
	UniformRefinement(gmesh, Href);
	
#ifdef LOG4CXX
	{
		//	Print Geometrical refined Base Mesh
        std::ofstream argument("RefinedGeometricMesh.txt");
		gmesh->Print(argument);
        std::ofstream Dummyfile("RefinedGeometricMesh.vtk");
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);	
	}	
#endif
	
	TPZCompMesh * ComputationalMeshElasticity = ComputationalElasticityMesh(gmesh, PElasticity);
	//	Print First computational mesh
    std::ofstream ArgumentElasticity("DumpFolder/ComputationalMeshForElasticity.txt");
	ComputationalMeshElasticity->Print(ArgumentElasticity);
	
	
	// Visualization of computational meshes
      bool mustOptimizeBandwidth = false;
      TPZAnalysis * ElasticAnalysis = new TPZAnalysis(ComputationalMeshElasticity,mustOptimizeBandwidth);
    
      SolveSist(ElasticAnalysis, ComputationalMeshElasticity);
	
	std::string ElasticityOutput;
	ElasticityOutput = "ComputationalMeshElasticity";
	std::stringstream ElasticityOutputfiletemp;
	ElasticityOutputfiletemp << ElasticityOutput << ".vtk";
	std::string ElasticityPlotfile = ElasticityOutputfiletemp.str();
	PostProcessElasticity(*ElasticAnalysis, ElasticityPlotfile);
	
	
    std::cout << "Check:: Calculation finished successfully" << std::endl;
	return EXIT_SUCCESS;
}


TPZCompMesh * ComputationalElasticityMesh(TPZGeoMesh * gmesh,int pOrder)
{
    
    // Getting mesh dimension
    int dim = 2;
    int matId1 = 1;
    
   TPZMatElasticity2D *material;
   material = new TPZMatElasticity2D(matId1);

    
    // Setting up paremeters
    // Plane strain assumption
    material->SetPlaneStrain();
    REAL lamelambda = 0.0e9,lamemu = 0.5e9, fx= 0, fy = 0;
    material->SetParameters(lamelambda,lamemu, fx, fy);
    //material->SetElasticParameters(40.0,0.0);
    REAL Sigmaxx = 0.0, Sigmayx = 0.0, Sigmayy = 0.0, Sigmazz = 0.0;
    material->SetPreStress(Sigmaxx,Sigmayx,Sigmayy,Sigmazz);
//    REAL Alpha = 1.0;
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    TPZAutoPointer<TPZFunction<STATE> > Pressure;
    Pressure = new TPZDummyFunction<STATE>(ReservoirPressure, 5);
    material->SetForcingFunction(Pressure);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond2 = material->CreateBC(material,2,1, val1, val2);

    val2(0,0) = 1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond3 = material->CreateBC(material,3,1, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material,4,1, val1, val2);
    
    val2(0,0) = -1.0*1000.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material,5,1, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond6 = material->CreateBC(material,6,0, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond7 = material->CreateBC(material,7,8, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BCond6);
    cmesh->InsertMaterialObject(BCond7);
    cmesh->AutoBuild();
    return cmesh;
    
}


#define VTK
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh)
{			
	TPZSkylineStructMatrix skymat(fCmesh);
	an->SetStructuralMatrix(skymat);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	an->SetSolver(step);
    int Iterations = 10;
    IterativeProcess(an, std::cout, Iterations);
}

void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(0), vecnames(0);
    
        
    if (oldmat) {
        scalnames.Resize(2);
        vecnames.Resize(1);
        scalnames[0] = "SigmaX";
        scalnames[1] = "SigmaY";
        vecnames[0]= "Displacement";
    }else{
        scalnames.Resize(4);
        vecnames.Resize(1);
        scalnames[0] = "TotStressXX";
        scalnames[1] = "TotStressYY";
        scalnames[2] = "EffStressXX";
        scalnames[3] = "EffStressYY";
        vecnames[0]= "DisplacementTotal";
    }
	
	const int dim = 2;
	int div = 2;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gMesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
}

void UniformRefinement(TPZGeoMesh * gMesh, int nh, int MatId)
{
	for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int64_t n = gMesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
			TPZGeoEl * gel = gMesh->ElementVec() [i];
			if (gel->Dimension() == 2 || gel->Dimension() == 1){
				if (gel->MaterialId()== MatId){
					gel->Divide (filhos);
				}
			}
		}//for i
	}//ref
}

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl)
{
	
	TPZVec<int64_t > subindex; 
	int64_t nel = cMesh->ElementVec().NElements(); 
	for(int64_t el=0; el < nel; el++){
		TPZCompEl * compEl = cMesh->ElementVec()[el];
		if(!compEl) continue;
		int64_t ind = compEl->Index();
		if(ind==indexEl){
			compEl->Divide(indexEl, subindex, 1);
		}
	}	
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
	
	TPZVec<int64_t > subindex;
	for (int64_t iref = 0; iref < ndiv; iref++) {
		TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
		int64_t nel = elvec.NElements(); 
		for(int64_t el=0; el < nel; el++){
			TPZCompEl * compEl = elvec[el];
			if(!compEl) continue;
			int64_t ind = compEl->Index();
			compEl->Divide(ind, subindex, 0);
		}
	}
}

void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter)
{
    int iter = 0;
    REAL error = 1.e10, NormResLambdaLast = 1.e10;;
    const REAL tol = 1.e-5;
    
    int numeq = an->Mesh()->NEquations();
    
    TPZFMatrix<STATE> Uatk0(an->Solution());
    TPZFMatrix<STATE> Uatk(Uatk0),DeltaU(Uatk0);
    if(Uatk0.Rows() != numeq) Uatk0.Redim(numeq,1);
    
    an->Assemble();
    an->Rhs() *= -1.0; //- [R(U0)];
    
    TPZAutoPointer< TPZMatrix<STATE> > matK; // getting X(Uatn)
    
    bool converged = false;
    while(!converged && iter < numiter) {
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            matK=an->Solver().Matrix();
            matK->Print("matK = ", sout,EMathematicaInput);
            an->Rhs().Print("Rhs = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        // Computing Uatk = Uatn + DeltaU;
        an->Solve();
        DeltaU= an->Solution();
        Uatk = Uatk0 + DeltaU;
        
        //Computing ||DeltaU||
        REAL NormOfDeltaU = Norm(DeltaU);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaU.Print("DeltaU = ", sout,EMathematicaInput);
            Uatk.Print("Uatk = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        an->LoadSolution(Uatk); // Loading Uatk
        an->Assemble();
        an->Rhs() *= -1.0; //- [R(U0)];
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            an->Rhs().Print("Res = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        // Computing ||[R(Uatk)]||
        double ResidualNorm = Norm(an->Rhs());
        double norm = NormOfDeltaU; //ResidualNorm;
        out << "Iteration n : " << (iter+1) << " : norms ||DeltaU|| e ||[R(Uatk)]|| : " << NormOfDeltaU << " / " << ResidualNorm << std::endl;
        
        if(norm < tol /*|| NormResLambda < tol*/) {
            out << "\nNewton Converged! Tolerance Of Norm(DeltaU) at n : " << (iter+1) << std::endl;
            out << "Norm ||DeltaU|| - USED : " << NormOfDeltaU << std::endl;
            out << "Norm ||[R(Uatk)]||  : " << ResidualNorm << std::endl;
            converged = true;
        }
        else if( (ResidualNorm - NormResLambdaLast) > 1.e-4 ) {
            out << "\nDivergent Method\n" << "Implement Line Search Please!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
        
        NormResLambdaLast = ResidualNorm;
        error = norm;
        iter++;
        Uatk0 = Uatk;
        out.flush();
    }
    
    if (error > tol) {
        DebugStop(); // Something is very wrong
    }
    
}

void ReservoirPressure(const TPZVec<REAL> &x, TPZVec<STATE> &p,  TPZFMatrix<STATE> &gradp)
{
  p[0] = 1.0e7;
}
