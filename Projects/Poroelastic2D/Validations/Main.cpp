	// Finished

	#ifdef HAVE_CONFIG_H
	#include <config.h>
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
	#include "../MeshGeneration.h"
	#include "../AnalyticalFunctions.h"
	#include "../ElasticMatInterface2D.h"
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

	#include "pzpoisson3d.h"
	#include "pzpoisson3dreferred.h"

	#include "pzelasmat.h"
	#include "pzmultiphysicselement.h"
	#include "pzmultiphysicscompel.h"
	#include "pzbuildmultiphysicsmesh.h"
	#include "TPZSpStructMatrix.h"
	#include "../pzporoelastic2d.h"
	#include "pzlog.h"
	#include <iostream>
	#include <string>
	#include "TPZVTKGeoMesh.h"
	#include "pzfunction.h"


	#include <cmath>
	#include <set>


	// Setting 1D  Validation:
	// This Validation is give on: Finite element analysis of poro-elastic consolidation in porous media: Standard and mixed approaches
	// Authors: Johannes Korsawe, Gerhard Starke, Wenqing Wang.

	// Task to do
	// Document all the code
	// Include all BC combinations
	// Include Finite elasticity
	// Conservation mass mini benchmark (tiago's suggestion)
	// How I can divide the code in cases for the realistic analysis?
	// Write Dimensionless Poroelastic Formulation


	// Using Log4cXX as logging tool
	//
	#ifdef LOG4CXX
	static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
	#endif

	#ifdef LOG4CXX
	static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
	#endif
	//
	// End Using Log4cXX as logging tool

	// Defintions of names space
	using namespace std;
	using namespace pzgeom;	

	// Defitions of Material index
	// Internal Materials
	// In this case just one material is used
	const int TotalMats = 1; 
	const int matId = 1;

	// Boundary Conditions Materials
	// 2D Rectangular Domain
	const int bcBottom = -1;
	const int bcRight = -2;
	const int bcTop = -3;
	const int bcLeft = -4;
	const int yfixedPoints = -6;
	const int xfixedPoints = -7;
	const int WellLine = -5;


	// Boundary Conditions Definitions
	// This Definitions are related with the Boudary conditions convention for multphysics simulation (in Poroelastic U/P formulation problems 3 boundary conditions are included) 
	// not defined
	const int dirichlet = 0;
	const int neumann = 1;


	// Defintions of Implemented Methods

	//	This Create Computational Meshes
	TPZCompMesh *MalhaCompPressao(TPZGeoMesh * gmesh,int pOrder);
	TPZCompMesh *MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool Initialstress);
	TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial);

	//	This Solve Different analysis
	void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
	void SolveSistTransient(REAL delt,REAL maxtime,TPZFMatrix<REAL> matK1, TPZAutoPointer < TPZMatrix<REAL> > matK2, TPZFMatrix<REAL> fvec, TPZFMatrix<REAL> &Initialsolution, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZAnalysis &an,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);

	//	This are tools for spatial and polinomyal refinement and Postprocess of solutions 
	void PosProcess(TPZAnalysis &an, std::string plotfile);
	void PosProcess2(TPZAnalysis &an, std::string plotfile);
	void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
	void RefinamentoUniforme(TPZGeoMesh  *gMesh, int nh);
	void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh, int MatId, int indexEl);
	void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
	void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
	void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);
	void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file);

	int main(int argc, char *argv[])
	{
	#ifdef LOG4CXX
		InitializePZLOG("../mylog4cxx.cfg");
	#endif		
		
		
		// GEOMETRICAL MESH CREATION	
		
		//	polynomial base degree approach
		int p=2;
		
		//	Rectangular geometric mesh
		
		MeshGeneration mygenerator;
		TPZVec <int> matIdlist(8,0);
		matIdlist[0]= matId;
		matIdlist[1]= bcBottom;
		matIdlist[2]= bcRight;
		matIdlist[3]= bcTop;
		matIdlist[4]= bcLeft;
		matIdlist[5]= WellLine;
		matIdlist[6]= xfixedPoints;
		matIdlist[7]= yfixedPoints;	
		REAL ModelLenght = 1.0;
		REAL ModelHight = 1.0;
		REAL Modelthickness = 1.0;
		mygenerator.Setdimensions(ModelLenght, ModelHight, Modelthickness);		
		TPZGeoMesh * gmesh = mygenerator.MalhaGeom(matIdlist);
		
		//	End Rectangular geometric mesh
		
		// END GEOMETRICAL MESH CREATION	
		
		
		//	Print Geometrical Base Mesh
		ofstream arg1("BaseGeoMesh.txt");
		gmesh->Print(arg1);
		ofstream file1("BaseGeoMesh.vtk");	
		//	In this option true -> let you use shrink paraview filter
		//	shrink filter in paraview just use for "amazing" visualization
		//	PrintGMeshVTK(gmesh,file1);	
		TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file1, true);	
		
		
		//	First computational mesh
		
		// initialization problem ...
		bool Initialstress = false;
		TPZCompMesh * cmesh1 = MalhaCompElast(gmesh, p+1, Initialstress);
		//	Print First computational mesh
		ofstream arg2("cmeshElasticity.txt");
		cmesh1->Print(arg2);
		
		//	Second computational mesh
		TPZCompMesh * cmesh2= MalhaCompPressao(gmesh, p);
		//	Print Second computational mesh
		ofstream arg3("cmeshDiffusion.txt");
		cmesh2->Print(arg3);	
		
		//	Clear reference of the geometric mesh to cmesh1
		gmesh->ResetReference();
		cmesh1->LoadReferences();
		
		//	Using Uniform Refinement
		RefinUniformElemComp(cmesh1,3);
		cmesh1->AdjustBoundaryElements();
		cmesh1->CleanUpUnconnectedNodes();
		
		ofstream arg4("cmeshElasticityRef.txt");
		cmesh1->Print(arg4);
		ofstream arg5("GeoMeshElasticity.txt");
		gmesh->Print(arg5);
		// Vtk visualization Elasticity Geometric Mesh
		ofstream file3("GeoMeshElasticity.vtk");
		PrintGMeshVTK(gmesh, file3);
		
		// Clear reference of the geometric mesh to cmesh2
		gmesh->ResetReference();
		cmesh2->LoadReferences();
		
		//	Using Uniform Refinement for second computational problem	
		//	Using Uniform Refinement
		RefinUniformElemComp(cmesh2,3);
		cmesh2->AdjustBoundaryElements();
		cmesh2->CleanUpUnconnectedNodes();
		cmesh2->ExpandSolution();
		
		ofstream arg6("cmeshDiffusionRef.txt");
		cmesh2->Print(arg6);
		ofstream arg7("GeoMeshDiffusion.txt");
		gmesh->Print(arg7);
		ofstream file5("GeoMeshDiffusion.vtk");
		PrintGMeshVTK(gmesh, file5);
		
		
		//	Solving First Computational problem
		//	TPZAnalysis an1(cmesh1);
		//	SolveSist(an1, cmesh1);
		//	std::string plotfile("SolutionElasticity.vtk");		
		//	PosProcess2(an1, plotfile);
		
		
		//	Solving Second Computational problem
		TPZAnalysis an2(cmesh2);
		SolveSist(an2, cmesh2);
		//	std::string plotfile2("SolutionDiffusion.vtk");
		//	PosProcess(an2, plotfile2);
		
		
		//	Solving Multiphysic Computational problem
		//	TPZVec<TPZCompMesh *> meshvec(2);
		//	meshvec[0] = cmesh1;
		//	meshvec[1] = cmesh2;
		//	TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,1);
		//	 
		//	ofstream file6("GeoMeshMultphysics.vtk");
		//	PrintGMeshVTK(gmesh, file6);
		//	 
		//	TPZAnalysis an(mphysics);
		//	SolveSist(an, mphysics);
		//	std::string plotfile3("SolutionMultphysics.vtk");
		//	PosProcessMultphysics(meshvec,mphysics, an, plotfile3);
		
		
		// This part is for Poroelastic analysis	
		
		//	Set initial conditions for pressure
		// constan initial pressure
		REAL InitialPressure = 0.0;
		int nrs = an2.Solution().Rows();
		TPZFMatrix<REAL> solucao1(nrs,1,InitialPressure);
		cmesh2->Solution() = solucao1;
		
		TPZVec<TPZCompMesh *> meshvec(2);
		meshvec[0] = cmesh1;
		meshvec[1] = cmesh2;
		TPZVec <TPZPoroElastic2d *>  materialist(TotalMats,0) ;
		TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh,meshvec,materialist);
		
		
		// time control
		REAL hour = 3600;
		REAL day = 86400;
		REAL month = 30*day;
		REAL year = 365*day;
		REAL delta = 0.99999;
		REAL MaxTime = 100.0;
		
		// Improve this
		for(int imat = 0; imat < TotalMats; imat++)
		{	
			materialist[imat]->SetTimeStep(delta);		
			materialist[imat]->SetTimeValue(0.0);
			//Criando matriz K2
			materialist[imat]->SetLastState();
		}
		
		TPZAnalysis an(mphysics);
		
		//Setting initial coditions 
		TPZFMatrix<REAL> Initialsolution = an.Solution();
		std::string output;
		output = "TransientSolution";
		std::stringstream outputfiletemp;
		outputfiletemp << output << ".vtk";
		std::string plotfile = outputfiletemp.str();
		//Print Initial conditions
		PosProcessMultphysics(meshvec,mphysics,an,plotfile);	
		
		TPZSpStructMatrix matsp(mphysics);
		//TPZSkylineStructMatrix matsp(mphysics);	
		std::set< int > materialid;
		materialid.clear();	
		//	int matid = 1;
		//	materialid.insert(imat);
		//	matsp.SetMaterialIds (materialid);
		// Improve this
		for(int imat = 0; imat < TotalMats; imat++)
		{	
			materialid.insert(imat+1);	
		}
		// Inserting Material ids for matrix K2 creation
		matsp.SetMaterialIds(materialid);		
		
		TPZAutoPointer<TPZGuiInterface> guiInterface;
		TPZFMatrix<REAL> Un;
		TPZAutoPointer <TPZMatrix<REAL> > matK2 = matsp.CreateAssemble(Un,guiInterface);
		
	#ifdef LOG4CXX
		if(logdata->isDebugEnabled())
		{
			//		std::stringstream sout;
			//		matK2->Print("K2 = ", sout,EMathematicaInput);
			//		LOGPZ_DEBUG(logdata,sout.str())
		}
	#endif
		
		//Criando matriz K1
		// Improve this
		for(int imat = 0; imat < TotalMats; imat++)
		{	
			materialist[imat]->SetCurrentState();
		}	
		TPZSkylineStructMatrix matsk(mphysics);
		TPZFMatrix<REAL> matK1;	
		TPZFMatrix<REAL> fvec; //vetor de carga
		an.SetStructuralMatrix(matsk);//	Set the structtural sky line matrix to our analysis
		
		
		TPZStepSolver<REAL> step; //Create Solver object
		step.SetDirect(ELDLt); //	Symmetric case
		//step.SetDirect(ELU);
		an.SetSolver(step); //	Set solver
		an.Run(); //	Excecute an analysis to obtain the Rhs vector (In this case we start with zero initial values)
		
		matK1 = an.StructMatrix(); //Storage the global matrix and load vector
		fvec = an.Rhs();
		
		// This code identify singular blocks
		TPZStepSolver<REAL> & temp = dynamic_cast<TPZStepSolver<REAL> &> (an.Solver());
		std::list <int> & zeropivot = temp.Singular(); 
		if (zeropivot.size()) 
		{
			int eq = * zeropivot.begin();
			an.Rhs().Zero();
			an.Rhs()(eq,0) = -10000.0;
			an.Solve();
			TPZFMatrix<REAL> TempSolution = an.Solution();
			
	#ifdef LOG4CXX
			// Print the temporal solution
			if(logdata->isDebugEnabled())
			{
				std::stringstream sout;
				TempSolution.Print("SingularNodes = ", sout,EMathematicaInput);
				LOGPZ_DEBUG(logdata,sout.str())
			}
	#endif	
			std::string output;
			output = "SingularNodes";
			std::stringstream outputfiletemp;
			outputfiletemp << output << ".vtk";
			std::string plotfile = outputfiletemp.str();
			PosProcessMultphysics(meshvec,mphysics,an,plotfile);
			
		}
		
	#ifdef LOG4CXX
		if(logdata->isDebugEnabled())
		{
			// base system to invert
			// just one for checking purpose
			//		std::stringstream sout;
			//		an.Solver().Matrix()->Print("K1 = ", sout,EMathematicaInput);
			//		fvec.Print("fvec = ", sout,EMathematicaInput);		
			//		//Print the temporal solution
			//		Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
			//		TPZFMatrix<REAL> Temp;
			//		TPZFMatrix<REAL> Temp2;
			//		matK2->Multiply(Initialsolution,Temp);
			//		Temp.Print("Temp K2 = ", sout,EMathematicaInput);	
			//		LOGPZ_DEBUG(logdata,sout.str())
		}
	#endif
		
		///start transient problem
		SolveSistTransient(delta,MaxTime,matK1, matK2, fvec, Initialsolution, materialist, an, meshvec,  mphysics);
		
		
		return EXIT_SUCCESS;
	}


	TPZCompMesh*MalhaCompElast(TPZGeoMesh * gmesh,int pOrder, bool Initialstress)
	{
		int dim = 2;
		int planestress = 0;	
		// All material with equal properties
		// Needed include diferents materials and properties via keyword file	
		REAL gravity = 0.0; // SI system
		REAL overburdendepth = 2000.0; // SI system
		REAL layerthickness = 10.0;  // SI system	
		TPZVec <REAL> E(TotalMats,100.0);
		TPZVec <REAL> poisson(TotalMats,0.35);
		TPZVec <REAL> rockrho(TotalMats,2330.0);	
		TPZVec <REAL> force(2,0.);
		TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);	
		TPZVec <TPZElasticityMaterial *> materialist(TotalMats,0);
		
		TPZCompEl::SetgOrder(pOrder);
		TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDimModel(dim);
		cmesh->SetAllCreateFunctionsContinuous();	
		
		
		
		for(int imat = 0; imat < TotalMats; imat++)
		{	
			
			// No Gravity field
			//		forceVec[imat][1]=-gravity*rockrho[imat];
			materialist[imat] = new TPZElasticityMaterial(imat+1, E[imat], poisson[imat], forceVec[imat][0], forceVec[imat][1], planestress); 
			
			TPZMaterial * mat(materialist[imat]);
			cmesh->InsertMaterialObject(mat);
			
			// This is a dummy boundary conditions
			// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
			TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
			
			// Top Boundary			
			TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTop,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCTOP);
			
			// Bottom Boundary			
			TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottom,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCBOTT);
			
			// Left Boundary			
			TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeft,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCLEFT);
			
			// Right Boundary			
			TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRight,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCRIGHT);	
			
		}
		cmesh->AutoBuild();
		
		
		return cmesh;
		
	}

	TPZCompMesh*MalhaCompPressao(TPZGeoMesh * gmesh, int pOrder)
	{
		/// criar materiais
		int dim = 2;
		// All material with equal properties
		// Needed include diferents materials and properties via keyword file
		TPZVec <REAL> convdir(3,0.);
		TPZVec <REAL> diff(TotalMats,0.1);
		TPZVec <REAL> conv(TotalMats,0.0);
		TPZVec <TPZVec< REAL > > convdirVec(TotalMats,convdir);
		TPZVec <REAL> flux(TotalMats,0.0);
		
		TPZVec <TPZMatPoisson3d *> materialist(TotalMats,0);	
		TPZCompEl::SetgOrder(pOrder);
		TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
		cmesh->SetDimModel(dim);
		cmesh->SetAllCreateFunctionsContinuous();	
		
		for(int imat = 0; imat < TotalMats; imat++)
		{
			materialist[imat] = new TPZMatPoisson3d(imat+1,dim); 
			TPZMaterial * mat(materialist[imat]);
			
			materialist[imat]->SetParameters(diff[imat], conv[imat], convdirVec[imat]);
			materialist[imat]->SetInternalFlux(flux[imat]);
			materialist[imat]->NStateVariables();
			
			cmesh->InsertMaterialObject(mat);
			
			// This is a dummy boundary conditions
			// Here is necessary to include the same set of computational information for all computational meshes -> Elastic, Pressure, Biot's Poroelasticity
			TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
			
			// Top Boundary			
			TPZMaterial * BCTOP = materialist[imat]->CreateBC(mat, bcTop,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCTOP);
			
			// Bottom Boundary			
			TPZMaterial * BCBOTT = materialist[imat]->CreateBC(mat, bcBottom,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCBOTT);
			
			// Left Boundary			
			TPZMaterial * BCLEFT = materialist[imat]->CreateBC(mat, bcLeft,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCLEFT);
			
			// Right Boundary			
			TPZMaterial * BCRIGHT = materialist[imat]->CreateBC(mat, bcRight,dirichlet, val1, val2);
			cmesh->InsertMaterialObject(BCRIGHT);	
		}
		
		cmesh->AutoBuild();
		
		return cmesh;
	}

	TPZCompMesh *MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZVec <TPZPoroElastic2d  * > &mymaterial)
	{
		//Creating computational mesh for multiphysic elements
		int dim = 2;
		int planestress = 0; // This is a Plain strain problem
		
		gmesh->ResetReference();
		TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
		mphysics->SetAllCreateFunctionsMultiphysicElem();
		TPZAutoPointer<TPZFunction<STATE> > TimeDepForcingF;
		TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact;
		TimeDepFExact = new TPZDummyFunction<STATE>(ExactSolution1D);
		
		// Definitions
		REAL lamb = 8.3e9;						//	[Pa]
		REAL lambu = 13.4e9;					//	[Pa]
		REAL alpha = 27;						//	[-]
		REAL G= 5.5e9;							//	[Pa]
		REAL rhof = 2000.0;						//	[kg/m3]
		REAL poisson = 0.4;						//	[-]
		REAL Eyoung = 6.0e3;					//	[Pa]			
		REAL Ssig = ((pow(alpha,2))/((lambu-lamb)))*((lambu+2.0*G)/(lamb+2.0*G));				//	[-]
		REAL Se = 40.0;
		REAL rockrho = 2500.0;					//	[kg/m3]
		REAL gravity = 0.0;						//	[m/s2]
		REAL c = 0.083;							//	[m2/s]	
		REAL visc = 1.0;						//	[Pa.s]
		REAL perm = 0.2;						//	[m2]
		REAL qo = -0.001;
		REAL Bo = 1.0;
		REAL PI = atan(1.)*4.;
		REAL Yforce = -1000.0;
		TPZVec <REAL> force(2,0.);	
		REAL fx=0.0;
		REAL fy=0.0;
		
		
		// Definitions for paramter vectors Al materials have the same porperties
		TPZVec <REAL> lambVec(TotalMats,lamb);
		TPZVec <REAL> lambuVec(TotalMats,lambu);
		TPZVec <REAL> alphaVec(TotalMats,alpha);
		TPZVec <REAL> GVec(TotalMats,G);
		TPZVec <REAL> rhofVec(TotalMats,rhof);
		TPZVec <REAL> poissonVec(TotalMats,poisson);
		TPZVec <REAL> EyoungVec(TotalMats,Eyoung);
		TPZVec <REAL> SsigVec(TotalMats,Ssig);
		TPZVec <REAL> SeVec(TotalMats,Se);
		TPZVec <REAL> rockrhoVec(TotalMats,rockrho);
		TPZVec <REAL> cVec(TotalMats,c);
		TPZVec <REAL> viscVec(TotalMats,visc);
		TPZVec <REAL> permVec(TotalMats,perm);
		TPZVec <REAL> qoVec(TotalMats,qo);
		TPZVec <REAL> BoVec(TotalMats,Bo);
		TPZVec <REAL> YforceVec(TotalMats,Yforce);
		TPZVec <TPZVec< REAL > > forceVec(TotalMats,force);
		for(int imat = 0; imat < TotalMats; imat++)
		{	
			
			// No Gravity field
			//		forceVec[imat][1]=gravity*rockrhoVec[imat];		
			
			mymaterial[imat] = new TPZPoroElastic2d (imat+1, dim);
			mymaterial[imat]->SetParameters(EyoungVec[imat], poissonVec[imat],forceVec[imat][0], forceVec[imat][1]);
			mymaterial[imat]->SetParameters(permVec[imat],viscVec[imat]);
			mymaterial[imat]->SetfPlaneProblem(planestress);
			mymaterial[imat]->SetBiotParameters(alpha,Se);
			//	mymaterial[imat]->SetTimeDependentForcingFunction(TimeDepForcingF);
			mymaterial[imat]->SetTimeDependentFunctionExact(TimeDepFExact);		
			
			ofstream argm("mymaterial.txt");
			mymaterial[imat]->Print(argm);
			TPZMaterial * mat(mymaterial[imat]);
			mphysics->InsertMaterialObject(mat);
			
			///--- --- Boundary Conditions
			// Squared Model 1D compactation
			
			// Elastic problem -> 1 ,  Diffsion problem -> 2
			
			// Top Boundary	-3
			// Setting free Boundary condition N1N2 means: Elastic -> Neumman, Diffusion -> Dirichlet
			int N1D2 = 10;
			TPZFMatrix<REAL> val13(3,2,0.), val23(3,1,0.);
			REAL SigTopx	= 0.0;
			REAL SigTopy	= 0.0;			
			REAL pTop		= 1000.0;
			val23(0,0)=SigTopx;
			val23(1,0)=SigTopy;			
			val23(2,0)=pTop;
			TPZMaterial * BCTOP = mymaterial[imat]->CreateBC(mat,bcTop, N1D2, val13, val23);
			mphysics->InsertMaterialObject(BCTOP);			
			
			// Bottom Boundary -1					
			/// Setting free displacement non permeable surface Boundary condition DYFX1N2 means: Elastic -> Dirichlet Y component, Diffusion -> Neumman
			TPZFMatrix<REAL> val11(3,2,0.), val21(3,1,0.);
			int DYFX1N2 = 200; 			
			REAL uDBottx	=0.0;
			REAL uDBotty	=0.0;
			REAL fluxBott	=0.0;
			val21(0,0)=uDBottx;
			val21(1,0)=uDBotty;
			val21(2,0)=fluxBott;
			TPZMaterial * BCBOTT = mymaterial[imat]->CreateBC(mat,bcBottom, DYFX1N2, val11, val21);
			mphysics->InsertMaterialObject(BCBOTT);
			
			// Letf Boundary -4		
			/// Setting free displacement non permeable surface Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman
			int DXFY1N2 = 300;
			TPZFMatrix<REAL> val14(3,2,0.), val24(3,1,0.);
			REAL uDLeftx	= 0.0;
			REAL uDLefty	= 0.0;			
			REAL fluxLeft	= 0.0;
			val24(0,0)=uDLeftx;
			val24(1,0)=uDLefty;			
			val24(2,0)=fluxLeft;
			TPZMaterial * BCLEFT = mymaterial[imat]->CreateBC(mat, bcLeft, DXFY1N2, val14, val24);
			mphysics->InsertMaterialObject(BCLEFT);
			
			// Right Boundary -2			
			/// Setting free displacement non permeable surface  Boundary condition DXFY1N2 means: Elastic -> Dirichlet X component, Pressure -> Neumman			
			TPZFMatrix<REAL> val12(3,2,0.), val22(3,1,0.);
			
			REAL uDRightx	= 0.0;
			REAL uDRighty	= 0.0;			
			REAL fluxRight	= 0.0;
			val22(0,0)=uDRightx;
			val22(1,0)=uDRighty;			
			val22(2,0)=fluxRight;
			TPZMaterial * BCRIGHT = mymaterial[imat]->CreateBC(mat, bcRight, DXFY1N2, val12, val22);
			mphysics->InsertMaterialObject(BCRIGHT);
			
		}
		
		mphysics->AutoBuild();
		mphysics->AdjustBoundaryElements();
		mphysics->CleanUpUnconnectedNodes();
		
		// Creating multiphysic elements into mphysics computational mesh
		TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
		TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
		TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
		
	#ifdef LOG4CXX
		{
			//        std::stringstream sout;
			//        mphysics->Print(sout);
			//        LOGPZ_DEBUG(logger, sout.str())
		}
	#endif
		ofstream arg("Mphysic.txt");
		mphysics->Print(arg);
		
		return mphysics;
	}

	#define VTK
	void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
	{			
		//TPZBandStructMatrix full(fCmesh);
		TPZSkylineStructMatrix full(fCmesh); //caso simetrico
		an.SetStructuralMatrix(full);
		TPZStepSolver<REAL> step;
		step.SetDirect(ELDLt); //caso simetrico
		//step.SetDirect(ELU);
		an.SetSolver(step);
		an.Run();
		
		//	ofstream file("Solution.out");
		//	an.Solution().Print("solution", file); 
	}

	void SolveSistTransient(REAL deltime,REAL maxtime,TPZFMatrix<REAL> matK1, TPZAutoPointer < TPZMatrix<REAL> > matK2, TPZFMatrix<REAL> fvec, TPZFMatrix<REAL> &Initialsolution, TPZVec <TPZPoroElastic2d  * > &mymaterial, TPZAnalysis &an,TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics){
		
		int nrows;
		nrows = matK2->Rows();
		TPZFMatrix<REAL> TotalRhs(nrows,1,0.0);
		TPZFMatrix<REAL> TotalRhstemp(nrows,1,0.0);
		TPZFMatrix<REAL> Lastsolution = Initialsolution;
		
		std::string outputfile;
		outputfile = "TransientSolution";
		
		REAL delt = deltime;
		REAL Maxtime = maxtime;
		REAL TimeValue = 0.0;
		int cent = 1;
		TimeValue = cent*delt; 
		while (TimeValue < Maxtime)
		{	
			// This time solution i for Transient Analytic Solution
			for(int imat = 0; imat < TotalMats; imat++)
			{
				mymaterial[imat]->SetTimeValue(TimeValue);
			}
			matK2->Multiply(Lastsolution,TotalRhstemp);
			TotalRhs = fvec + TotalRhstemp;
			an.Rhs() = TotalRhs;
			
			
	#ifdef LOG4CXX
			if(logdata->isDebugEnabled())
			{
				//	Print the temporal solution
				//			std::stringstream sout;
				//			matK1.Print("Temporal Solution = ", sout,EMathematicaInput);
				//			TotalRhstemp.Print("Temporal Solution = ", sout,EMathematicaInput);
			}
	#endif	
			
			an.Solve(); //	Solve the current Linear system
			Lastsolution = an.Solution(); //	Save the current solution
			//Lastsolution.Print("Sol = ");
			
	#ifdef LOG4CXX
			//Print the temporal solution
			//		if(logdata->isDebugEnabled()){
			//			std::stringstream sout;
			//			Lastsolution.Print("Temporal Solution = ", sout,EMathematicaInput);
			//			LOGPZ_DEBUG(logdata,sout.str())
			//		}
	#endif
			
			//General post-processing
			//TPZBuildMultiphysicsMesh * Objectdumy;
			//Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
			std::stringstream outputfiletemp;
			outputfiletemp << outputfile << ".vtk";
			std::string plotfile = outputfiletemp.str();
			PosProcessMultphysics(meshvec,mphysics,an,plotfile);		
			
			// Total mass calculation
			//
			
			//		int totalel=mphysics->NElements();
			//		int vartointegrate = 13;
			//		REAL TotalMass=0;
			//		
			//		
			//		
			//		for ( int el = 0; el < totalel; el++ )
			//		{
			//		
			//			TPZVec <REAL> result(1,0);
			//			TPZCompEl * celvar = mphysics->ElementVec()[el];
			//			if (celvar->Reference()->MaterialId() > 0 )// avoid bc conditions
			//			{
			//				celvar->Integrate(vartointegrate, result);
			//				TotalMass += result[0];
			//				
			//			}
			//			
			//		}
			//		
			//		cout << "Numerical total mass  " << TotalMass << "  at time " << TimeValue << endl;
			//		cout << "Theorical total mass  " << 0.000001*TimeValue << "  at time " << TimeValue << endl;		
			
			//
			// End TotalMass Calculation
			
			// Next Calculation
			cent++;
			TimeValue = cent*delt;
		}
		
		
	}

	void PosProcess(TPZAnalysis &an, std::string plotfile){
		TPZManVector<std::string,10> scalnames(1), vecnames(1);
		scalnames[0] = "Solution";
		vecnames[0]= "MinusKGradU";
		
		const int dim = 2;
		int div = 2;
		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);
		std::ofstream out("malha.txt");
		an.Print("nothing",out);
	}

	void PosProcess2(TPZAnalysis &an, std::string plotfile){
		TPZManVector<std::string,10> scalnames(3), vecnames(1);
		scalnames[0] = "SigmaX";
		scalnames[1] = "SigmaY";
		scalnames[2] = "Pressure";	
		vecnames[0]= "displacement";
		
		
		const int dim = 2;
		int div = 2;
		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);
		std::ofstream out("malha.txt");
		an.Print("nothing",out);
	}

	void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile)
	{
		TPZBuildMultiphysicsMesh * Objectdumy;
		Objectdumy->TransferFromMultiPhysics(meshvec, mphysics);
		TPZManVector<std::string,10> scalnames(12), vecnames(2);
		scalnames[0] = "SigmaX";
		scalnames[1] = "SigmaY";
		scalnames[2] = "TauXY";	
		scalnames[3] = "DisplacementX";
		scalnames[4] = "DisplacementY";
		scalnames[5] = "SolidPressure";
		scalnames[6] = "FluidPressure";
		scalnames[7] = "EDisplacementX";
		scalnames[8] = "EDisplacementY";
		scalnames[9] = "EPressure";
		scalnames[10] = "ESIGX";
		scalnames[11] = "ESIGY";
		vecnames[0]  = "Displacement";
		vecnames[1]  = "FluxVector";	
		
		
		
		const int dim = 2;
		int div =2;
		an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
		an.PostProcess(div,dim);
		std::ofstream out("malha.txt");
		an.Print("nothing",out);
	}

	void RefinamentoUniforme(TPZGeoMesh *gMesh, int nh){
		for ( int ref = 0; ref < nh; ref++ ){
			TPZVec<TPZGeoEl *> filhos;
			int n = gMesh->NElements();
			for ( int i = 0; i < n; i++ ){
				TPZGeoEl * gel = gMesh->ElementVec() [i];
				if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
			}//for i
		}//ref
	}

	void RefinamentoUniforme(TPZGeoMesh * gMesh, int nh, int MatId, int indexEl){
		for ( int ref = 0; ref < nh; ref++ ){
			TPZVec<TPZGeoEl *> filhos;
			int n = gMesh->NElements();
			for ( int i = 0; i < n; i++ ){
				TPZGeoEl * gel = gMesh->ElementVec() [i];
				if (gel->Dimension() == 2 || gel->Dimension() == 1){
					if (gel->MaterialId()== MatId && gel-> Index()==indexEl){
						gel->Divide (filhos);
					}
				}
			}//for i
		}//ref
	}

	void RefinElemComp(TPZCompMesh  *cMesh, int indexEl){
		
		TPZVec<int > subindex; 
		int nel = cMesh->ElementVec().NElements(); 
		for(int el=0; el < nel; el++){
			TPZCompEl * compEl = cMesh->ElementVec()[el];
			if(!compEl) continue;
			int ind = compEl->Index();
			if(ind==indexEl){
				compEl->Divide(indexEl, subindex, 1);
			}
		}	
	}

	void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv){
		
		TPZVec<int > subindex;
		for (int iref = 0; iref < ndiv; iref++) {
			TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
			int nel = elvec.NElements(); 
			for(int el=0; el < nel; el++){
				TPZCompEl * compEl = elvec[el];
				if(!compEl) continue;
				int ind = compEl->Index();
				compEl->Divide(ind, subindex, 0);
			}
		}
	}

	void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
	{
		file.clear();
		int nelements = gmesh->NElements();
		
		std::stringstream node, connectivity, type;
		
		//Header
		file << "# vtk DataFile Version 3.0" << std::endl;
		file << "TPZGeoMesh VTK Visualization" << std::endl;
		file << "ASCII" << std::endl << std::endl;
		
		file << "DATASET UNSTRUCTURED_GRID" << std::endl;
		file << "POINTS ";
		
		int actualNode = -1, size = 0, nVALIDelements = 0;
		
		for(int el = 0; el < nelements; el++)
		{        
			if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
			{
				continue;
			}
			if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
			{
				continue;
			}
			if(gmesh->ElementVec()[el]->HasSubElement())
			{
				continue;
			}
			
			int elNnodes = gmesh->ElementVec()[el]->NNodes();
			size += (1+elNnodes);
			connectivity << elNnodes;
			
			for(int t = 0; t < elNnodes; t++)
			{
				for(int c = 0; c < 3; c++)
				{
					double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
					node << coord << " ";
				}            
				node << std::endl;
				
				actualNode++;
				connectivity << " " << actualNode;
			}
			connectivity << std::endl;
			
			int elType = -1;
			switch (gmesh->ElementVec()[el]->Type())
			{
				case (ETriangle):
				{
					elType = 5;
					break;                
				}
				case (EQuadrilateral ):
				{
					elType = 9;
					break;                
				}
				case (ETetraedro):
				{
					elType = 10;
					break;                
				}
				case (EPiramide):
				{
					elType = 14;
					break;                
				}
				case (EPrisma):
				{
					elType = 13;
					break;                
				}
				case (ECube):
				{
					elType = 12;
					break;                
				}
				default:
				{
					//ElementType NOT Found!!!
					DebugStop();
					break;    
				}
			}
			
			type << elType << std::endl;
			nVALIDelements++;
		}
		node << std::endl;
		actualNode++;
		file << actualNode << " float" << std::endl << node.str();
		
		file << "CELLS " << nVALIDelements << " ";
		
		file << size << std::endl;
		file << connectivity.str() << std::endl;
		
		file << "CELL_TYPES " << nVALIDelements << std::endl;
		file << type.str();
		
		file.close();
	}


	void PrintRefPatternVTK(TPZAutoPointer<TPZRefPattern> refp, std::ofstream &file)
	{
		TPZGeoMesh *gmesh;
		
		//    RefPatternMesh();
		//TPZGeoMesh * gmesh = refp->Mesh();
		PrintGMeshVTK(gmesh, file);
	}
