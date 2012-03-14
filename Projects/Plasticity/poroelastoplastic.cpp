//$Id: poroelastoplastic.cpp,v 1.53 2010-06-11 22:13:02 diogo Exp $

#include <iostream>
#include <string>
#include "pzmatrix.h"
#include "TPZTensor.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include <math.h>
#include "pzelasmat.h" 
#include "pzelast3d.h"
#include "TPZMohrCoulomb.h"

#include "TPZTensor.h"
#include "TPZYCDruckerPrager.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic.h"
#include "TPZSandlerDimaggio.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"


#include "TPZDruckerPrager.h"

//#include "TPZPlasticityTest.h"
#include "pzlog.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZYCWillamWarnke.h"

#include "TPZSpStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZTimer.h"


#include "pzelastoplastic2D.h"

#include "pzlog.h"

#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"

#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"

#include "pzsbstrmatrix.h"
#include "pzstepsolver.h"

#include "pzadmchunk.h"
#include "gmres.h"
#include "pzbndcond.h"
#include "pzelast3d.h"
#include "pzblockdiag.h"
#include "pzvisualmatrix.h"
#include "cg.h"
#include "pzseqsolver.h"

#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "TPZParSkylineStructMatrix.h"


#include "TPZParFrontMatrix.h"
#include "TPZFrontMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "BrazilianTestGeoMesh.h"
#include "GeoMeshClass.h"

#ifdef LOG4CXX // LOG4CXX may be defined alone or with LOG4CXX_PLASTICITY. The latter shall not be used alone.
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#endif

//#define BOOST_TEST_DYN_LINK
//
#ifdef LOG4CXX_PLASTICITY
static LoggerPtr logger(Logger::getLogger("plasticity.main"));
static LoggerPtr logger2(Logger::getLogger("plasticity2.main"));
#endif



void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SolutionGrafic(TPZAnalysis& an_p1);
void SetUPPostProcessVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void SetUPPostProcessElasticVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
							int BCId,int BCId2, int nsteps, REAL PGRatio,
							TPZFMatrix & val1Begin, TPZFMatrix & val1End,
							TPZFMatrix & val2Begin, TPZFMatrix & val2End,
							TPZPostProcAnalysis * ppAnalysis, int res);

//thetaRad > 0 -> giro no sentido anti-horário
void RotationMatrix(TPZFMatrix &R, double thetaRad, int axis);
void RotateMatrix(TPZFMatrix &Mat, double thetaRad,int rotateaboutaxes);
void RotateMesh(TPZGeoMesh &geomesh, REAL angleDegree,int rotateabout);
void LinearMeshAnalisys();
void BrazilianPlasticAnalysis(TPZGeoMesh mesh);
void BlendMeshAnalisys();
void BrazilianElasticAnalysis();
//void CheckConvTests();
void SimpleMeshPlasticAnalysis();
void ComplexMeshPlasticAnalysis();
void SolveSistLin(TPZAnalysis &an, TPZCompMesh *fCmesh);
void CMeshBC(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);
void BrazilianElasticAnalysis2D();
void CMesh2DBC(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);
void SolveSistLin2D(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SetUPPostProcessElasticVariables2D(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void CheckConvTests(int a);
void CMesh2DSlopeStability(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);
void SlopeStability();
void BrazilianPlasticAnalysis2D();
void PressureCilinder();
void CMeshPressuredCilinder(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);
void Bean();
void CMeshBean(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);



//
//#define BOOST_TEST_DYN_LINK
//#include <boost/test/unit_test.hpp>
//
//
//BOOST_AUTO_TEST_SUITE(int_quad_2d_test);
//
//BOOST_AUTO_TEST_CASE(test_npoints)
//{
//	
//	int npoints = 1;
//	
//	BOOST_CHECK_EQUAL(npoints,1);
//	
//
//	npoints = 4;
//	
//	BOOST_CHECK_EQUAL(npoints,4);
//
//	npoints = 9;
//	
//	BOOST_CHECK_EQUAL(npoints,9);
//}
//

//#include "TPZPlasticityTest.h"

int main28(int argc, char * const argv[]) 
{
	
	InitializePZLOG("log4cxx.cfg");
	
	//	BrazilianElasticAnalysis2D();
	//  ComplexMeshPlasticAnalysis();
	//	BrazilianPlasticAnalysis2D();
	//  PressureCilinder();
	//	TPZVec< TPZTensor<REAL> > a;
	//	TPZPlasticTest::PlastifiedTensors(10,30.,120.,a);
	//	a.Print(cout);
	//	TPZPlasticTest::UndocumentedTest2();
	//	TPZPlasticTest::DruckerTest();
	 CheckConvTests(1);
	//	SlopeStability();
	//  Bean();
	
	//  TPZPlasticTest::VonMisesTest();
	//  TPZPlasticTest::DruckerPrager();
	
	//	int intervals =5;
	
	//	
	//	TPZDruckerPrager DP;
	//	DP.SetUp();
	//	
	//	TPZVonMises VM;
	//	VM.SetUp();
	//	
	//	TPZTensor<REAL> epst,epsp;
	//	
	//	epst.fData[_XX_] =  0.001;
	//	epst.fData[_XY_] =  0.0000001;
	//	epst.fData[_XZ_] =  0.0000001;
	//	epst.fData[_YY_] =  0.0000001;
	//	epst.fData[_YZ_] =  0.0000001;
	//	epst.fData[_ZZ_] =  0.0000001;
	//	
	//	epsp.fData[_XX_] =  0.;
	//	epsp.fData[_XY_] =  0.;
	//	epsp.fData[_XZ_] =  0.;
	//	epsp.fData[_YY_] =  0.;
	//	epsp.fData[_YZ_] =  0.;
	//	epsp.fData[_ZZ_] =  0.;
	//	
	//	
	//	
	//TPZPlasticTest::PlasticIntegratorCheck(intervals,VM);
	//TPZPlasticTest::PlasticIntegratorCheck(intervals,DP);
	
	//	TPZPlasticTest::VonMisesTest();
	//	TPZPlasticTest::DruckerTest();
	//TPZPlasticTest::MohrCoulombTest();
//	TPZPlasticTest::WillamWarnkeTest();
	
    return 0;
}

void Bean()
{
	
	TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
	
	int BC1,BC2,nsteps,taxa,nnewton;
	int order = 2;
	REAL tol = 1.e-5;
	nnewton = 10;
	BC1 = -3;
	nsteps = 10;
	taxa = 1;
	BeginForce(1,0) = 0.05;
	EndForce(1,0) = 0.07;
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	
    MESH = GeoMeshClass::Beangm();
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZVonMises * Pstep = new TPZVonMises();
	
	REAL Yield = 0.24;//GPA
	Pstep->fTFA.SetUp(Yield,1000.);
	Pstep->fER.SetUp(/*young GPA */ 210., /*poisson*/ 0.30);
	
	TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);//(ID,PlaneStrain)
	EPMat2.SetPlasticity(*Pstep);
	TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
	plastic->Print(cout);
	CMESH->InsertMaterialObject(plastic);
	
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
	
	CMeshBean(CMESH,plastic);	
	
	ofstream arg("CMESHPLASTIC2D.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	SolveSist(EPAnalysis,CMESH);
	
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"BEAN.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	//	ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	
	EPAnalysis.IterativeProcess(cout, tol, nnewton);
	cout << EPAnalysis.Solution() << endl;
	EPAnalysis.AcceptSolution();
	EPAnalysis.TransferSolution(PPAnalysis);
	PPAnalysis.PostProcess(0);
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"BEAN.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
}	

void CMeshBean(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	
	REAL BIG = 1.e12;
	
	TPZFMatrix kk(2,2,0.);
	TPZFMatrix ff(2,1,0.);
	kk(0,0)=BIG;
	kk(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -1,2, kk, ff);
	CMESH->InsertMaterialObject(ContBC);
	
	TPZFMatrix k1(2,2,0.);
	TPZFMatrix f1(2,1,0.);
	k1(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2,2, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
	
	TPZFMatrix k2(2,2,0.);
	TPZFMatrix f2(2,1,0.);
	f2(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -3, 5, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
	
	CMESH->AutoBuild();
}


void CheckConvTests(int a)
{
	
	
	
	//TPZPlasticTest::WillamWarnkeTest();
	//TPZPlasticTest::MohrCoulombTest();
	//	TPZPlasticTest::VonMisesTest();
	
	TPZYCWillamWarnke WW;
	TPZYCMohrCoulomb  MC;
//	TPZYCModifiedMohrCoulomb MMC;
	//	TPZYCVonMises  VM;
	TPZFNMatrix<6> input(6,1), Range(6,1);
	input(_XX_) = -2;
	input(_YY_) = -2.;
	input(_ZZ_) =  2.;
	input(_XY_) =  2.;
	input(_XZ_) =  2.;
	input(_YZ_) =  0.;
	Range = input * (1./19.);
	TPZVec< REAL > Coefs(1,1.);
//	CheckConvergence(MMC, input, Range, Coefs);
	
	//StressAtPoint.DruckerPragerTest();
}

//void RotateMesh(TPZGeoMesh &geomesh, REAL angleDegree,int rotateabout)
//{
//	REAL pi = M_PI;
//	REAL th = angleDegree/180. * pi;
//	
//	for(int node = 0; node < geomesh.NodeVec().NElements(); node++)
//	{
//		TPZFMatrix nodeCoord(3,1,0.);
//		for(int c = 0; c < 3; c++)
//		{
//			nodeCoord.Put(c, 0, geomesh.NodeVec()[node].Coord(c));
//		}
//		RotateMatrix(nodeCoord,th,rotateabout);
//		for(int c = 0; c < 3; c++)
//		{
//			geomesh.NodeVec()[node].SetCoord(c,nodeCoord(c,0));
//		}
//	}
//}

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	//TPZParSkylineStructMatrix full(fCmesh,2);
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver step;
	step.SetDirect(ELDLt);
	an.SetSolver(step);
	
}

void SolveSistCG(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	int numiter = 30;
	REAL tol = 1.e-5;
	//	TPZSpStructMatrix StrMatrix(an.Mesh());
	//	an.SetStructuralMatrix(StrMatrix);
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZMatrix * mat = full.Create();
	//	TPZMatrix * mat = StrMatrix.Create();
	TPZBlockDiagonalStructMatrix strBlockDiag(an.Mesh());
	TPZStepSolver Pre;
	TPZBlockDiagonal * block = new TPZBlockDiagonal();
	
    strBlockDiag.AssembleBlockDiagonal(*block); // just to initialize structure
	Pre.SetMatrix(block);
    Pre.SetDirect(ECholesky);
    TPZStepSolver Solver;
 	Solver.SetBiCGStab(numiter, Pre, tol, 0);
    Solver.SetMatrix(mat);
    an.SetSolver(Solver);
	
	an.BuildPreconditioner(TPZAnalysis::ENodeCentered, true);
	
}


void SolveSistLin(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	
	TPZParSkylineStructMatrix sky(fCmesh,2);
	an.SetStructuralMatrix(sky);
	TPZStepSolver step;
	step.SetDirect(ECholesky);
	an.SetSolver(step);
	an.SetCompMesh(fCmesh);
	//an.Run(cout);
	
}
void SetUPPostProcessElasticVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	vecnames.Resize(5);
	vecnames[0] = "Displacement";
	vecnames[1] = "NormalStress";
	vecnames[2] = "NormalStrain";
	vecnames[3] = "Strain";
	vecnames[4] = "Stress";
	postprocvars.Resize(scalnames.NElements()+vecnames.NElements());
	int i, k=0;
	for(i = 0; i < scalnames.NElements(); i++)
	{
		postprocvars[k] = scalnames[i];
		k++;
	}
	for(i = 0; i < vecnames.NElements(); i++)
	{
		postprocvars[k] = vecnames[i];
		k++;
	}
	
}

void SetUPPostProcessVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	scalnames.Resize(10);
	scalnames[0] = "Alpha";
	scalnames[1] = "PlasticSteps";
	scalnames[2] = "VolElasticStrain";
	scalnames[3] = "VolPlasticStrain";
	scalnames[4] = "VolTotalStrain";
	scalnames[5] = "I1Stress";
	scalnames[6] = "J2Stress";
	scalnames[7] = "YieldSurface";
	scalnames[8] = "EMisesStress";
	scalnames[9] = "TotalPlasticStrain";
	
	vecnames.Resize(5);
	vecnames[0] = "Displacement";
	vecnames[1] = "NormalStress";
	vecnames[2] = "ShearStress";
	vecnames[3] = "NormalStrain";
	vecnames[4] = "ShearStrain";
	//vecnames[5] = "NormalPlasticStrain";
	
	
	postprocvars.Resize(scalnames.NElements()+vecnames.NElements());
	int i, k=0;
	for(i = 0; i < scalnames.NElements(); i++)
	{
		postprocvars[k] = scalnames[i];
		k++;
	}
	for(i = 0; i < vecnames.NElements(); i++)
	{
		postprocvars[k] = vecnames[i];
		k++;
	}
}

void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
							int BCId,int BCId2, int nsteps, REAL PGRatio,
							TPZFMatrix & val1Begin, TPZFMatrix & val1End,
							TPZFMatrix & val2Begin, TPZFMatrix & val2End,
							TPZPostProcAnalysis * ppAnalysis, int res)
{
	if(!analysis.Mesh())return;
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() ***";
		sout << "\nWith parameters:\n";
		sout << "\ntol = " << tol;
		sout << "\nnumiter = " << numiter;
		sout << "\nBCId = " << BCId;
		sout << "\nBCId2 = " << BCId2;
		sout << "\nnsteps = " << nsteps;
		sout << "\nPGRatio = " << PGRatio;
		sout << "\nval1Begin = " << val1Begin;
		sout << "\nval1End = " << val1End;
		sout << "\nval2Begin = " << val2Begin;
		sout << "\nval2End = " << val2End;
		if(ppAnalysis)
		{
			sout << "\nppanalysis set";
		}else
		{
			sout << "\nppanalysis NOT set";
		}
		//LOGPZ_INFO(logger,sout.str().c_str());
	}
#endif
	
	// computing the initial value for the PG progression such that its sum equals one;
	REAL a0;
	
	if(fabs(PGRatio - 1.) < 1.e-3)
	{
	    a0 = 1. / REAL(nsteps);
	}
	
	else
	{
		a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
	}
	TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
	TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
	
	deltaVal1 = val1End;
	deltaVal1.ZAXPY(-1., val1Begin);
	deltaVal2 = val2End;
	deltaVal2.ZAXPY(-1., val2Begin);	
	
	
	//-19
	TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
	if(!pBC)return;
	
	//-20
	//COMENTAR NO CILINDRO
	//	TPZAutoPointer<TPZMaterial> mat2 = analysis.Mesh()->FindMaterial(BCId2);
	//	TPZBndCond * pBC2 = dynamic_cast<TPZBndCond *>(mat2.operator->());
	//	if(!pBC2)return;
	
	
    int i;
	for(i = 0; i < nsteps; i++)
	{
		REAL stepLen;
		if(fabs(PGRatio - 1.) < 1.e-3)
		{
			stepLen = REAL(i+1) / REAL(nsteps);
		}
		else
		{
		    stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
		}
		
		val1 = val1Begin;
		val1.ZAXPY(stepLen, deltaVal1);
		val2 = val2Begin;
		val2.ZAXPY(stepLen, deltaVal2);
		
		pBC->Val1() = val1;
		pBC->Val2() = val2;
		
		cout <<  "PRESSURE = " << val2 << endl;
		
		//COMENTAR NO CILINDRO
		//	pBC2->Val1() = -1.*val1;
		//	pBC2->Val2() = -1.*val2;
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i;
			sout << "\n stepLen = " << stepLen;
			sout << "\n PBC1 ";
			pBC->Print(sout);
			//sout << "\n PBC2 ";
			//pBC2->Print(sout);
			//LOGPZ_INFO(logger2,sout.str().c_str());
		}
#endif
		
		
		analysis.IterativeProcess(out, tol, numiter);
		
		
		
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** load step " << i << " ended";
			//LOGPZ_INFO(logger2,sout.str().c_str());
		}
#endif
		
		analysis.AcceptSolution();
		
		
		if(ppAnalysis)
		{
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "*** TPZElastoPlasticAnalysis::ManageIterativeProcess() *** PostProcessing ";
				//	LOGPZ_INFO(logger2,sout.str().c_str());
			}
#endif
			analysis.TransferSolution(*ppAnalysis);
			ppAnalysis->PostProcess(res);
		}
	}
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() *** Exiting";
		//	LOGPZ_INFO(logger2,sout.str().c_str());
	}
#endif
}
/*
 #include <math.h>
 void RotationMatrix(TPZFMatrix &R, double thetaRad, int axis)
 {
 R.Resize(3,3);
 
 switch (axis) 
 {
 
 case 0://ROTATE ABOUT X
 
 R.Put(0,0,1.);
 R.Put(1,1,cos(thetaRad));R.Put(1,2,sin(thetaRad));
 R.Put(2,1,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
 
 break;
 
 case 1://ROTATE ABOUT Y
 
 R.Put(1,1,1.);
 R.Put(0,0,cos(thetaRad));R.Put(0,2,sin(thetaRad));
 R.Put(2,0,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
 
 break;
 
 case 2://ROTATE ABOUT Z
 
 R.Put(0,0,cos(thetaRad)); R.Put(0, 1,sin(thetaRad));
 R.Put(1,0,-sin(thetaRad)); R.Put(1, 1,cos(thetaRad));
 R.Put(2,2,1.);
 
 break;
 
 default:
 
 std::cout << " NON SPECIFIED AXIS "<<std::endl;
 break;
 }
 
 }
 
 void RotateMatrix(TPZFMatrix &Mat, double thetaRad,int rotateaboutaxes)
 {
 TPZFMatrix R;
 TPZFMatrix temp;
 RotationMatrix(R, thetaRad,rotateaboutaxes);
 if(R.Cols() != Mat.Rows())
 {
 cout << " \n -- MATRIX WITH INCOMPATIBLE SHAPES -- \n";
 DebugStop();
 }
 #ifdef LOG4CXX
 {
 std::stringstream sout;
 sout << "\n" <<std::endl;
 sout << "\n ROTATION MATRIX "<< R << std::endl;
 sout << "\n MATRIX TO BE ROTATED " << Mat;
 //	LOGPZ_DEBUG(logger,sout.str().c_str());
 }
 #endif
 
 int matcol = Mat.Cols();
 TPZFMatrix RT;
 if(matcol == 1)
 {
 R.Transpose(&RT);
 temp = RT*Mat;
 Mat = temp;
 }
 else
 {
 R.Transpose(&RT);
 temp = RT*Mat;
 Mat = temp*R;
 }
 //	if(matcol == 1)
 //	{
 //		temp = R*Mat;
 //		Mat = temp;
 //	}
 //	else
 //	{
 //		temp = R*Mat;
 //		R.Transpose();
 //		Mat = temp*R;
 //	}
 #ifdef LOG4CXX
 {
 std::stringstream sout;
 sout << "\n" <<std::endl;
 sout << "\n ROTATED MATRIX " << Mat;
 //LOGPZ_DEBUG(logger,sout.str().c_str());
 }
 #endif
 
 }
 
 */

void BrazilianElasticAnalysis()
{
	int h=4;
	int order = 2;
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	BrazilianTestGeoMesh::TransformBlendToLinearMesh(MESH,h);
	ofstream arg1("ElasticMaterialGeoMesh.txt");
	MESH->Print(arg1);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZVec<REAL> force(3,1);
	force[0]=0.;
	force[1]=0.;
	force[2]=0.;
	
	TPZElasticity3D *elast = new TPZElasticity3D(1,20000.,0.,force);
	
	CMESH->SetDefaultOrder(order);
	TPZAutoPointer<TPZMaterial> elastt = elast;
	CMESH->InsertMaterialObject(elastt);
	
	CMeshBC(CMESH,elastt);
	
	ofstream arg("CMESHELASTIC.txt");
	CMESH->Print(arg);
	
	TPZAnalysis an;
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	
	SetUPPostProcessElasticVariables(PostProcVars,scalNames,vecNames);
	
	SolveSistLin(an,CMESH);
	
	cout << " \n ELASTIC-SOL \n";
	cout << an.Solution() << endl;
	an.DefineGraphMesh(3, scalNames,vecNames,"ELASTIC1.vtk");
	an.PostProcess(0);
	an.CloseGraphMesh();
	
}

void ComplexMeshPlasticAnalysis()
{
	TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
	TPZVec<int> TopolPoint(1);
	int h=2;
	int order = 2;
	TPZGeoMesh * MESH = new TPZGeoMesh;
	BrazilianTestGeoMesh::TransformBlendToLinearMesh2(MESH,h);
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZDruckerPrager * Pstep = new TPZDruckerPrager();
	Pstep->fYC.SetUp(/*phi=20*/ 20./180. * M_PI ,/*innerMCFit*/0);
	Pstep->fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 9.2376, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep->fER.SetUp(/*young*/ 20000., /*poisson*/ 0.);
	typedef TPZMatElastoPlastic <TPZDruckerPrager> EPMat;
	EPMat *phil = new EPMat(1);
	phil->SetPlasticity(*Pstep);
	CMESH->SetDefaultOrder(order);
	TPZAutoPointer<TPZMaterial> plastic = phil;
	CMESH->InsertMaterialObject(plastic);
	
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
	
	CMeshBC(CMESH,plastic);	
	
	ofstream arg("CMESHPLASTIC.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	SolveSist(EPAnalysis,CMESH);
	//	SolveSistCG(EPAnalysis,CMESH);
	
	int dimension =3;
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	//	TPZParSkylineStructMatrix structmatrix(PPAnalysis.Mesh(),2);
	TPZSkylineStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"3DH2P2.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	EPAnalysis.IterativeProcess(cout, 1.e-5, 30);
	cout << " \n PPPLASTIC-SOL \n";
	cout << EPAnalysis.Solution() << endl;
	EPAnalysis.AcceptSolution();
	EPAnalysis.TransferSolution(PPAnalysis);
	PPAnalysis.PostProcess(0);
	
	//	BeginForce(0,0) = 0.;
	//	EndForce(0,0) = 70.;
	//	ManageIterativeProcess(EPAnalysis,cout, 1.e-5, 30,
	//						   -5/*BCId*/,-6/*BCId2*/, 1 /*nsteps*/, 1/*PGRatio*/,
	//						   BeginStress/*val1Begin*/, EndStress/*val1End*/,
	//						   BeginForce/*val2Begin*/, EndForce/*val2End*/,
	//						   &PPAnalysis, 0);
	
	//ManageIterativeProcess(EPAnalysis,cout,1.e-5,30,-5,-6,10,1,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"3DH2P2.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
	
}

void CMeshBC(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	REAL BIG = 1.e-2;
	TPZFMatrix k1(3,3,0.);
	TPZFMatrix f1(3,1,0.);
	f1(0,0) = 50.;
	TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -5, 5, k1, f1);
	CMESH->InsertMaterialObject(ContBC);
	
	
	TPZFMatrix k2(3,3,0.);
	TPZFMatrix f2(3,1,0.);
	f2(0,0) = 50.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -6, 5, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
	
	
	////// MIXED CONDITION /////
	
	//CONDICAO MISTA
	TPZFMatrix k3(3,3,0.);
	TPZFMatrix f3(3,1,0.);
	//O BIG VAI NA DIRECAO QUE SE DESEJA IMPEDIR NO CASO EM X, Y e Z
	k3(0,0)=BIG,k3(1,1)=BIG,k3(2,2)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -4, 2, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
	
	//CONDICAO MISTA
	TPZFMatrix k4(3,3,0.);
	TPZFMatrix f4(3,1,0.);
	k4(1,1)=BIG,k4(2,2)=BIG;//O BIG VAI NA DIRECAO QUE SE DESEJA IMPEDIR NO CASO EM Y E Z
	TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -3, 2, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
	
	//CONDICAO MISTA
	TPZFMatrix k5(3,3,0.);
	TPZFMatrix f5(3,1,0.);
	k5(1,1)=BIG;
	//O BIG VAI NA DIRECAO QUE SE DESEJA IMPEDIR NO CASO EM Z
	TPZAutoPointer<TPZMaterial> ContBC5 = mat->CreateBC(mat, -2, 2, k5, f5);
	CMESH->InsertMaterialObject(ContBC5);
	
	
	/*	
	 //Z=0
	 TPZFMatrix k6(3,3,0.);
	 TPZFMatrix f6(3,1,0.);
	 f6(0,0) = 5.;
	 TPZAutoPointer<TPZMaterial> ContBC6 = mat->CreateBC(mat, -7, 5, k6, f6);
	 CMESH->InsertMaterialObject(ContBC6);
	 
	 //Z=30
	 TPZFMatrix k7(3,3,0.);
	 TPZFMatrix f7(3,1,0.);
	 f7(0,0) = 5.;
	 TPZAutoPointer<TPZMaterial> ContBC7 = mat->CreateBC(mat, -8, 5, k7, f7);
	 CMESH->InsertMaterialObject(ContBC7);
	 
	 //LATERAIS
	 TPZFMatrix k8(3,3,0.);
	 TPZFMatrix f8(3,1,0.);
	 f8(0,0) = 5.;
	 TPZAutoPointer<TPZMaterial> ContBC8 = mat->CreateBC(mat, -9, 5, k8, f8);
	 CMESH->InsertMaterialObject(ContBC8);	
	 */
	
	CMESH->AutoBuild();
	
}

void BrazilianPlasticAnalysis2D()
{
	TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
	
	int BC1,BC2,nsteps,taxa,nnewton;
	int h=3;
	int order = 2;
	REAL tol = 1.e-3;
	nnewton = 30;
	BC1 = -2;//down line
	BC2 = -3;//upper line
	nsteps = 10;
	taxa = 1;
	BeginForce(1,0) = 30.;
	EndForce(1,0) = 40.;
	
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	//void BrazilianTestGeoMesh::TransformBlendToLinearMesh2(TPZGeoMesh *newlinearmesh, int h)
	//BrazilianTestGeoMesh::TransformBlendToLinearMesh2(MESH,h);
	int refdir = 5;
	MESH = BrazilianTestGeoMesh::TwoDMesh(h,refdir);
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZDruckerPrager * Pstep = new TPZDruckerPrager();
	//	TPZMohrCoulomb * Pstep = new TPZMohrCoulomb();
	Pstep->fYC.SetUp(/*phi=20*/ 20./180. * M_PI ,/*innerMCFit*/0);
	//	Pstep->fYC.SetUp(/*phi=20*/ 20./180. * M_PI);
	Pstep->fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 9.2376, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep->fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);
	
	TPZMatElastoPlastic2D<TPZDruckerPrager> EPMat2(1,1);
	EPMat2.SetPlasticity(*Pstep);
	TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
	
	plastic->Print(cout);
	CMESH->InsertMaterialObject(plastic);
	
	
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
	
	CMesh2DBC(CMESH,plastic);	
	
	ofstream arg("CMESHPLASTIC2D.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	SolveSist(EPAnalysis,CMESH);
	
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"Br_Test3085H3P2.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	
	ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	
	//	cout << "\nInitial Solution Exported. Solving Problem\n";
	//	EPAnalysis.IterativeProcess(cout, 1.e-5, 30);
	//	cout << " \n PPPLASTIC-SOL \n";
	//	cout << EPAnalysis.Solution() << endl;
	//	EPAnalysis.AcceptSolution();
	//	EPAnalysis.TransferSolution(PPAnalysis);
	//	PPAnalysis.PostProcess(0);
	//	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"Br_Test3085H2P2.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
}

void CMesh2DBC(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	
	
	int neumannPressure = 5;
	int neumann= 1;
	int dirichelet = 0;
	int mixed =2;
	REAL BIG = 1.e-2;
	
	
	//down line
	TPZFMatrix k1(2,2,0.);
	TPZFMatrix f1(2,1,0.);
	f1(1,0) = 17.;
	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2, neumann, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
	
	//upper line
	TPZFMatrix k2(2,2,0.);
	TPZFMatrix f2(2,1,0.);
	f2(1,0) = -17.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -3, neumann, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
	
	//CONDICAO MISTA
	TPZFMatrix k3(2,2,0.);
	TPZFMatrix f3(2,1,0.);
	
	k3(0,0)=BIG,k3(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -4, mixed, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
	
	//CONDICAO MISTA
	TPZFMatrix k4(3,3,0.);
	TPZFMatrix f4(3,1,0.);
	k4(1,1)=BIG;//,k4(2,2)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -5, mixed, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
	
	
	CMESH->AutoBuild();
}

void SolveSistLin2D(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	TPZSkylineStructMatrix sky(fCmesh);
	an.SetStructuralMatrix(sky);
	TPZStepSolver step;
	step.SetDirect(ECholesky);
	an.SetSolver(step);
	an.SetCompMesh(fCmesh);
	an.Run(cout);
}

void SetUPPostProcessElasticVariables2D(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	vecnames.Resize(1);
	vecnames[0] = "displacement";
	
	scalnames.Resize(4);
	scalnames[0]="SigmaX";
	scalnames[1]="SigmaY";
	scalnames[2]="tau_xy";
	scalnames[3]="Pressure";
	postprocvars.Resize(scalnames.NElements()+vecnames.NElements());
	int i, k=0;
	for(i = 0; i < scalnames.NElements(); i++)
	{
		postprocvars[k] = scalnames[i];
		k++;
	}
	for(i = 0; i < vecnames.NElements(); i++)
	{
		postprocvars[k] = vecnames[i];
		k++;
	}
	
}

void SlopeStability()
{
	
	
	int nnewton;
	int h=1;
	int order = 2;
	REAL tol = 1.e-5;
	nnewton = 30;
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	
	int refdir = 3;
	MESH = BrazilianTestGeoMesh::TwoDMeshSlopeStability45(h,refdir);
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZDruckerPrager * Pstep = new TPZDruckerPrager();
	//		TPZMohrCoulomb * Pstep = new TPZMohrCoulomb();
	Pstep->fYC.SetUp(/*phi=20*/ 20./180. * M_PI ,/*innerMCFit*/0);
	//		Pstep->fYC.SetUp(/*phi=20*/ 20./180. * M_PI);
	Pstep->fTFA.SetUp(/*Kpa */ 50.,1.);
	//		Pstep->fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 50., /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
	Pstep->fER.SetUp(/*young*/ 20000., /*poisson*/ 0.49);
	
	TPZMatElastoPlastic2D<TPZDruckerPrager> EPMat2(1,1);
	EPMat2.SetBulkDensity(6.);
	EPMat2.SetPlasticity(*Pstep);
	TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
	plastic->Print(cout);
	CMESH->InsertMaterialObject(plastic);
	
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
	
	CMesh2DSlopeStability(CMESH,plastic);	
	
	ofstream arg("CMESHPLASTIC2D.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	//	SolveSist(EPAnalysis,CMESH);
	SolveSistLin(EPAnalysis,CMESH);	
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"Nslop9H3.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	//	ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	EPAnalysis.IterativeProcess(cout, tol, nnewton);
	//	cout << EPAnalysis.Solution() << endl;
	EPAnalysis.AcceptSolution();
	EPAnalysis.TransferSolution(PPAnalysis);
	PPAnalysis.PostProcess(0);
	
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"Nslop9H3.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
}

void CMesh2DSlopeStability(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	
	
	int dirichelet = 0;
	//int mixed =2;
	//int neumann =1;
	REAL BIG = 1.e15;
	
	TPZFMatrix k(2,2,0.);
	TPZFMatrix f(2,1,0.);
	//	k(0,0)=BIG,k(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -1, dirichelet, k, f);
	CMESH->InsertMaterialObject(ContBC);
	
	TPZFMatrix k2(2,2,0.);
	TPZFMatrix f2(2,1,0.);
	//k2(0,0)=BIG;
	f2(0,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -2, 3, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
	
	
	TPZFMatrix k3(2,2,0.);
	TPZFMatrix f3(2,1,0.);
	//O BIG VAI NA DIRECAO QUE SE DESEJA IMPEDIR NO CASO EM X, Y e Z
	//k3(0,0)=BIG;
	f3(0,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -3, 3, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
	
	
	CMESH->AutoBuild();
}

void PressureCilinder()
{
	
	TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
	
	int BC1,BC2,nsteps,taxa,nnewton;
	int h=3;
	int order = 2;
	REAL tol = 1.e-5;
	nnewton = 10;
	BC1 = -3;
	nsteps =10;
	taxa = 1;
	BeginForce(0,0) = 0.;
	EndForce(0,0) = 0.25;
	
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	
	int refdir = 3;
	MESH = BrazilianTestGeoMesh::MisesPressure(h,refdir);
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	//TPZDruckerPrager * Pstep = new TPZDruckerPrager();
	TPZVonMises * Pstep = new TPZVonMises();
	
	REAL Yield = 0.24;//GPA
	Pstep->fTFA.SetUp(Yield,1000.);
	Pstep->fER.SetUp(/*young GPA */ 210., /*poisson*/ 0.30);
	
	TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);
	EPMat2.SetPlasticity(*Pstep);
	TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
	plastic->Print(cout);
	CMESH->InsertMaterialObject(plastic);
	
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
	
	CMeshPressuredCilinder(CMESH,plastic);	
	
	ofstream arg("CMESHPLASTIC2D.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	SolveSist(EPAnalysis,CMESH);
	
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	//		EPAnalysis.IterativeProcess(cout, tol, nnewton);
	//		cout << EPAnalysis.Solution() << endl;
	//		EPAnalysis.AcceptSolution();
	//		EPAnalysis.TransferSolution(PPAnalysis);
	//		PPAnalysis.PostProcess(0);
	
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess3.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
	//	TPZFMatrix BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	//	TPZFMatrix val1(3,1,0.);TPZFMatrix val2(3,1,0.);TPZFMatrix BeginForce(3,1,0.);TPZFMatrix EndForce(3,1,0.);
	//	
	//	int BC1,BC2,nsteps,taxa,nnewton;
	//	int h=3;
	//	int order = 2;
	//	REAL tol = 1.e-5;
	//	nnewton = 10;
	//	BC1 = -3;
	//	nsteps =2;
	//	taxa = 1;
	//	BeginForce(0,0) = 0.;
	//	EndForce(0,0) = 0.2;
	//	
	//	
	//	TPZGeoMesh * MESH = new TPZGeoMesh;
	//	
	//	int refdir = 3;
	//	MESH = BrazilianTestGeoMesh::MisesPressure(h,refdir);
	//	
	//	ofstream arg1("GeoMesh.txt");
	//	MESH->Print(arg1);
	//	TPZCompEl::SetgOrder(order);
	//	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	//	
	//	//TPZDruckerPrager * Pstep = new TPZDruckerPrager();
	//	TPZVonMises * Pstep = new TPZVonMises();
	//	
	//	REAL Yield = 0.24;//GPA
	//	Pstep->fTFA.SetUp(Yield,1.);
	//	Pstep->fER.SetUp(/*young GPA */ 210., /*poisson*/ 0.30);
	//	
	//	TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);
	//	EPMat2.SetPlasticity(*Pstep);
	//	TPZAutoPointer<TPZMaterial> plastic(&EPMat2);
	//	//plastic->SetBulkDensity(20.);
	//	plastic->Print(cout);
	//	CMESH->InsertMaterialObject(plastic);
	//	
	//	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem();
	//	
	//	CMeshPressuredCilinder(CMESH,plastic);	
	//	
	//	ofstream arg("CMESHPLASTIC2D.txt");
	//	CMESH->Print(arg);
	//	
	//	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	//	
	//	SolveSist(EPAnalysis,CMESH);
	//	
	//	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	//	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	//	PPAnalysis.SetStructuralMatrix(structmatrix);
	//	TPZVec<int> PostProcMatIds(1,1);
	//	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	//	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
	//	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	//	
	//	EPAnalysis.TransferSolution(PPAnalysis);
	//	
	//	cout << "\nDefining Graph Mesh\n";
	//	int dimension =2;
	//	
	//	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess.vtk");
	//	
	//	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	//	
	//	PPAnalysis.PostProcess(0/*pOrder*/);
	//	
	//	ManageIterativeProcess(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);
	//	
	//	cout << "\nInitial Solution Exported. Solving Problem\n";
	//	//	EPAnalysis.IterativeProcess(cout, tol, nnewton);
	//	//	cout << EPAnalysis.Solution() << endl;
	//	//	EPAnalysis.AcceptSolution();
	//	//	EPAnalysis.TransferSolution(PPAnalysis);
	//	//	PPAnalysis.PostProcess(0);
	//	
	//	
	//	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"TuboNewPostProcess.vtk");
	//	PPAnalysis.PostProcess(0);
	//	PPAnalysis.CloseGraphMesh();
	
}

void CMeshPressuredCilinder(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	
	
	//	int dirichelet = 0;
	//	int mixed =2;
	//	int neumann =1;
	//int pressure =5;
	REAL BIG = 1.e12;
	
	
	TPZFMatrix kk(2,2,0.);
	TPZFMatrix ff(2,1,0.);
	ff(0,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -3, 5, kk, ff);
	CMESH->InsertMaterialObject(ContBC);
	
	TPZFMatrix k1(2,2,0.);
	TPZFMatrix f1(2,1,0.);
	//		k1(0,0)=BIG;
	//		k1(1,1)=BIG;
	f1(0,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2,3, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
	
	
	TPZFMatrix k2(2,2,0.);
	TPZFMatrix f2(2,1,0.);
	//		k2(0,0)=BIG;
	//		k2(1,1)=BIG;
	f2(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -1, 3, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
	
	
	
	CMESH->AutoBuild();
	
	//	//	int dirichelet = 0;
	//	//	int mixed =2;
	//	//	int neumann =1;
	//	int pressure =5;
	//	REAL BIG = 1.e12;
	//	
	//	
	//	TPZFMatrix k(2,2,0.);
	//	TPZFMatrix f(2,1,0.);
	//	f(0,0)=1.;
	//	TPZAutoPointer<TPZMaterial> ContBC = mat->CreateBC(mat, -3, pressure, k, f);
	//	CMESH->InsertMaterialObject(ContBC);
	//	
	//	TPZFMatrix k1(2,2,0.);
	//	TPZFMatrix f1(2,1,0.);
	//	//	k1(0,0)=BIG;
	//	//	k1(1,1)=BIG;
	//	f1(0,0)=1.;
	//	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2, 3, k1, f1);
	//	CMESH->InsertMaterialObject(ContBC1);
	//	
	//	
	//	TPZFMatrix k2(2,2,0.);
	//	TPZFMatrix f2(2,1,0.);
	//	//	k2(0,0)=BIG;
	//	//	k2(1,1)=BIG;
	//	f2(1,0)=1.;
	//	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -1, 3, k2, f2);
	//	CMESH->InsertMaterialObject(ContBC2);
	//	
	//	
	//	
	//	/*	
	//	 TPZFMatrix k3(2,2,0.);
	//	 TPZFMatrix f3(2,1,0.);
	//	 k3(0,0)=BIG;
	//	 k3(1,1)=BIG;
	//	 TPZAutoPointer<TPZMaterial> ContBC3= mat->CreateBC(mat, -9, mixed, k3, f3);
	//	 CMESH->InsertMaterialObject(ContBC3);
	//	 
	//	 TPZFMatrix k4(2,2,0.);
	//	 TPZFMatrix f4(2,1,0.);
	//	 k4(0,0)=BIG;
	//	 k4(1,1)=BIG;
	//	 TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -10, mixed, k4, f4);
	//	 CMESH->InsertMaterialObject(ContBC4);
	//	 */
	//	
	//	
	//	
	//	
	//	CMESH->AutoBuild();
	//	
}





