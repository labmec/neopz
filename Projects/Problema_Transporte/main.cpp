
#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>

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
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzl2projection.h"

#include "pzconvectionproblem.h"
#include "pzgradientreconstruction.h"

#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

using namespace std;

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material"));
#endif

const int matId = 1;
const int inflow = 0;
const int outflow = 3;

int matIdL2Proj = 2;

int neumann = 1;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;

TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly);
TPZGeoMesh *MalhaGeom2();
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder,TPZMatConvectionProblem * &material);
TPZCompMesh *SetCondicaoInicial(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void CreatInterface(TPZCompMesh *cmesh);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix);
void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh, TPZGradientReconstruction *gradreconst);
void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh);

void PosProcessSolution(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile);
TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZMatConvectionProblem * mymaterial, TPZCompMesh* cmesh);
void StiffMatrixLoadVec(TPZMatConvectionProblem *mymaterial, TPZCompMesh*cmesh, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<STATE> > &matK1, TPZFMatrix<STATE> &fvec);

void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

void FilterEquation(TPZMatConvectionProblem *mymaterial, TPZCompMesh *cmesh, TPZAnalysis &an, bool currentstate);

int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    int p = 1;
    TPZGeoMesh *gmesh = MalhaGeom(1.,1.);
    UniformRefine(gmesh,2);
    ofstream arg1("gmesh.txt");
    gmesh->Print(arg1);
    
    TPZMatConvectionProblem *material;
    TPZCompMesh *cmesh = MalhaComp(gmesh, p, material);
   // CreatInterface(cmesh);
    ofstream arg2("cmesh_inicial.txt");
    cmesh->Print(arg2);
    
    //Refinar elemnto: malha adaptada
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,3, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    //Refinar elemnto: malha adaptada
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,12, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,2, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,13, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,7, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,8, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,6, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompEl(cmesh,9, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    
//    gmesh->ResetReference();
//    cmesh->LoadReferences();
//    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh,2, false);
//    gmesh->ResetConnectivities();
//	gmesh->BuildConnectivity();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();

    
    //Set initial conditions
    TPZAnalysis an(cmesh);
	int nrs = an.Solution().Rows();
    TPZVec<STATE> solini(nrs,0.);
    TPZCompMesh  * cmeshL2 = SetCondicaoInicial(gmesh, p, solini);
    
    TPZAnalysis anL2(cmeshL2);
    ResolverSistema(anL2, cmeshL2, true);
    an.LoadSolution(anL2.Solution());
    //an.Solution().Print("sol_S0");
    
    REAL deltaT=0.01; //second
    material->SetTimeStep(deltaT);
    REAL maxTime = 1.;
    material->SetTimeValue(maxTime);
    
    gmesh->ResetReference();
    cmesh->LoadReferences();
    CreatInterface(cmesh);
    ofstream arg3("gmesh_final.txt");
    gmesh->Print(arg3);
    
    ofstream arg4("cmesh_final.txt");
    cmesh->Print(arg4);
    
    TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(false,1.);
    SolveSistTransient(deltaT, material, cmesh, gradreconst);
}


TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly)
{
	
	int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <long> TopolQuad(4);
	TPZVec <long> TopolLine(2);
	
	//indice dos nos
	long id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id = 0;
    
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
	id++;
    
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    
	gmesh->BuildConnectivity();
    
	return gmesh;
	
}

TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh, int pOrder,TPZMatConvectionProblem * &material)
{
	/// criar materiais
	int dim = 2;
	
	material = new TPZMatConvectionProblem(matId,dim);
	TPZMaterial * mat(material);
	
	TPZVec<REAL> convdir(dim,0.);
    convdir[0] = 1.;
	REAL flux = 0.;
    REAL rho = 1.;
	
	material->SetParameters(rho,convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
	cmesh->InsertMaterialObject(mat);
    
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,neumann, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,neumann, val1, val2);
    REAL uD =1.;
    val2(0,0) = uD;
	TPZMaterial * BCond3 = material->CreateBC(mat, bc3,inflow, val1, val2);
    val2(0,0) = 1.;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,outflow, val1, val2);
    	
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond2);
    
    
    TPZVec<STATE> sol(1,0.);
    TPZL2Projection *matl2proj = new TPZL2Projection(matIdL2Proj,dim,material->NStateVariables(),sol);
    cmesh->InsertMaterialObject(matl2proj);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    //TPZCompElDisc::SetTotalOrderShape(cmesh);
	
	return cmesh;
}

#include "pzl2projection.h"
TPZCompMesh *SetCondicaoInicial(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
	int dim = 2;
	TPZL2Projection *material;
	material = new TPZL2Projection(matId, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(ForcingInicial);
    material->SetForcingFunction(forcef);
    
    
	cmesh->SetAllCreateFunctionsDiscontinuous();
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    
    //TPZCompElDisc::SetTotalOrderShape(cmesh);
    
	return cmesh;
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
	// Re-constructing connectivities
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}

void CreatInterface(TPZCompMesh *cmesh){
    
    for(long el = 0; el < cmesh->ElementVec().NElements(); el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        long index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension() || compEl->Dimension() == cmesh->Dimension()-1)
        {
            TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(cmesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces(false);
        }
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
}

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh, bool symmetric_matrix)
{
    if(symmetric_matrix ==true){
        TPZSkylineStructMatrix skmat(fCmesh);
        an.SetStructuralMatrix(skmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELDLt);
        an.SetSolver(direct);
    }
    else{
        TPZBandStructMatrix bdmat(fCmesh);
        an.SetStructuralMatrix(bdmat);
        TPZStepSolver<STATE> direct;
        direct.SetDirect(ELU);
        an.SetSolver(direct);
    }
	an.Run();
	
	//Saida de Dados: solucao txt
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);
}

void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh, TPZGradientReconstruction *gradreconst){
	
    TPZAnalysis an(cmesh);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "SolutionComRecGrad";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    //gradient reconstruction
    TPZFMatrix<REAL> datagradients;
    gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
    Initialsolution = cmesh->Solution();
    
    PosProcessSolution(cmesh,an,plotfile);

    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mymaterial, cmesh);
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
            std::stringstream sout;
        	matM->Print("matM = ", sout,EMathematicaInput);
        	LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvec;
    //StiffMatrixLoadVec(mymaterial, cmesh, an, matK, fvec);
    FilterEquation(mymaterial, cmesh, an, true);
    matK = an.Solver().Matrix();
    fvec = an.Rhs();
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{

        std::stringstream sout;
        matK->Print("matK = ", sout,EMathematicaInput);
        fvec.Print("fvec = ", sout,EMathematicaInput);
        //Print the temporal solution
        Initialsolution.Print("Initial conditions = ", sout,EMathematicaInput);
        TPZFMatrix<STATE> Temp;
        matM->Multiply(Initialsolution,Temp);
        Temp.Print("Temp matM = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
    
	long nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
    REAL maxTime;
    mymaterial->GetTimeValue(maxTime);
	while (TimeValue <= maxTime)
	{
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
#ifdef LOG4CXX
        if(logdata->isDebugEnabled())
        {
            std::stringstream sout;
            sout<< " tempo = " << cent;
            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logdata,sout.str())
        }
#endif
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
        Lastsolution = an.Solution();
        
        
        gradreconst-> ProjectionL2GradientReconstructed(cmesh, matIdL2Proj);
        an.LoadSolution(cmesh->Solution());
        Lastsolution = an.Solution();
        
        if(cent%1==0){
            std::stringstream outputfiletemp;
            outputfiletemp << outputfile << ".vtk";
            std::string plotfile = outputfiletemp.str();
            PosProcessSolution(cmesh,an,plotfile);
        }
        
        cent++;
		TimeValue = cent*deltaT;
	}
 
}

void SolveSistTransient(REAL deltaT, TPZMatConvectionProblem * &mymaterial, TPZCompMesh* cmesh){
	
    
    TPZAnalysis an(cmesh);
	TPZFMatrix<STATE> Initialsolution = an.Solution();
    
    std::string outputfile;
	outputfile = "SolutionMultiphysics";
    
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    
    PosProcessSolution(cmesh,an,plotfile);
    
    //Criando matriz de massa (matM)
    TPZAutoPointer <TPZMatrix<STATE> > matM = MassMatrix(mymaterial, cmesh);
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        matM->Print("matM = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
    
    //Criando matriz de rigidez (matK) e vetor de carga
	TPZAutoPointer< TPZMatrix<STATE> > matK;
	TPZFMatrix<STATE> fvec;
    StiffMatrixLoadVec(mymaterial, cmesh, an, matK, fvec);
    
//#ifdef LOG4CXX
//	if(logdata->isDebugEnabled())
//	{
//        
//        std::stringstream sout;
//        matK->Print("matK = ", sout,EMathematicaInput);
//        fvec.Print("fvec = ", sout,EMathematicaInput);
//        //Print the temporal solution
//        Initialsolution.Print("Intial conditions = ", sout,EMathematicaInput);
//        TPZFMatrix<STATE> Temp;
//        matM->Multiply(Initialsolution,Temp);
//        Temp.Print("Temp matM = ", sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str())
//	}
//#endif
    
    
	long nrows;
	nrows = matM->Rows();
	TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
	TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
	TPZFMatrix<STATE> Lastsolution = Initialsolution;
	
	REAL TimeValue = 0.0;
	int cent = 1;
	TimeValue = cent*deltaT;
    REAL maxTime;
    mymaterial->GetTimeValue(maxTime);
	while (TimeValue <= maxTime)
	{
		// This time solution i for Transient Analytic Solution
		mymaterial->SetTimeValue(TimeValue);
		matM->Multiply(Lastsolution,TotalRhstemp);
        
//#ifdef LOG4CXX
//        if(logdata->isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout<< " tempo = " << cent;
//            Lastsolution.Print("\nIntial conditions = ", sout,EMathematicaInput);
//            TotalRhstemp.Print("Mat Mass x Last solution = ", sout,EMathematicaInput);
//            LOGPZ_DEBUG(logdata,sout.str())
//        }
//#endif
        
		TotalRhs = fvec + TotalRhstemp;
		an.Rhs() = TotalRhs;
		an.Solve();
        Lastsolution = an.Solution();
        
        if(cent%10==0){
                    std::stringstream outputfiletemp;
                    outputfiletemp << outputfile << ".vtk";
                    std::string plotfile = outputfiletemp.str();
                    PosProcessSolution(cmesh,an,plotfile);
        }
        
        cent++;
		TimeValue = cent*deltaT;
	}
}


void PosProcessSolution(TPZCompMesh* cmesh, TPZAnalysis &an, std::string plotfile)
{
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "Solution";
	scalnames[1] = "ConvDirGradU";
   
	vecnames[0]  = "Derivative";
	vecnames[1]  = "Flux";
		
	const int dim = 2;
	int div =0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

#include "TPZSkylineNSymStructMatrix.h"
TPZAutoPointer <TPZMatrix<STATE> > MassMatrix(TPZMatConvectionProblem * mymaterial, TPZCompMesh* cmesh){
    
    mymaterial->SetLastState();
    TPZSkylineNSymStructMatrix matsp(cmesh);
	//TPZSpStructMatrix matsp(cmesh);
    
	std::set< int > materialid;
	int matid = mymaterial->MatId();
	materialid.insert(matid);
	matsp.SetMaterialIds (materialid);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
	TPZFMatrix<STATE> Un;
    
    TPZAutoPointer <TPZMatrix<STATE> > matK2 = matsp.CreateAssemble(Un,guiInterface);
    
    return matK2;
}

void StiffMatrixLoadVec(TPZMatConvectionProblem *mymaterial, TPZCompMesh*cmesh, TPZAnalysis &an, TPZAutoPointer< TPZMatrix<STATE> > &matK1, TPZFMatrix<STATE> &fvec){
    
	mymaterial->SetCurrentState();
    TPZSkylineStructMatrix matsk(cmesh);
    //TPZBandStructMatrix matsk(cmesh);
	an.SetStructuralMatrix(matsk);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt);
	an.SetSolver(step);
	an.Run();
	
	//matK1 = an.StructMatrix();
    matK1 = an.Solver().Matrix();
	fvec = an.Rhs();
}

void ForcingInicial(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    //double y = pt[1];
    if(x<0.0)disp[0] = 1.;
    if(x>0.0)disp[0] = 0.;
}

//Ativar apenas a ultima equacao, que corresponde a funcao de base constante
void FilterEquation(TPZMatConvectionProblem *mymaterial, TPZCompMesh *cmesh, TPZAnalysis &an, bool currentstate)
{
    int ncon_saturation = cmesh->NConnects();
    TPZManVector<long> active(0);
    for(int i = 0; i<ncon_saturation; i++)
    {
        TPZConnect &con = cmesh->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        if(seqnum==-1) continue;
        int pos = cmesh->Block().Position(seqnum);
        int blocksize = cmesh->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+1);
        int ieq = blocksize-1;
        active[vs] = pos+ieq;
        
    }
    if(currentstate==true)
    {
        mymaterial->SetCurrentState();
        TPZSkylineStructMatrix matsk(cmesh);
        matsk.SetNumThreads(4);
        matsk.EquationFilter().SetActiveEquations(active);
        an.SetStructuralMatrix(matsk);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();
    }
    else
    {
        mymaterial->SetLastState();
        TPZSkylineNSymStructMatrix matst(cmesh);
        //TPZSpStructMatrix matst(mphysics);
        //matsp.SetNumThreads(30);
        matst.EquationFilter().SetActiveEquations(active);
        an.SetStructuralMatrix(matst);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();
    }
}


