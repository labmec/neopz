//
//  main_tese.cpp
//  PZ
//
//  Created by Agnaldo Farias on 2/12/13.
//
//

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

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"
#include "TPZVTKGeoMesh.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "mixedpoisson.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int MatId = 1;
const int bcdirichlet = 0;
const int bcneumann = 1;

const int bc0 = -1;
const int bc1 = -2;
const int bc2 = -3;
const int bc3 = -4;
bool fTriang = true;
const int eps=1000;
const REAL Pi=M_PI;
TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly,bool triang_elements);

TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh);
TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh);
void UniformRefine2(TPZGeoMesh* gmesh, int nDiv);

void PrintGMeshVTK2(TPZGeoMesh * gmesh, std::ofstream &file);

void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SaidaSolucao(TPZAnalysis &an, std::string plotfile);

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out);

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);

void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC3(const TPZVec<REAL> &loc, TPZVec<STATE> &result);
void NeumannBC4(const TPZVec<REAL> &loc, TPZVec<STATE> &result);

///Metodos para o ArcTangente
void RefineGeoElements2(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined);
void RegularizeMesh2(TPZGeoMesh *gmesh);
void GetPointsOnCircunference2(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);
void RefiningNearCircunference2(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
void SolArcTan2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux);
void ForcingTang2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);


int main3(int argc, char *argv[])
{
    //InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    TPZVec<REAL> erros;
    ofstream saidaerros("Erros.txt");
    saidaerros << "\nCalculo do Erro\n";
	
	
	
	
    for(int p=1; p<4; p++){
		
		//TPZGeoMesh * gmesh = MalhaGeom(1,1,fTriang);
		//UniformRefine2(gmesh,1);
		//std::ofstream filemesh("MalhaInicial.vtk");
		//TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);

		
		
        saidaerros << "Ordem "<<p<<" \n";
        
        int pq = p;
        int pp=p;

        
     //   saidaerros<<"\n CALCULO DO ERRO, COM ORDEM POLINOMIAL pq = " << pq << " e pp = "<< pp <<endl;
        for(int h=1; h<6; h++){
         //   saidaerros << "\nRefinamento: Ndiv = "<< h <<"\n";

            
            /*------------ Etapa 1 ------------*/
            
            // Criando a malha geométrica
            //TPZGeoMesh * gmesh = MalhaGeom(1,1,fTriang);
//			std::ofstream filemesh("MalhaInicial.vtk");
//			TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
//            UniformRefine2(gmesh,h);
			TPZGeoMesh * gmesh = MalhaGeom(1,1,fTriang);
			std::ofstream filemesh("MalhaGeoHPInicial.vtk");
			TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh, true);
			
			
			RefiningNearCircunference2(2,gmesh,h,1);
			
			

//			std::ofstream filemesh2("MalhaGeoHPRef.vtk");
//			 TPZVTKGeoMesh::PrintGMeshVTK(gmesh,filemesh2, true);
            
//            // Criando a primeira malha computacional
            TPZCompMesh * cmesh1= CMeshFlux(pq,gmesh);
//            
//            // Criando a segunda malha computacional
            TPZCompMesh * cmesh2 = CMeshPressure(pp,gmesh);

			
            // Criando a malha computacional multifísica
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(2);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            TPZCompMesh * mphysics = MalhaCompMultifisica(meshvec,gmesh);
			int nDofs;
            nDofs=mphysics->NEquations();
            
            saidaerros<< "\nNRefinamento "<<h<< "   NDofs "<<nDofs<<std::endl;
			
			
            /*---------- Etapa 4 ----------*/
            
            // Resolvendo o sistema linear
            TPZAnalysis an(mphysics);
            ResolverSistema(an, mphysics);
            
            /*---------- Etapa 5 ----------*/
			
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            // Arquivo de saida para plotar a solução
			string plotfile("Solution_ArcTanHP.vtk");
			SaidaSolucao(an, plotfile);
            
            saidaerros<<"\n\nErro da simulacao multifisica  para o flux";
            //TPZAnalysis an1(cmesh1);
            //an1.SetExact(*SolSuave);
            //an1.PostProcessError(erros, saidaerros);
            ErrorHDiv2(cmesh1,saidaerros);
            
            saidaerros<<"\n\nErro da simulacao multifisica  para a pressao";
			//            TPZAnalysis an2(cmesh2);
			//            an2.SetExact(*SolSuave);
			//            an2.PostProcessError(erros, saidaerros);
            ErrorH1(cmesh2, saidaerros);
			
        }
    }
    
	return EXIT_SUCCESS;
}

#include "pzgengrid.h"
TPZGeoMesh *MalhaGeom(REAL Lx, REAL Ly, bool triang_elements)
{
    TPZManVector<int,2> nx(2,1);
    TPZManVector<REAL,3> x0(3,0.),x1(3,0);
    x1[0] = Lx;
    x1[1] = Ly;
    TPZGenGrid gengrid(nx,x0,x1);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    if(fTriang)
    {
        gengrid.SetElementType(ETriangle);
    }
    gengrid.Read(gmesh);
    
    //elementos de contorno
    TPZManVector<REAL,3> firstpoint(3,0.),secondpoint(3,0.);
    
    secondpoint[0] = Lx;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc0);
    
    firstpoint = secondpoint;
    secondpoint[1] = Ly;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc1);
    
    firstpoint = secondpoint;
    secondpoint[0] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc2);
    
    firstpoint = secondpoint;
    secondpoint[1] = 0.;
    gengrid.SetBC(gmesh,firstpoint,secondpoint,bc3);
    
    return gmesh;
}


TPZCompMesh *CMeshFlux(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
    int dim = 2;
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(MatId,dim);
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    cmesh->InsertMaterialObject(mat);
    
    ///Criar condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_1 Fluxo\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    return cmesh;
}

TPZCompMesh *CMeshPressure(int pOrder,TPZGeoMesh *gmesh)
{
    /// criar materiais
    int dim = 2;
    TPZMatPoisson3d *material;
    material = new TPZMatPoisson3d(MatId,dim);
    material->NStateVariables();
	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
	
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
    TPZMaterial * BCond0 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
	
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
	
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(fTriang==true) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
        
    }
    
    return cmesh;
}


TPZCompMesh *MalhaCompMultifisica(TPZVec<TPZCompMesh *> meshvec,TPZGeoMesh * gmesh){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim =2;
    TPZMixedPoisson *material = new TPZMixedPoisson(MatId,dim);
    
    //incluindo os dados do problema
    REAL coefk = 1.;
    material->SetPermeability(coefk);
    REAL coefvisc = 1.;
    material->SetViscosity(coefvisc);
    
    //permeabilidade
    TPZFMatrix<REAL> Ktensor(dim,dim,0.);
    TPZFMatrix<REAL> InvK(dim,dim,0.);
    Ktensor(0,0)=1.; Ktensor(1,1)=1.;
    InvK=Ktensor;
    material->SetPermeabilityTensor(Ktensor,InvK);
    
    //solucao exata
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolArcTan2);
    material->SetForcingFunctionExact(solexata);
	
    //funcao do lado direito da equacao do problema
    TPZAutoPointer<TPZFunction<STATE> > force;
    TPZDummyFunction<STATE> *dum;
    dum = new TPZDummyFunction<STATE>(ForcingTang2);
    dum->SetPolynomialOrder(20);
    force = dum;
    material->SetForcingFunction(force);
    
    //inserindo o material na malha computacional
    TPZMaterial *mat(material);
    mphysics->InsertMaterialObject(mat);
    
    //Criando condicoes de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCond1 = material->CreateBC(mat, bc0,bcdirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc1,bcdirichlet, val1, val2);
    TPZMaterial * BCond3 = material->CreateBC(mat, bc2,bcdirichlet, val1, val2);
    TPZMaterial * BCond4 = material->CreateBC(mat, bc3,bcdirichlet, val1, val2);
    
    //Set force function
	//    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAbaixo;
	//    bcmatNeumannAbaixo = new TPZDummyFunction<STATE>(NeumannBC1);
	//    BCond1->SetForcingFunction(bcmatNeumannAbaixo);
	//    
	//    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannDirichlet;
	//    bcmatNeumannDirichlet = new TPZDummyFunction<STATE>(NeumannBC2);
	//    BCond2->SetForcingFunction(bcmatNeumannDirichlet);
	//
	//    TPZAutoPointer<TPZFunction<STATE> > bcmatNeumannAcima;
	//    bcmatNeumannAcima = new TPZDummyFunction<STATE>(NeumannBC3);
	//    BCond3->SetForcingFunction(bcmatNeumannAcima);
    
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    return mphysics;
}


#define VTK
void ResolverSistema(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<STATE> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
    
	//Saida de Dados: solucao e  grafico no VTK
	ofstream file("Solution.out");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	
}

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void SolSuave(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    u.Resize(1, 0.);
    du.Resize(3, 1.);
    du(0,0)=du(1,0)=du(2,0)=0.;
    
    const REAL sol = sin(M_PI*x)*sin(M_PI*y);
    u[0] = sol;
    
    du(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y);//-k*Grad(u)[[0]]
    du(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x);////-k*Grad(u)[[1]]
    du(2,0) = 2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}

void NeumannBC1(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,-1.};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannBC2(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {1.,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannBC3(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {0,1.};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void NeumannBC4(const TPZVec<REAL> &loc, TPZVec<STATE> &result){
    REAL normal[2] = {-1.,0};
    
    TPZManVector<REAL> u(1);
    TPZFNMatrix<10> du(2,1);
    SolSuave(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}


void SaidaSolucao(TPZAnalysis &an, std::string plotfile){
    
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
    
	scalnames[0] = "Pressure";
	scalnames[1] = "ExactPressure";
	vecnames[0]= "Flux";
	vecnames[1]= "ExactFlux";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

void PrintGMeshVTK2(TPZGeoMesh * gmesh, std::ofstream &file)
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

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out)
{
    long nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<STATE,10> elerror(10,0.);
        cel->EvaluateError(SolArcTan2, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
    
}

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out)
{
    long nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<STATE,10> globerrors(10,0.);
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<STATE,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(SolArcTan2, elerror, NULL);
        
        int nerr = elerror.size();
        //globerrors.resize(nerr);
#ifdef LOG4CXX
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "\nH1 Norm = "    << sqrt(globerrors[0]) << endl;
    out << "\nL2 Norm = "    << sqrt(globerrors[1]) << endl;
    out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
    out << "\n=============================\n"<<endl;
}


void UniformRefine2(TPZGeoMesh* gmesh, int nDiv)
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
    //	gmesh->BuildConnectivity();
}
///Metodos para o problema ArcTangente
void RefineGeoElements2(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}
void GetPointsOnCircunference2(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points) {
	Points.Resize(npoints);
	TPZManVector<REAL> point(3,0.);
	REAL angle = (2*M_PI)/npoints;
	for(int i=0;i<npoints;i++) {
		point[0] = center[0]+radius*cos(i*angle);
		point[1] = center[1]+radius*sin(i*angle);
		Points[i] = point;
	}
}
void RefiningNearCircunference2(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs) {
    
	int i;
	bool isdefined = false;
	
	// Refinando no local desejado
	int npoints = 1000;
	TPZVec<REAL> point(3);
	point[0] = point[1] = 0.5; point[2] = 0.0;
	REAL r = 0.25;
	TPZVec<TPZManVector<REAL> > Points(npoints);
	GetPointsOnCircunference2(npoints,point,r,Points);
	
	if(ntyperefs==2) {
		REAL radius = 0.19;
		for(i=0;i<nref;i+=2) {
			// To refine elements with center near to points than radius
			RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
			RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
			if(i < 5) radius *= 0.35;
			else if(i < 7) radius *= 0.2;
			else radius *= 0.1;
		}
		if(i==nref) {
			RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
		}
	}
	else {
		REAL radius = 0.2;
		for(i=0;i<nref+1;i++) {
			// To refine elements with center near to points than radius
			RefineGeoElements2(dim,gmesh,point,r,radius,isdefined);
			if(i < 8) radius *= 0.6;
			else if(i < 7) radius *= 0.3;
			else radius *= 0.15;
		}
	}
	// Constructing connectivities
	//gmesh->ResetConnectivities();
    RegularizeMesh2(gmesh);
	gmesh->BuildConnectivity();
}

void RegularizeMesh2(TPZGeoMesh *gmesh)
{
    bool changed = true;
    while (changed)
    {
        changed = false;
        int nel = gmesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (gel->HasSubElement()) {
                continue;
            }
            int dim = gel->Dimension();
            if (dim != 2) {
                continue;
            }
            int nsides = gel->NSides();
            int nrefined = 0;
            int nsidedim = 0;
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                nsidedim++;
            }
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                TPZGeoElSide thisside(gel,is);
                TPZGeoElSide neighbour = thisside.Neighbour();
                if (neighbour != thisside) {
                    TPZStack<TPZGeoElSide> subelements;
                    neighbour.GetSubElements2(subelements);
                    int nsub = subelements.size();
                    if (nsub > 0) {
                        nrefined++;
                    }
                    for (int isub=0; isub<nsub; isub++) {
                        TPZGeoElSide sub = subelements[isub];
                        if (sub.Dimension() != dim-1) {
                            continue;
                        }
                        if (sub.HasSubElement()) {
                            TPZManVector<TPZGeoEl *> newsub;
                            gel->Divide(newsub);
                            changed = true;
                            break;
                        }
                    }
                }
                if (gel->HasSubElement()) {
                    break;
                }
            }
            if (nrefined >= nsidedim-1) {
                TPZManVector<TPZGeoEl *> newsub;
                gel->Divide(newsub);
                changed = true;
            }
        }
    }
}

#define Power pow
#define ArcTan atan
#define Sqrt sqrt

void SolArcTan2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux){
    REAL x = pt[0];
    REAL y = pt[1];
    p[0]=0;
    flux(0,0)=0;
    flux(1,0)=0;
    //flux(2,0)=0;
    
    p[0]= 5.*(-1. + x)*x*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
	
	
	//px
    flux(0,0)=(12.732395447351628*Sqrt(eps)*(-1. + x)*(-0.5 + x)*x*(-1. + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) - 
	5.*(-1. + x)*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) - 
	5.*x*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
	
	
    //py    
    flux(1,0)= (12.732395447351628*Sqrt(eps)*(-1. + x)*x*(-1. + y)*(-0.5 + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) - 
	5.*(-1. + x)*x*(-1. + y)*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) - 
	5.*(-1. + x)*x*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))));
	
	
    
    
}


void ForcingTang2(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    //void ForcingTang(const TPZVec<REAL> &pt, TPZVec<REAL> &res,TPZFMatrix<STATE> &disp){
    //   disp.Redim(2,1);
    double x = pt[0];
    double y = pt[1];
    
	disp[0]= (-1.)*((3.183098861837907*Sqrt(eps)*(-1. + x)*x*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) - 
													   64.*eps*Power(-0.5 + x,2)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))*(-1. + y)*y)/
    Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) - 
	(5.092958178940651*Sqrt(eps)*(-0.5 + x)*(-5. + 10.*x)*(-1. + y)*y)/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) + 
	10.*(-1. + y)*y*(1. + 0.6366197723675814*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))) + 
	5.*(-1. + x)*x*(2. + (0.6366197723675814*Sqrt(eps)*(-4.*(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) - 
														64.*eps*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2))*Power(-0.5 + y,2))*(-1. + y)*y)/
					Power(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2),2) - 
					(5.092958178940651*Sqrt(eps)*(-0.5 + y)*(-1. + 2*y))/(1 + 4.*eps*Power(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2),2)) + 
					1.2732395447351628*ArcTan(2.*Sqrt(eps)*(0.0625 - 1.*Power(-0.5 + x,2) - 1.*Power(-0.5 + y,2)))));


}

