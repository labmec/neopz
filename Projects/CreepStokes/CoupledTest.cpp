/*
 *  CoupledTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "CoupledTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"

CoupledTest::CoupledTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    fmatIdS = 1;
    fmatIdD = 2;
    
    //Materiais das condições de contorno (Stokes)
    fmatSBCbott = -1;
    fmatSBCtop = -2;
    fmatSBCleft = -3;
    fmatSBCright = -4;
    
    //Materiais das condições de contorno (Darcy)
    fmatDBCbott = -21;
    fmatDBCtop = -22;
    fmatDBCleft = -23;
    fmatDBCright = -24;
    
    //Material do elemento de interface
    fmatSInterface = 4;
    fmatDInterface = 14;
    fmatInterfaceDS= 5;
    
    //Materiais das condições de contorno (elementos de interface - Stokes)
    fmatIntSBCbott=-11;
    fmatIntSBCtop=-12;
    fmatIntSBCleft=-13;
    fmatIntSBCright=-14;
    
    //Materiais das condições de contorno (elementos de interface - Darcy)
    fmatIntDBCbott=-31;
    fmatIntDBCtop=-32;
    fmatIntDBCleft=-33;
    fmatIntDBCright=-34;
    
    //Materia de um ponto
    fmatPoint=-5;
    fmatPoint2 =-6;
    
    //Condições de contorno do problema
    fdirichlet=0;
    fneumann=1;
    fpenetration=2;
    fpointtype=5;
    fdirichletvar=4;
    
    
    fquadmat1=1; //Parte inferior do quadrado
    fquadmat2=2; //Parte superior do quadrado
    fquadmat3=3; //Material de interface
    
    fviscosity=1.;
    fpermeability=1.;
    ftheta=-1.;
    
    
}

CoupledTest::~CoupledTest()
{
    
}

void CoupledTest::Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE permeability, STATE theta)
{
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(nx, ny, hx, hy); //Função para criar a malha geometrica
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, Space, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, Space, pOrder); //Função para criar a malha computacional da pressão
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, Space, pOrder, visco, permeability, theta); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    
    AddMultiphysicsInterfaces(*cmesh_m,fmatSInterface,fmatIdS);
    AddMultiphysicsInterfaces(*cmesh_m,fmatDInterface,fmatIdD);
    
    //AddInterfaceCoupllingDS(gmesh,cmesh_m,fmatInterfaceDS,fmatIdS,fmatIdD);
    AddInterfaceCoupllingDS(gmesh,cmesh_m,fmatInterfaceDS,fmatIdD,fmatIdS);
    
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntSBCbott,fmatSBCbott);
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntSBCtop,fmatSBCtop);
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntSBCleft,fmatSBCleft);
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntSBCright,fmatSBCright);
    
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntDBCbott,fmatDBCbott);
    //AddMultiphysicsInterfaces(*cmesh_m,fmatIntDBCtop,fmatDBCtop);
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntDBCleft,fmatDBCleft);
    AddMultiphysicsInterfaces(*cmesh_m,fmatIntDBCright,fmatDBCright);
    
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZFStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    //TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
    std::cout << "Solving Matrix " << std::endl;
    
    an.Solve();
    
    
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,3> Errors;
    ofstream ErroOut("Erro.txt");
    an.SetExact(Sol_exact);
    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, ErroOut);
    
    
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("CouplingDS.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    vecnames.Push("P_exact");
    //        vecnames.Push("V_exactBC");
    
    
    int postProcessResolution = 0; //  keep low as possible
    
    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}

TPZGeoMesh *CoupledTest::CreateGMesh(int nx, int ny, double hx, double hy)
{

    int i,j;
    int64_t id, idno;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    
    
    //Inicialização dos nós:
    
    int64_t nnodes = nx*ny; //numero de nos do problema
    gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
    
    
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            idno = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = -hy/2.+(i)*hy/(ny - 1);
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[idno].SetCoord(coord);
            //Get the index in the mesh nodes vector for the new node
            //Set the value of the node in the mesh nodes vector
            //gmesh->NodeVec()[idno] = TPZGeoNode(id,coord,*gmesh);
            gmesh->NodeVec()[idno].SetNodeId(idno);
        }
    }
    
    //    Ponto 1
    TPZVec<int64_t> pointtopology(1);
    //pointtopology[0] = nx*(ny-1)/2;
    pointtopology[0] = 0;
    
    gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);
    
    //    Ponto 2
    TPZVec<int64_t> pointtopology2(1);
    pointtopology2[0] = nx*ny-1;
    //pointtopology2[0] = nx*((ny+1)/2);
    
    gmesh->CreateGeoElement(EPoint,pointtopology2,fmatPoint2,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <int64_t> connect(4,0);
    
    int64_t idD, idS;
    
    //Conectividade dos elementos (Parte Darcy):
    
    for(i = 0; i < (ny - 1)/2.; i++){
        for(j = 0; j < (nx - 1); j++){
            idD = (i)*(nx - 1)+ (j);
            connect[0] = (i)*nx + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,fmatIdD,id);
        }
    }
    
    //Conectividade dos elementos (Parte Stokes):
    
    for(i = (ny - 1)/2.; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            idS = (i)*(nx - 1)+ (j);
            connect[0] = (i)*nx + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,fmatIdS,id);
            
        }
    }
    
    
    //Gerando informação da vizinhança:
    
    gmesh->BuildConnectivity();
    
    {
        TPZCheckGeom check(gmesh);
        check.CheckUniqueId();
    }
    int64_t el, numelements = gmesh->NElements();
    
    TPZManVector <int64_t> TopolPlate(4);
    
    for (el=0; el<numelements; el++)
    {
        int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
        TPZGeoEl *plate = gmesh->ElementVec()[el];
        for (int i=0; i<4; i++){
            TopolPlate[i] = plate->NodeIndex(i);
        }
        
        //Colocando as condicoes de contorno:
        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
        TPZManVector <REAL,3> nodecoord(3);
        
        //Na face x = 1
        TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
        TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
        TPZVec<int64_t> ncoordzleftupVec(0); int64_t sizeOfleftupVec = 0;
        TPZVec<int64_t> ncoordzleftdownVec(0); int64_t sizeOfleftdownVec = 0;
        TPZVec<int64_t> ncoordzrightupVec(0); int64_t sizeOfrightupVec = 0;
        TPZVec<int64_t> ncoordzrightdownVec(0); int64_t sizeOfrightdownVec = 0;
        
        for (int64_t i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == -hy/2.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            
            if (nodecoord[2] == 0. & nodecoord[1] == hy/2.)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.& nodecoord[1] <= 0)
            {
                sizeOfleftdownVec++;
                ncoordzleftdownVec.Resize(sizeOfleftdownVec);
                ncoordzleftdownVec[sizeOfleftdownVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.& nodecoord[1] >= 0)
            {
                sizeOfleftupVec++;
                ncoordzleftupVec.Resize(sizeOfleftupVec);
                ncoordzleftupVec[sizeOfleftupVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx & nodecoord[1] <= 0)
            {
                sizeOfrightdownVec++;
                ncoordzrightdownVec.Resize(sizeOfrightdownVec);
                ncoordzrightdownVec[sizeOfrightdownVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx & nodecoord[1] >= 0)
            {
                sizeOfrightupVec++;
                ncoordzrightupVec.Resize(sizeOfrightupVec);
                ncoordzrightupVec[sizeOfrightupVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2) {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,fmatDBCbott);
            TPZGeoElBC(platesidebott,fmatIntDBCbott);
        }
    
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,fmatSBCtop);
            TPZGeoElBC(platesidetop,fmatIntSBCtop);
        }
        
        if (sizeOfleftdownVec == 2) {
            int sidesleftdown = plate->WhichSide(ncoordzleftdownVec);
            TPZGeoElSide platesideleftdown(plate, sidesleftdown);
            TPZGeoElBC(platesideleftdown,fmatDBCleft);
            TPZGeoElBC(platesideleftdown,fmatIntDBCleft);
        }
        
        if (sizeOfleftupVec == 2) {
            int sidesleftup = plate->WhichSide(ncoordzleftupVec);
            TPZGeoElSide platesideleftup(plate, sidesleftup);
            TPZGeoElBC(platesideleftup,fmatSBCleft);
            TPZGeoElBC(platesideleftup,fmatIntSBCleft);
        }
        
        if (sizeOfrightdownVec == 2) {
            int sidesrightdown = plate->WhichSide(ncoordzrightdownVec);
            TPZGeoElSide platesiderightdown(plate, sidesrightdown);
            TPZGeoElBC(platesiderightdown,fmatDBCright);
            TPZGeoElBC(platesiderightdown,fmatIntDBCright);
        }
        
        if (sizeOfrightupVec == 2) {
            int sidesrightup = plate->WhichSide(ncoordzrightupVec);
            TPZGeoElSide platesiderightup(plate, sidesrightup);
            TPZGeoElBC(platesiderightup,fmatSBCright);
            TPZGeoElBC(platesiderightup,fmatIntSBCright);
        }
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;

        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftdownVec.Resize(0);
        sizeOfleftdownVec = 0;
        ncoordzleftupVec.Resize(0);
        sizeOfleftupVec = 0;
        ncoordzrightdownVec.Resize(0);
        sizeOfrightdownVec = 0;
        ncoordzrightupVec.Resize(0);
        sizeOfrightupVec = 0;
        
    }

    
    //Criando interface (Geralizado):
    
    TPZVec<int64_t> nodint(2);
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            
            if(j>0&&j<(nx-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*(i+1);
                if(i<(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, fmatDInterface, id); //Criando elemento de interface (GeoElement)
                }
                if(i>=(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, fmatSInterface, id); //Criando elemento de interface (GeoElement)
                }
                
            }
            if(i>0&&i<(ny-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*i+1;
                
                if(i==(ny-1)/2.){
                    gmesh->CreateGeoElement(EOned, nodint, fmatInterfaceDS, id); //Criando elemento de interface (GeoElement)
                }
                if(i<(ny-1)/2.){
                    
                    gmesh->CreateGeoElement(EOned, nodint, fmatDInterface, id); //Criando elemento de interface (GeoElement)
                    
                }
                if(i>(ny-1)/2.){
                    
                    gmesh->CreateGeoElement(EOned, nodint, fmatSInterface, id); //Criando elemento de interface (GeoElement)
                }
                
                
            }
        }
    }
    
    
    gmesh->AddInterfaceMaterial(fmatIdS, fmatIdD, fmatInterfaceDS);
    gmesh->AddInterfaceMaterial(fmatIdD, fmatIdS, fmatInterfaceDS);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    
    
}

TPZCompEl *CoupledTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void CoupledTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    const REAL Pi=M_PI;
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    //Stokes
    if(yv>0.){
        
        STATE v_x =  (2.*cos(xv)*cos(Pi*yv)*sin(Pi*yv))/Pi;
        STATE v_y =  sin(xv)*(-2. + pow(sin(Pi*yv),2)/(Pi*Pi));
        STATE p =  sin(xv)*sin(yv);
        
        sol[0]=v_x;
        sol[1]=v_y;
        sol[2]=p;
        
        // vx direction
        dsol(0,0)= -(2.*cos(Pi*yv)*sin(xv)*sin(Pi*yv))/Pi;
        dsol(0,1)= 2.*cos(xv)*cos(Pi*yv)*cos(Pi*yv) + 2.*cos(xv)*sin(Pi*yv)*sin(Pi*yv);
        
        // vy direction
        dsol(1,0)= cos(xv)*(-2. + sin(Pi*yv)*sin(Pi*yv)/(Pi*Pi));
        dsol(1,1)= (2.*cos(Pi*yv)*sin(xv)*sin(Pi*yv))/Pi;
        
        // Gradiente pressão
        
        dsol(2,0)= cos(xv)*sin(yv);
        dsol(2,1)= cos(yv)*sin(xv);
        
        
    }
    //Darcy
    if(yv<0.){
        
        STATE v_x = -(-exp(-yv) + exp(yv))*cos(xv);
        STATE v_y = -(exp(-yv) + exp(yv))*sin(xv);
        STATE p = (-exp(-yv) + exp(yv))*sin(xv);
        
        sol[0]=v_x;
        sol[1]=v_y;
        sol[2]=p;
        
        // vx direction
        dsol(0,0)= -(-exp(-yv) + exp(yv))*sin(xv);
        dsol(0,1)= (exp(-yv) + exp(yv))*cos(xv);
        
        // vy direction
        dsol(1,0)= (exp(-yv) + exp(yv))*cos(xv);
        dsol(1,1)= (-exp(-yv) + exp(yv))*sin(xv);
        
        // Gradiente pressão
        
        dsol(2,0)= v_x;
        dsol(2,1)= v_y;
        
        
    }
    
    
}

void CoupledTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(2);
    
    REAL xv = x[0];
    REAL yv = x[1];
    const REAL Pi=M_PI;
    
    if(yv>0.){
        
        STATE f_x = (cos(xv)*(Pi*sin(yv) + (1. + 4.*Pi*Pi)*sin(2.*Pi*yv)))/Pi;
        STATE f_y = -((-1. + 4.*Pi*Pi - 2*Pi*Pi*cos(yv) + (1. + 4.*Pi*Pi)*cos(2.*Pi*yv))*
                      sin(xv))/(2.*Pi*Pi);
        
        f[0] = f_x;
        f[1] = f_y;
    }else{
        
        f[0] = 0.;
        f[1] = 0.;
    }
 
}


TPZCompMesh *CoupledTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{

    if (Space==2||Space==3) {
        DebugStop();
    }

    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    //cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    
    
    //Criando elementos com graus de liberdade differentes para cada elemento (descontínuo):
    
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Darcy:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialDarcy = new TPZMat2dLin(fmatIdD); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialDarcy); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
    //materialDarcy->SetMaterial(xkin, xcin, xfin);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    materialDarcy->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * BCDond0 = materialDarcy->CreateBC(materialDarcy, fmatDBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCDond0); //Insere material na malha
    
    // TPZMaterial * BCDond1 = materialDarcy->CreateBC(materialDarcy, fmatDBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    // cmesh->InsertMaterialObject(BCDond1); //Insere material na malha
    
    TPZMaterial * BCDond2 = materialDarcy->CreateBC(materialDarcy, fmatDBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCDond2); //Insere material na malha
    
    TPZMaterial * BCDond3 = materialDarcy->CreateBC(materialDarcy, fmatDBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCDond3); //Insere material na malha
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialStokes = new TPZMat2dLin(fmatIdS); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialStokes); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin2(2,2,0.), xcin2(2,2,0.), xfin2(2,2,0.);
    //materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin2(1,1,0.), xcin2(1,1,0.), xfin2(1,1,0.);
    materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    
    //Condições de contorno:
    
    //TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
//    TPZMaterial * BCSond0 = materialStokes->CreateBC(materialStokes, fmatSBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
//    cmesh->InsertMaterialObject(BCSond0); //Insere material na malha
    
    TPZMaterial * BCSond1 = materialStokes->CreateBC(materialStokes, fmatSBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCSond1); //Insere material na malha
    
    TPZMaterial * BCSond2 = materialStokes->CreateBC(materialStokes, fmatSBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCSond2); //Insere material na malha
    
    TPZMaterial * BCSond3 = materialStokes->CreateBC(materialStokes, fmatSBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCSond3); //Insere material na malha
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Coupling:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat1dLin *materialCoupling = new TPZMat1dLin(fmatInterfaceDS); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialCoupling); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    //TPZFMatrix<STATE> xkin3(2,2,0.), xcin3(2,2,0.),xbin3(2,2,0.), xfin3(2,2,0.);
    //materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
    //Dimensões do material (para HDiv):
    TPZFMatrix<STATE> xkin3(1,1,0.), xcin3(1,1,0.), xbin3(1,1,0.), xfin3(1,1,0.);
    materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    
    return cmesh;
    
}


TPZCompMesh *CoupledTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    if (Space==2||Space==3) {
        DebugStop();
    }
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    // @omar::
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Darcy:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialDarcy = new TPZMat2dLin(fmatIdD);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialDarcy); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin1(1,1,0.), xcin1(1,1,0.), xfin1(1,1,0.);
    materialDarcy->SetMaterial(xkin1, xcin1, xfin1);
    
    //Condições de contorno:
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //    val2(0,0) = 0.0; // px -> 0
    //    val2(1,0) = 0.0; // py -> 0
    //
    //    TPZMaterial * BCond0 = material->CreateBC(material, matDBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //    TPZMaterial * BCond1 = material->CreateBC(material, matDBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //    TPZMaterial * BCond2 = material->CreateBC(material, matDBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //    TPZMaterial * BCond3 = material->CreateBC(material, matDBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    
    TPZMaterial * BCPointD = materialDarcy->CreateBC(materialDarcy, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    cmesh->InsertMaterialObject(BCPointD); //Insere material na malha
    //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *materialStokes = new TPZMat2dLin(fmatIdS);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(materialStokes); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin2(1,1,0.), xcin2(1,1,0.), xfin2(1,1,0.);
    materialStokes->SetMaterial(xkin2, xcin2, xfin2);
    
    //Condições de contorno:
    
    
    //TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //    val2(0,0) = 0.0; // px -> 0
    //    val2(1,0) = 0.0; // py -> 0
    //
    //    TPZMaterial * BCond0 = material->CreateBC(material, matSBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //    TPZMaterial * BCond1 = material->CreateBC(material, matSBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //    TPZMaterial * BCond2 = material->CreateBC(material, matSBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //    TPZMaterial * BCond3 = material->CreateBC(material, matSBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //Ponto de pressao:
    
    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    
    TPZMaterial * BCPointS = materialStokes->CreateBC(materialStokes, fmatPoint2, fpointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    cmesh->InsertMaterialObject(BCPointS); //Insere material na malha
    
    
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Coupling:
    
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat1dLin *materialCoupling = new TPZMat1dLin(fmatInterfaceDS); //Criando material que implementa a formulação fraca do problema modelo
    
    //cmesh->InsertMaterialObject(materialCoupling); //Insere material na malha
    
    //Dimensões do material (para H1 e descontinuo):
    TPZFMatrix<STATE> xkin3(1,1,0.), xcin3(1,1,0.), xbin3(1,1,0.), xfin3(1,1,0.);
    materialCoupling->SetMaterial(xkin3, xcin3, xbin3, xfin3);
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatIdD);
    materialids.insert(fmatIdS);
    materialids.insert(fmatInterfaceDS);
    //materialids.insert(matInterfaceDS);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    std::set<int> materialids2;
    materialids2.insert(fmatPoint);
    //materialids2.insert(fmatpoint2);
    cmesh->AutoBuild(materialids2);
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    //   cmesh->AdjustBoundaryElements();
    //   cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;

    
}

TPZCompMesh *CoupledTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE permeability, STATE theta)
{

    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    // Criando material Darcy:
    
    TPZDarcyPMaterial *materialDarcy = new TPZDarcyPMaterial(fmatIdD,fdim,Space,visco,permeability,theta);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, 5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact, 5);
    
    materialDarcy->SetForcingFunction(fp);
    materialDarcy->SetForcingFunctionExact(solp);
    cmesh->InsertMaterialObject(materialDarcy);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCDond0 = materialDarcy->CreateBC(materialDarcy, fmatDBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    BCDond0->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond0); //Insere material na malha
    
    //    TPZMaterial * BCDond1 = materialDarcy->CreateBC(materialDarcy, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //    BCDond1->SetForcingFunction(sol_exact,bc_inte_order);
    //    cmesh->InsertMaterialObject(BCDond1); //Insere material na malha
    
    TPZMaterial * BCDond2 = materialDarcy->CreateBC(materialDarcy, fmatDBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    BCDond2->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond2); //Insere material na malha
    
    TPZMaterial * BCDond3 = materialDarcy->CreateBC(materialDarcy, fmatDBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    BCDond3->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCDond3); //Insere material na malha
    
    //Ponto
    
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    val4(0,0)=0.0;
    //
    TPZMaterial * BCPointD = materialDarcy->CreateBC(materialDarcy, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    cmesh->InsertMaterialObject(BCPointD); //Insere material na malha
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Stokes:
    
    TPZStokesMaterial *materialStokes = new TPZStokesMaterial(fmatIdS,fdim,Space,visco,theta,0.);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp2 = new TPZDummyFunction<STATE> (F_source, 5);
    TPZAutoPointer<TPZFunction<STATE> > solp2 = new TPZDummyFunction<STATE> (Sol_exact, 5);
    
    materialStokes->SetForcingFunction(fp2);
    materialStokes->SetForcingFunctionExact(solp2);
    cmesh->InsertMaterialObject(materialStokes);
    
    
    //Condições de contorno:
    
    //TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    //    TPZMaterial * BCSond0 = materialStokes->CreateBC(materialStokes, fmatSBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //    BCSond0->SetForcingFunction(solucaoS_exact,bc_inte_order);
    //    cmesh->InsertMaterialObject(BCSond0); //Insere material na malha
    
    TPZMaterial * BCSond1 = materialStokes->CreateBC(materialStokes, fmatSBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    BCSond1->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond1); //Insere material na malha
    
    TPZMaterial * BCSond2 = materialStokes->CreateBC(materialStokes, fmatSBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    BCSond2->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond2); //Insere material na malha
    
    TPZMaterial * BCSond3 = materialStokes->CreateBC(materialStokes, fmatSBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    BCSond3->SetForcingFunction(Sol_exact,bc_inte_order);
    cmesh->InsertMaterialObject(BCSond3); //Insere material na malha
    
    //Ponto
    
    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    val6(0,0)=0.0;
    
    TPZMaterial * BCPointS = materialStokes->CreateBC(materialStokes, fmatPoint2, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    cmesh->InsertMaterialObject(BCPointS); //Insere material na malha
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Criando material Acoplamento:
    
    TPZCouplingDSMaterial *materialCoupling = new TPZCouplingDSMaterial(fmatInterfaceDS,fdim,visco,permeability,theta);
    cmesh->InsertMaterialObject(materialCoupling);
    
    TPZAutoPointer<TPZFunction<STATE> > solp3 = new TPZDummyFunction<STATE> (Sol_exact, 5);
    
    //materialCoupling->SetForcingFunction(fp3);
    materialCoupling->SetForcingFunctionExact(solp3);
    cmesh->InsertMaterialObject(materialCoupling);
    
    
    
#ifdef PZDEBUG
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
#endif
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}


void CoupledTest::AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{

    TPZGeoMesh *gmesh = cmesh.Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != matfrom) {
            continue;
        }
        
        int nsides= gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        
        if (celstack.size() != 2) {
#ifdef PZDEBUG
            std::ofstream filecm("MalhaCPARCIAL_m.txt"); //Impressão da malha computacional multifísica (formato txt)
            cmesh.Print(filecm);
#endif
            
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        int64_t index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}

void CoupledTest::AddInterfaceCoupllingDS(TPZGeoMesh *gmesh,TPZCompMesh *cmesh, int matInterfaceDS ,int matleft, int matright){
    
    const int nel = gmesh->NElements();
    int64_t index;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        const int gelMatId = gel->MaterialId();
        if (gelMatId != matInterfaceDS)
            continue;
        if (gel->HasSubElement())
            continue;
        TPZCompElSide left, right;
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        int debug = 0;
        while (gelside != neigh) {
            if (neigh.Element()->MaterialId() == matleft) {
                TPZCompEl *cel = neigh.Element()->Reference();
                if (!cel) {
                    DebugStop();
                }
                left = TPZCompElSide(cel, 6);
                debug++;
            }
            else if (neigh.Element()->MaterialId() == matright) {
                TPZCompEl *cel = neigh.Element()->Reference();
                if (!cel) {
                    DebugStop();
                }
                right = TPZCompElSide(cel, 4);
                debug++;
            }
            neigh = neigh.Neighbour();
        }
        if (debug != 2)
            DebugStop();
        new TPZMultiphysicsInterfaceElement(*cmesh, gel, index, left, right);
    }
    
}



