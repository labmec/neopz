/*
 *  StokesTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "HStokesTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZHStokesMaterial.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"

#ifdef USING_MKL
#include "TPZPardisoControl.h"
#include "mkl_pardiso.h"
#include "pzsysmp.h"
#include "pzysmp.h"
#endif

const REAL Pi=M_PI;

HStokesTest::HStokesTest()
{
    
}

HStokesTest::~HStokesTest()
{
    
}

void HStokesTest::Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE Sigma)
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
    
    TPZCompMesh *cmesh_v = this->CMesh_v(gmesh, Space, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = this->CMesh_p(gmesh, Space, pOrder); //Função para criar a malha computacional da pressão
    TPZCompMesh *cmesh_m = this->CMesh_m(gmesh, Space, pOrder, visco, theta, Sigma); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
        
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    AddMultiphysicsInterfaces(*cmesh_m);
    
#ifdef PZDEBUG
    {
        std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
        gmesh->Print(fileg1);
        
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
#endif
#ifdef PZDEBUG
    {
        std::ofstream fileg("MalhaGeoFinal.txt"); //Impressão da malha geométrica (formato txt)
        std::ofstream filegvtk("MalhaGeoFinal.vtk"); //Impressão da malha geométrica (formato vtk)
        gmesh->Print(fileg);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
    }
#endif

    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    
    //TPZFStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix matskl(cmesh_m); //OK para Hdiv
#else
    TPZSkylineNSymStructMatrix matskl(cmesh_m);
#endif

    matskl.SetNumThreads(numthreads);
    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatPoint);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);
    matids.insert(fmatIntBCbott);
    matids.insert(fmatIntBCright);
    matids.insert(fmatIntBCtop);
    matids.insert(fmatIntBCleft);
    matskl.SetMaterialIds(matids);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    
    
//#ifdef USING_MKL
//
//
//#endif
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    if(0)
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
//    std::cout << "Computing Error " << std::endl;
//    TPZManVector<REAL,3> Errors;
//    ofstream ErroOut("Erro.txt");
//    an.SetExact(Sol_exact);
//    an.PostProcessError(Errors,ErroOut);
    
    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors(6);
    ofstream ErroOut("Error_HStokes.txt", std::ofstream::app);
    an.SetExact(Sol_exact);
    
    cmesh_m->ElementSolution().Redim(cmesh_m->NElements(),6);


    
    an.PostProcessError(Errors);
    
    ErroOut <<"Sigma = "<< Sigma/(pOrder*pOrder*(nx-1)) << "  //  Ordem = "<< pOrder << "  //  Tamanho da malha = "<< nx-1 <<" x "<< ny-1 << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[1] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[2] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("Stokes.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    //        vecnames.Push("V_exactBC");
    
    
    int postProcessResolution = 3; //  keep low as possible
    
    int dim = gmesh->Dimension();
//    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
//    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}

TPZGeoMesh *HStokesTest::CreateGMesh(int nx, int ny, double hx, double hy)
{
    
    
    int i,j;
    int64_t id, index;
    
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    
    //Vetor auxiliar para armazenar coordenadas:
    
    TPZVec <REAL> coord (3,0.);
    
    
    //Inicialização dos nós:
    
    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = (i)*hy/(ny - 1);
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }
    
    //Ponto 1
    TPZVec<int64_t> pointtopology(1);
    pointtopology[0] = 0;
    gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);
    
    
    //Vetor auxiliar para armazenar as conecções entre elementos:
    
    TPZVec <int64_t> connect(4,0);
    
    
    //Conectividade dos elementos:
    
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            connect[0] = (i)*ny + (j);
            connect[1] = connect[0]+1;
            connect[2] = connect[1]+(nx);
            connect[3] = connect[0]+(nx);
            gmesh->CreateGeoElement(EQuadrilateral,connect,fmatID,id);
        }
    }
    
    
    //Gerando informação da vizinhança:

    gmesh->BuildConnectivity();
    {
        TPZCheckGeom check1(gmesh);
        check1.CheckUniqueId();
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
        TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
        TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
        
        for (int64_t i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[2] == 0. & nodecoord[1] == 0.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[1] == hy)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == 0.)
            {
                sizeOfleftVec++;
                ncoordzleftVec.Resize(sizeOfleftVec);
                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
            }
            if (nodecoord[2] == 0. & nodecoord[0] == hx)
            {
                sizeOfrightVec++;
                ncoordzrightVec.Resize(sizeOfrightVec);
                ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
            }
        }
        
        if (sizeOfbottVec == 2) {
            int sidesbott = plate->WhichSide(ncoordzbottVec);
            TPZGeoElSide platesidebott(plate, sidesbott);
            TPZGeoElBC(platesidebott,fmatBCbott);
            TPZGeoElBC(platesidebott,fmatIntBCbott);
        }
        
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,fmatBCtop);
            TPZGeoElBC(platesidetop,fmatIntBCtop);
        }
        
        if (sizeOfleftVec == 2) {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,fmatBCleft);
            TPZGeoElBC(platesideleft,fmatIntBCleft);
        }
        
        if (sizeOfrightVec == 2) {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,fmatBCright);
            TPZGeoElBC(platesideright,fmatIntBCright);
        }
        
        
        ncoordzbottVec.Resize(0);
        sizeOfbottVec = 0;
        ncoordztopVec.Resize(0);
        sizeOftopVec = 0;
        ncoordzleftVec.Resize(0);
        sizeOfleftVec = 0;
        ncoordzrightVec.Resize(0);
        sizeOfrightVec = 0;
        
    }
    
    // Criando e inserindo elemento de interfação:
    //    TPZVec<int64_t> nodind3(2);
    //
    //    nodind3[0]=1;
    //    nodind3[1]=4;
    //
    //    gmesh->CreateGeoElement(EOned, nodind3, matInterface, index); //Criando elemento de interface (GeoElement)
    
    
    //Criando interface (Geralizado):
    
    TPZVec<int64_t> nodint(2);
    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            if(j>0&&j<(nx-1)){
                nodint[0]=j+nx*i;
                nodint[1]=j+nx*(i+1);
                gmesh->CreateGeoElement(EOned, nodint, ftangentVelocity, index); //Criando elemento de interface (GeoElement)
                
            }
            if(i>0&&j<(ny-1)){
                nodint[0]=j+ny*i;
                nodint[1]=j+ny*i+1;
                gmesh->CreateGeoElement(EOned, nodint, ftangentVelocity, index); //Criando elemento de interface (GeoElement)
                
            }
            
        }
    }
    
    
    //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
    id++;
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    {
        ofstream bf("before.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    }
    return gmesh;
    

}

TPZCompEl *HStokesTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}



//TPZGeoMesh *StokesTest::GMeshDeformed(int dim, bool ftriang, int ndiv)
//{
//
//    DebugStop();
//
//}

void HStokesTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
    dsol.Resize(3,2);
    sol.Resize(3);
    
    REAL xv = x[0];
    REAL yv = x[1];
    
    STATE v_x =  cos(2*Pi*yv)*sin(2*Pi*xv);
    STATE v_y =  -(cos(2*Pi*xv)*sin(2*Pi*yv));
    STATE pressure= xv*xv+yv*yv;
    
    sol[0]=v_x;
    sol[1]=v_y;
    sol[2]=pressure;
    
    // vx direction
    dsol(0,0)= 2*Pi*cos(2*Pi*xv)*cos(2*Pi*yv);
    dsol(0,1)= 2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
    
    // vy direction
    dsol(1,0)= -2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
    dsol(1,1)= -2*Pi*cos(2*Pi*xv)*cos(2*Pi*yv);
    
    // Gradiente pressão
    
    dsol(2,0)= 2*xv;
    dsol(2,1)= 2*yv;
    
}



void HStokesTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(2);
    
    REAL xv = x[0];
    REAL yv = x[1];
    //    STATE zv = x[2];
    
    STATE f_x = 2.0*xv + 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
    STATE f_y = 2.0*yv - 8.0*Pi*Pi*cos(2.0*Pi*xv)*sin(2.0*Pi*yv);
    
    f[0] = f_x; // x direction
    f[1] = f_y; // y direction
    
    
}




TPZCompMesh *HStokesTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aproximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID); //Criando material que implementa a formulação fraca do problema modelo
    TPZMat2dLin *material2 = new TPZMat2dLin(fmatIDFlux);
    //Criando material que implementa a formulação fraca do problema modelo

    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(material2); //Insere material na malha

    if (Space==1) {
        cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        
        //Dimensões do material (para HDiv):
        TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
        material->SetMaterial(xkin, xcin, xfin);
        material2->SetMaterial(xkin, xcin, xfin);
        
    }else if(Space==2){
        
        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
        
        //Dimensões do material (para H1 e descontinuo):
        TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
        material->SetMaterial(xkin, xcin, xfin);
        
        
    }else if(Space==3){
        
        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
        //Criando elementos com graus de liberdade differentes para cada elemento (descontínuo):
        cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
        
        //Dimensões do material (para H1 e descontinuo):
        TPZFMatrix<STATE> xkin(2,2,0.), xcin(2,2,0.), xfin(2,2,0.);
        material->SetMaterial(xkin, xcin, xfin);
        
    }
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(2,1,0.);
    
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    }
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatIntBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatIntBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatIntBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatIntBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
        
        TPZMaterial * BCond4 = material->CreateBC(material, ftangentVelocity, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    std::set<int> matids;
    matids.insert(fmatID);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCright);
    matids.insert(fmatBCtop);
    matids.insert(fmatBCleft);
    cmesh->AutoBuild(matids);

//    int64_t nel = cmesh->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = cmesh->Element(el);
//        TPZCompElHDiv<pzshape::TPZShapeQuad> *celhdiv = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(cel);
    
//        if(0 && celhdiv)
//        {
//            TPZGeoEl *gel = celhdiv->Reference();
//            TPZGeoElBC gelbc(gel,gel->NSides()-1,fmatID);  //oioioi IDFlux -> ID
//            int64_t index;
//            TPZCompElHDiv<pzshape::TPZShapeQuad> *newel = new TPZCompElHDiv<pzshape::TPZShapeQuad>(*cmesh,*celhdiv,index);
//            newel->SetReference(gelbc.CreatedElement()->Index());
//        }
//    }
    
    gmesh->ResetReference();
    matids.clear();
    matids.insert(ftangentVelocity);
    matids.insert(fmatIntBCbott);
    matids.insert(fmatIntBCright);
    matids.insert(fmatIntBCtop);
    matids.insert(fmatIntBCleft);
    cmesh->AutoBuild(matids);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
   
    
}


TPZCompMesh *HStokesTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    if (Space==2||Space==3) {
        pOrder--;
    }
    //pOrder--;
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    
    // @omar::
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    //Criando material:
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    
    //Condições de contorno
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //    val2(0,0) = 0.0; // px -> 0
    //    val2(1,0) = 0.0; // py -> 0
    //
    //    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    //    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    //
    //    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    //    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    //
    //    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    //
    //    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    //    Ponto de pressao:
    //
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    ////
    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressao
    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    //    //    Ponto de pressao2:
    //    //
    //    TPZFMatrix<STATE> val5(1,1,0.), val6(1,1,0.);
    //    ////
    //    TPZMaterial * BCPoint2 = material->CreateBC(material, matPoint2, pointtype, val5, val6); //Cria material que implementa um ponto para a pressao
    //    cmesh->InsertMaterialObject(BCPoint2); //Insere material na malha
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    std::set<int> materialids;
    materialids.insert(fmatID);
    materialids.insert(fpointtype);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *HStokesTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE Sigma)
{
    
    //Criando malha computacional:
    //int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    
    // Criando material que implementa a formulacao fraca do problema modelo
    TPZStokesMaterial *material = new TPZHStokesMaterial(fmatID,fdim,Space,visco,theta,Sigma);
    TPZStokesMaterial *material2 = new TPZHStokesMaterial(fmatIDFlux,fdim,Space,visco,theta,Sigma);
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, 5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact, 5);
    ((TPZDummyFunction<STATE>*)fp.operator->())->SetPolynomialOrder(1);
    ((TPZDummyFunction<STATE>*)solp.operator->())->SetPolynomialOrder(1);
    
    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(solp);

    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(material2);

    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
        //BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
        //BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
        //BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
        //BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    }
    {
        TPZMaterial * BCond0 = material->CreateBC(material, fmatIntBCbott, fneumann, val1, val2); //Cria material que implementa a condição de contorno inferior
                                                                                                   //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
                                                                                                   //BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatIntBCtop, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
                                                                                                  //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
                                                                                                  //BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatIntBCleft, fneumann, val1, val2); //Cria material que implementa a condicao de contorno esquerda
                                                                                                   //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
                                                                                                   //BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatIntBCright, fneumann, val1, val2); //Cria material que implementa a condicao de contorno direita
                                                                                                    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
                                                                                                    //BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
        TPZMaterial * BCond4 = material->CreateBC(material, ftangentVelocity, fneumann, val1, val2); //Cria material que implementa a condicao de contorno direita
                                                                                                       //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
                                                                                                       //BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    }

    
    //Ponto
    
    TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    val4(0,0)=0.0;
    
    TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
    
    

    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}



void HStokesTest::AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
{
    // create interface elements between the tangent velocity element and the volumetric elements
    {
        TPZGeoMesh *gmesh = cmesh.Reference();
        std::set<int> velmatid;
        velmatid.insert(ftangentVelocity);
        velmatid.insert(fmatPoint);
        velmatid.insert(fmatIntBCtop);
        velmatid.insert(fmatIntBCbott);
        velmatid.insert(fmatIntBCleft);
        velmatid.insert(fmatIntBCright);
        int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            int matid = gel->MaterialId();
            if(velmatid.find(matid) != velmatid.end())
            {
                int nsides = gel->NSides();
                TPZGeoElSide gelside(gel,nsides-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    
                    TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
                    LeftElIndices[0]=0;
                    RightElIndices[0]=0;
                    
                    if (neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() == fmatID && gelside.Element()->Dimension() == 1) { //oioioi IDFlux -> ID
                        // create an interface element
                        TPZCompElSide celside = gelside.Reference();
                        TPZCompElSide celneigh = neighbour.Reference();
                        if (!celside || !celneigh) {
                            DebugStop();
                        }
                        std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                        " side " << neighbour.Side() << 
                        " and interface element " << gelside.Element()->Index() << std::endl;
                        TPZGeoElBC gelbc(gelside,fmatID);
                        int64_t index;
                        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,celneigh,celside);
                        intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                        
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
        }
    }
    
}





