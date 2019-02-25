/*
 *  BrinkmanTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "BrinkmanTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZGmshReader.h"

#define TRIANGLEMESH

using namespace std;

const REAL Pi=M_PI;

const REAL phi_r = 0.;

BrinkmanTest::BrinkmanTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=-1;
    fmatBCtop=-2;
    fmatBCleft=-3;
    fmatBCright=-4;
    
    //Material do elemento de interface
    fmatInterface=4;
    
    //Materiais das condições de contorno (elementos de interface)
    fmatIntBCbott=-11;
    fmatIntBCtop=-12;
    fmatIntBCleft=-13;
    fmatIntBCright=-14;
    
    //Materia de um ponto
    fmatPoint=-5;
    
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
    
    fSpaceV=0;
    
    fphi_r=0;
    
    f_is_hdivFull = false;
}

BrinkmanTest::~BrinkmanTest()
{
    
}

void BrinkmanTest::Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE sigma)
{
    
    
    //Gerando malha geométrica:
    fSpaceV = Space;
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
    
    //ChangeExternalOrderConnects(cmesh_v);
    
    TPZCompMesh *cmesh_m = this->CMesh_m(gmesh, Space, pOrder, visco, theta, sigma); //Função para criar a malha computacional multifísica
    
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

    if(fSpaceV!=2){
        AddMultiphysicsInterfaces(*cmesh_m,fmatInterface,fmatID);
        AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCbott,fmatBCbott);
        AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCtop,fmatBCtop);
        AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCleft,fmatBCleft);
        AddMultiphysicsInterfaces(*cmesh_m,fmatIntBCright,fmatBCright);
    }
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    
//            TPZSpStructMatrix struct_mat(cmesh_m);
//            struct_mat.SetNumThreads(numthreads);
//            an.SetStructuralMatrix(struct_mat);
    
    
    //TPZParSkylineStructMatrix matskl(cmesh_m, numthreads);

    TPZSkylineNSymStructMatrix matskl(cmesh_m); //OK para Hdiv
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    
//    if (Space==3) {
//        TPZFStructMatrix matsklD(cmesh_m); //caso nao simetrico *** //OK para discont.
//        matsklD.SetNumThreads(numthreads);
//        an.SetStructuralMatrix(matsklD);
//    }


    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    
    
    
    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;
    
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    
    

    
    std::cout << "Solving Matrix " << std::endl;
    
    an.Solve();
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
    }
#endif
    
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
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_Brinkman.txt", std::ofstream::app);
    an.SetExact(Sol_exact);
    an.PostProcessError(Errors,false);
    
    ErroOut <<"Sigma = "<< sigma/(pOrder*pOrder*(nx-1)) << "  //  Ordem = "<< pOrder << "  //  Tamanho da malha = "<< nx-1 <<" x "<< ny-1 << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[1] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[2] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[4] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("Brinkman.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    scalnames.Push("Div");


    int postProcessResolution = 3; //  keep low as possible

    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
}


void BrinkmanTest::Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate){
    
    if (rotate==true) {
        //rotação +
        co_r[0] = co[0]*cos(phi_r) - co[1]*sin(phi_r);
        co_r[1] = co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }else{
        
        co_r[0] = co[0]*cos(phi_r) + co[1]*sin(phi_r);
        co_r[1] = - co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }
    
    
}

TPZGeoMesh *BrinkmanTest::CreateGMesh(int nx, int ny, double hx, double hy)
{
    
#ifdef TRIANGLEMESH

    //Criando malha geométrica, nós e elementos.
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = dirname + "/Projects/CreepStokes/gmsh_meshes/GeometryTri.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    Geometry.fPZMaterialId[1]["bottom"] = fmatBCbott;
    Geometry.fPZMaterialId[1]["right"] = fmatBCright;
    Geometry.fPZMaterialId[1]["top"] = fmatBCtop;
    Geometry.fPZMaterialId[1]["left"] = fmatBCleft;
    
    gmesh = Geometry.GeometricGmshMesh(grid);
    
    
    // Criando e inserindo elemento de interface:
    
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        
        if (gel->Dimension() == 1) {
            
            
            int nsides = gel->NSides();
            
            TPZGeoElSide gelside(gel,nsides-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                
                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                    
                    break;
                    
                }
                neighbour = neighbour.Neighbour();
                
            }
            
            int mat_id = gel->MaterialId();
            int mat_bc_inter_id = 0;
            
            if (mat_id==-1) mat_bc_inter_id = fmatIntBCbott;
            if (mat_id==-2) mat_bc_inter_id = fmatIntBCtop;
            if (mat_id==-3) mat_bc_inter_id = fmatIntBCleft;
            if (mat_id==-4) mat_bc_inter_id = fmatIntBCright;
            if (mat_bc_inter_id==0) DebugStop();
            
            // tagging interfaces materials on boundaries
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, mat_bc_inter_id);
            }
            
            
        }
        
        if (gel->Dimension() != gmesh->Dimension()) {
            
            continue;
            
        }
        
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
                continue;
            }
            
            TPZGeoElSide gelside(gel,is);
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                
                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                    
                    break;
                    
                }
                neighbour = neighbour.Neighbour();
                
            }
            
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, fmatInterface);
            }
        }
    }
    
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    
    int n_div = 0;
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    

    
    
    
//    int i,j;
//    int64_t id, index;
//
//
//    //Criando malha geométrica, nós e elementos.
//    //Inserindo nós e elementos no objeto malha:
//
//    TPZGeoMesh *gmesh = new TPZGeoMesh();
//    gmesh->SetDimension(2);
//
//    //Vetor auxiliar para armazenar coordenadas:
//
//    TPZVec<REAL> coord (3,0.), coord_r(3,0.);
//
//    //Inicialização dos nós:
//
//    for(i = 0; i < ny; i++){
//        for(j = 0; j < nx; j++){
//            id = i*nx + j;
//            coord[0] = (j)*hx/(nx - 1);
//            coord[1] = -1.+(i)*hy/(ny - 1);
//            //using the same coordinate x for z
//            coord[2] = 0.;
//            //cout << coord << endl;
//
//            //rottação phi
//            Rotate(coord, coord_r, true);
//
//            //Get the index in the mesh nodes vector for the new node
//            index = gmesh->NodeVec().AllocateNewElement();
//            //Set the value of the node in the mesh nodes vector
//            gmesh->NodeVec()[index] = TPZGeoNode(id,coord_r,*gmesh);
//        }
//    }
//
//
//    //Ponto 1
//    TPZVec<int64_t> pointtopology(1);
//    pointtopology[0] = nx-1;
//
//    gmesh->CreateGeoElement(EPoint,pointtopology,fmatPoint,id);
//
//
//    //Vetor auxiliar para armazenar as conecções entre elementos:
//
//    TPZVec <int64_t> connectD(3,0);
//    TPZVec <int64_t> connectU(3,0);
//
//
//    //Conectividade dos elementos:
//
//    for(i = 0; i < (ny - 1); i++){
//        for(j = 0; j < (nx - 1); j++){
//            index = (i)*(nx - 1)+ (j);
//            connectD[0] = (i)*ny + (j);
//            connectD[1] = connectD[0]+1;
//            connectD[2] = connectD[0]+nx+1;
//            gmesh->CreateGeoElement(ETriangle,connectD,fmatID,id);
//
//            connectU[0] = connectD[2];
//            connectU[1] = connectD[2]-1;
//            connectU[2] = connectD[0];
//            gmesh->CreateGeoElement(ETriangle,connectU,fmatID,id);
//
//            //   std::cout<<connectD<<std::endl;
//            //   std::cout<<connectU<<std::endl;
//
//            id++;
//        }
//    }
//
//
//    //Gerando informação da vizinhança:
//
//    gmesh->BuildConnectivity();
//
//    {
//        TPZCheckGeom check(gmesh);
//        check.CheckUniqueId();
//    }
//    int64_t el, numelements = gmesh->NElements();
//
//    TPZManVector <int64_t> TopolPlate(4);
//
//    for (el=0; el<numelements; el++)
//    {
//        int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
//        TPZGeoEl *plate = gmesh->ElementVec()[el];
//        for (int i=0; i<4; i++){
//            TopolPlate[i] = plate->NodeIndex(i);
//        }
//
//        //Colocando as condicoes de contorno:
//        TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
//        TPZManVector <REAL,3> nodecoord(3,0.),nodecoord_r(3,0.);
//
//        //Na face x = 1
//        TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
//        TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
//        TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
//        TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
//
//
//        for (int64_t i = 0; i < totalnodes; i++)
//        {
//            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
//            Nodefinder[i].GetCoordinates(nodecoord_r);
//            int id_node = Nodefinder[i].Id();
//
//            for (int64_t j = 0; j < ny; j++){
//
//
//                if (id_node==j)
//                {
//                    sizeOfbottVec++;
//                    ncoordzbottVec.Resize(sizeOfbottVec);
//                    ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
//                }
//                if (id_node==j+nx*(nx-1))
//                {
//                    sizeOftopVec++;
//                    ncoordztopVec.Resize(sizeOftopVec);
//                    ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
//                }
//
//
//                if (id_node==j*nx)
//                {
//                    sizeOfleftVec++;
//                    ncoordzleftVec.Resize(sizeOfleftVec);
//                    ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
//                }
//                if (id_node==(j+1)*nx-1)
//                {
//                    sizeOfrightVec++;
//                    ncoordzrightVec.Resize(sizeOfrightVec);
//                    ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
//                }
//
//            }
//        }
//
//
//
////        for (int64_t i = 0; i < totalnodes; i++)
////        {
////            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
////            Nodefinder[i].GetCoordinates(nodecoord_r);
////
////            // Desrotacionar:
////
////            Rotate(nodecoord_r, nodecoord, false);
////
////
////            if (nodecoord[2] == 0. & nodecoord[1] == -1.+0.)
////            {
////                sizeOfbottVec++;
////                ncoordzbottVec.Resize(sizeOfbottVec);
////                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
////            }
////            if (nodecoord[2] == 0. & nodecoord[1] == -1.+hy)
////            {
////                sizeOftopVec++;
////                ncoordztopVec.Resize(sizeOftopVec);
////                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
////            }
////            if (nodecoord[2] == 0. & nodecoord[0] == 0.)
////            {
////                sizeOfleftVec++;
////                ncoordzleftVec.Resize(sizeOfleftVec);
////                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
////            }
////            if (nodecoord[2] == 0. & nodecoord[0] == hx)
////            {
////                sizeOfrightVec++;
////                ncoordzrightVec.Resize(sizeOfrightVec);
////                ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
////            }
////        }
//
//        if (sizeOfbottVec == 2) {
//            int sidesbott = plate->WhichSide(ncoordzbottVec);
//            TPZGeoElSide platesidebott(plate, sidesbott);
//            TPZGeoElBC(platesidebott,fmatBCbott);
//            if(fSpaceV!=2){
//                TPZGeoElBC(platesidebott,fmatIntBCbott);
//            }
//        }
//
//        if (sizeOftopVec == 2) {
//            int sidestop = plate->WhichSide(ncoordztopVec);
//            TPZGeoElSide platesidetop(plate, sidestop);
//            TPZGeoElBC(platesidetop,fmatBCtop);
//            if(fSpaceV!=2){
//                TPZGeoElBC(platesidetop,fmatIntBCtop);
//            }
//        }
//
//        if (sizeOfleftVec == 2) {
//            int sidesleft = plate->WhichSide(ncoordzleftVec);
//            TPZGeoElSide platesideleft(plate, sidesleft);
//            TPZGeoElBC(platesideleft,fmatBCleft);
//            if(fSpaceV!=2){
//                TPZGeoElBC(platesideleft,fmatIntBCleft);
//            }
//        }
//
//        if (sizeOfrightVec == 2) {
//            int sidesright = plate->WhichSide(ncoordzrightVec);
//            TPZGeoElSide platesideright(plate, sidesright);
//            TPZGeoElBC(platesideright,fmatBCright);
//            if(fSpaceV!=2){
//                TPZGeoElBC(platesideright,fmatIntBCright);
//            }
//        }
//
//
//        ncoordzbottVec.Resize(0);
//        sizeOfbottVec = 0;
//        ncoordztopVec.Resize(0);
//        sizeOftopVec = 0;
//        ncoordzleftVec.Resize(0);
//        sizeOfleftVec = 0;
//        ncoordzrightVec.Resize(0);
//        sizeOfrightVec = 0;
//
//    }
//
//
//    //Criando interface (Geralizado):
//
//    if(fSpaceV!=2){
//        TPZVec<int64_t> nodint(2);
//        for(i = 0; i < (ny - 1); i++){
//            for(j = 0; j <= (nx - 1); j++){
//                if(j>0&&j<(nx-1)){
//                    nodint[0]=j+nx*i;
//                    nodint[1]=j+nx*(i+1);
//                    gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
//
//                }
//
//
//                if(i>0&&j<(ny-1)){
//                    nodint[0]=j+ny*i;
//                    nodint[1]=j+ny*i+1;
//                    gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
//
//                }
//
//                if(j<(nx-1)&&i<(ny-1)){
//                    nodint[0]=j+nx*i;
//                    nodint[1]=j+nx*i+nx+1;
//                    gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
//
//                }
//
//            }
//        }
//    }
//
//    //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
//    id++;
//
//    gmesh->AddInterfaceMaterial(fquadmat1, fquadmat2, fquadmat3);
//    gmesh->AddInterfaceMaterial(fquadmat2, fquadmat1, fquadmat3);
//
//    TPZCheckGeom check(gmesh);
//    check.CheckUniqueId();
//
//    gmesh->BuildConnectivity();
//
//    //Impressão da malha geométrica:
//
//    ofstream bf("before.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
//    return gmesh;
    
#else
    
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
            coord[1] = -1 + (i)*hy/(ny - 1);
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
    pointtopology[0] = nx-1;
    
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
        TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
        TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;
        
        for (int64_t i = 0; i < totalnodes; i++)
        {
            Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
            Nodefinder[i].GetCoordinates(nodecoord);
            if (nodecoord[1] == -1.)
            {
                sizeOfbottVec++;
                ncoordzbottVec.Resize(sizeOfbottVec);
                ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
            }
            if (nodecoord[1] == -1.+hy)
            {
                sizeOftopVec++;
                ncoordztopVec.Resize(sizeOftopVec);
                ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
            }
            if (nodecoord[0] == 0.)
            {
                sizeOfleftVec++;
                ncoordzleftVec.Resize(sizeOfleftVec);
                ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
            }
            if (nodecoord[0] == hx)
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
            if(fSpaceV!=2){
                TPZGeoElBC(platesidebott,fmatIntBCbott);
            }
        }
        
        if (sizeOftopVec == 2) {
            int sidestop = plate->WhichSide(ncoordztopVec);
            TPZGeoElSide platesidetop(plate, sidestop);
            TPZGeoElBC(platesidetop,fmatBCtop);
            if(fSpaceV!=2){
                TPZGeoElBC(platesidetop,fmatIntBCtop);
            }
        }
        
        if (sizeOfleftVec == 2) {
            int sidesleft = plate->WhichSide(ncoordzleftVec);
            TPZGeoElSide platesideleft(plate, sidesleft);
            TPZGeoElBC(platesideleft,fmatBCleft);
            if(fSpaceV!=2){
                TPZGeoElBC(platesideleft,fmatIntBCleft);
            }
        }
        
        if (sizeOfrightVec == 2) {
            int sidesright = plate->WhichSide(ncoordzrightVec);
            TPZGeoElSide platesideright(plate, sidesright);
            TPZGeoElBC(platesideright,fmatBCright);
            if(fSpaceV!=2){
                TPZGeoElBC(platesideright,fmatIntBCright);
            }
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
    if(fSpaceV!=2){
        TPZVec<int64_t> nodint(2);
        for(i = 0; i < (ny - 1); i++){
            for(j = 0; j < (nx - 1); j++){
                if(j>0&&j<(nx-1)){
                    nodint[0]=j+nx*i;
                    nodint[1]=j+nx*(i+1);
                        gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
                    
                }
                if(i>0&&j<(ny-1)){
                    nodint[0]=j+ny*i;
                    nodint[1]=j+ny*i+1;
                        gmesh->CreateGeoElement(EOned, nodint, fmatInterface, index); //Criando elemento de interface (GeoElement)
                    
                }
                
            }
        }
    }
    
    //new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (nodind3,matInterface,*gmesh); //Criando elemento de interface (RefPattern)
    id++;
    
 //   gmesh->AddInterfaceMaterial(fquadmat1, fquadmat2, fquadmat3);
 //   gmesh->AddInterfaceMaterial(fquadmat2, fquadmat1, fquadmat3);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    //Impressão da malha geométrica:
    
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

#endif

}

TPZCompEl *BrinkmanTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void BrinkmanTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
//    dsol.Resize(3,2);
//    sol.Resize(3);
//
//    REAL x1 = x[0];
//    REAL x2 = x[1];
//
//    REAL e = exp(1.);
//
//    STATE v_1 = (1.-2./e)*sin(x1)*sin(x2);
//    STATE v_2 = -1.*cos(x1)*cos(x2);
//    STATE pressure= cos(x1)*sin(x2);
//
//    sol[0]=v_1;
//    sol[1]=v_2;
//    sol[2]=pressure;
//
//    // vx direction
//    dsol(0,0)= (1.-2./e)*cos(x1)*sin(x2);
//    dsol(0,1)= cos(x2)*sin(x1);
//
//    // vy direction
//    dsol(1,0)= (1.-2./e)*cos(x2)*sin(x1);
//    dsol(1,1)= cos(x1)*sin(x2);
//
//    // Gradiente pressão
//    dsol(2,0)= -sin(x1)*sin(x2);
//    dsol(2,1)= cos(x1)*cos(x2);

    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
    dsol.Resize(3,2);
    sol.Resize(3);

    REAL x1 = x[0];
    REAL x2 = x[1];

    REAL e = exp(1.);

    STATE v_1 = -1.*sin(x1)*sin(x2);
    STATE v_2 = -1.*cos(x1)*cos(x2);
    STATE pressure= cos(x1)*sin(x2);

    sol[0]=v_1;
    sol[1]=v_2;
    sol[2]=pressure;

    // vx direction
    dsol(0,0)= -1.*cos(x1)*sin(x2);
    dsol(0,1)= cos(x2)*sin(x1);

    // vy direction
    dsol(1,0)= -1.*cos(x2)*sin(x1);
    dsol(1,1)= cos(x1)*sin(x2);

    // Gradiente pressão
    dsol(2,0)= -sin(x1)*sin(x2);
    dsol(2,1)= cos(x1)*cos(x2);
    

    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
//        dsol.Resize(3,2);
//        sol.Resize(3);
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//
//        STATE v_1 = sin(x1)*sin(x2);
//        STATE v_2 = -1.*cos(x1)*cos(x2);
//        STATE pressure= cos(x1)*sin(x2);
//
//        sol[0]=v_1;
//        sol[1]=v_2;
//        sol[2]=pressure;
//
//        // vx direction
//        dsol(0,0)= cos(x1)*sin(x2);
//        dsol(0,1)= cos(x2)*sin(x1);
//
//        // vy direction
//        dsol(1,0)= cos(x2)*sin(x1);
//        dsol(1,1)= cos(x1)*sin(x2);
//
//        // Gradiente pressão
//        dsol(2,0)= -sin(x1)*sin(x2);
//        dsol(2,1)= cos(x1)*cos(x2);

    
}

void BrinkmanTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    f.resize(3);
    REAL x1 = x[0];
    REAL x2 = x[1];
    STATE f_1 =0., f_2=0.;
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
//    REAL e = exp(1.);
//
//    f_1 = (-8./e+ 4.)*sin(x1)*sin(x2);
//    f_2 = (2./e- 4.)*cos(x1)*cos(x2);
//    STATE g_1 = 2.*(1.-1./e)*cos(x1)*sin(x2);
//
//    f[0] = f_1; // x direction
//    f[1] = f_2; // y direction
//
//    f[2] = g_1; // g source

    // Stokes : : Artigo Botti, Di Pietro, Droniou
    

    f_1 = -3.*sin(x1)*sin(x2);
    f_2 = -1.*cos(x1)*cos(x2);

    f[0] = f_1; // x direction
    f[1] = f_2; // y direction
    f[2] = 0.;
    
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
//        f_1 = 0.;
//        f_2 = 0.;
//
//        f[0] = f_1; // x direction
//        f[1] = f_2; // y direction
//        f[2] = 2.*cos(x1)*sin(x2);



}

void BrinkmanTest::ChangeExternalOrderConnects(TPZCompMesh *mesh){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    int cordermin = -1;
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        
        if(cel->Dimension()== dim){
            for (int icon=0; icon<ncon-1; icon++){
                TPZConnect &co  = cel->Connect(icon);
                corder = co.Order();
                nshape = co.NShape();
                if(corder!=cordermin){
                    cordermin = corder-1;
                    int64_t cindex = cel->ConnectIndex(icon);
                    co.SetOrder(cordermin,cindex);
                    co.SetNShape(nshape-1);
                    mesh->Block().Set(co.SequenceNumber(),nshape-1);
                }
            }
        }
    }
    mesh->ExpandSolution();
    mesh->CleanUpUnconnectedNodes();
}


TPZCompMesh *BrinkmanTest::CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder)
{
 
    //Criando malha computacional:
    //pOrder++;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(fmatID); //Criando material que implementa a formulação fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    if (Space==1) {
        cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
        
        if (f_is_hdivFull == true) {
           cmesh->ApproxSpace().CreateDisconnectedElements(true); //HDIV-Full:
        }

        
        //Dimensões do material (para HDiv):
        TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
        material->SetMaterial(xkin, xcin, xfin);
        
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
    
    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
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


TPZCompMesh *BrinkmanTest::CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder)
{
    
    if (Space==2||Space==3) {
        pOrder--;
    }

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

TPZCompMesh *BrinkmanTest::CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE sigma)
{
    
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(fdim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    
    TPZBrinkmanMaterial *material = new TPZBrinkmanMaterial(fmatID,fdim,Space,visco,theta,sigma);//criando material que implementa a formulacao fraca do problema modelo
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (F_source, 5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_exact, 5);
    ((TPZDummyFunction<STATE>*)fp.operator->())->SetPolynomialOrder(5);
    ((TPZDummyFunction<STATE>*)solp.operator->())->SetPolynomialOrder(5);
    
    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(solp);

    cmesh->InsertMaterialObject(material);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
  
        TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbott, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
        BCond0->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond0); //Insere material na malha
        
        TPZMaterial * BCond1 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
        BCond1->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond1); //Insere material na malha
        
        TPZMaterial * BCond2 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
        BCond2->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond2); //Insere material na malha
        
        TPZMaterial * BCond3 = material->CreateBC(material, fmatBCright, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
        BCond3->SetForcingFunction(Sol_exact,bc_inte_order);
        cmesh->InsertMaterialObject(BCond3); //Insere material na malha
        //Ponto
        
        TPZFMatrix<STATE> val3(1,1,0.), val4(1,1,0.);
    
        val4(0,0)=-cos(2.)*sin(1.);
        
        
        TPZMaterial * BCPoint = material->CreateBC(material, fmatPoint, fpointtype, val3, val4); //Cria material que implementa um ponto para a pressão
        cmesh->InsertMaterialObject(BCPoint); //Insere material na malha
   

  


    

    int ncel = cmesh->NElements();
    for(int i = 0; i<ncel; i++){
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


void BrinkmanTest::AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
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
            DebugStop();
        }
        gel->SetMaterialId(mattarget);
        int64_t index;
        new TPZMultiphysicsInterfaceElement(cmesh,gel,index,celstack[1],celstack[0]);
    }
    
}





