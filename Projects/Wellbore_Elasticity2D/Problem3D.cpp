//
//  Problem3D.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "Problem3D.hpp"

//*********************************************** PROBELMA 3D ******************************************************//
int Problem3D(){
    
    bool Is3DQ = false;
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/Wellbore_Elasticity2D/";
    FileName += "Wellbore_Elasticity2DLog.cfg";
    InitializePZLOG(FileName);
#endif
    
    std::string grid = dirname;
    grid = grid + "/Projects/Wellbore_Elasticity2D/";
    //grid += "SingleWellRef.dump";
    grid += "CirularHole.dump";
    TPZGeoMesh *gmesh = ReadGeoMesh(grid,2);
    
    const std::string nm("Single_Well");
    gmesh->SetName(nm);
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
    
    
    
    //******** Configura malha Computacional ***************/
    
    int p = 1;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CMesh3D(gmesh, p, Is3DQ); //funcao para criar a malha COMPUTACIONAL de todo o poco
    
    // Solving linear equations
    // Initial steps
    TPZAnalysis an (cmesh);
    int numthreads = 2;
    std::cout << "Entering into Assemble ..." << std::endl;
    std::cout << "number of dof = " << cmesh->NEquations() << std::endl;
    
    
    bool UseIterativeSolverQ = true;
    
    if (UseIterativeSolverQ) {
        TPZSkylineStructMatrix skylstr(cmesh); //caso simetrico
        skylstr.SetNumThreads(numthreads);
        an.SetStructuralMatrix(skylstr);
        
        TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = skylstr.Create();
        TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
        
        TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>(matClone);
        TPZStepSolver<STATE> *Solver = new TPZStepSolver<STATE>(matbeingcopied);
        precond->SetReferenceMatrix(matbeingcopied);
        precond->SetDirect(ECholesky);
        Solver->SetCG(10, *precond, 1.0e-10, 0);
        an.SetSolver(*Solver);
    }
    else {
        TPZSkylineStructMatrix strskyl(cmesh);
        strskyl.SetNumThreads(numthreads);
        an.SetStructuralMatrix(strskyl);
        TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
        direct->SetDirect(ECholesky);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
    }
    
    an.Assemble();
    
    
    std::cout << "Entering into Solver ..." << std::endl;
    an.Solve();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    
    std::cout << "Entering into Postprocess ..." << std::endl;
    //    TPZFMatrix<REAL> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
    //    solucao.Print("Sol",cout,EMathematicaInput);//imprime na formatacao do Mathematica
    
    
    // Post processing
    int ndiv = 1;
    int dimension = gmesh->Dimension();
    TPZStack<std::string> scalarnames, vecnames;
    std::string name;
    
    if(Is3DQ){
        scalarnames.Push("StressX");
        scalarnames.Push("StressY");
        scalarnames.Push("StressZ");
        vecnames.Push("Displacement");
        name = "ElasticitySolutions3D.vtk";
    }
    else{
        scalarnames.Push("SigmaX");
        scalarnames.Push("SigmaY");
        scalarnames.Push("SigmaZ");
        vecnames.Push("Displacement");
        name = "ElasticitySolutions2D.vtk";
    }
    
    an.DefineGraphMesh(dimension,scalarnames,vecnames,name);
    
    an.PostProcess(ndiv);
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
    
}

/******************************************************* MALHA COMPUTACIONAL 3D ****************************************************/

TPZCompMesh *CMesh3D(TPZGeoMesh *gmesh, int pOrder, bool Is3DQ){
    
    int matId = 1;
    int dim;
    if(Is3DQ){
        dim = 3; //dimensao do problema
    }
    else{
        dim = 2; //dimensao do problema
    }
    
    
    
    if(Is3DQ){
        
        //**************** Criando material  ********************************
        TPZElasticity3D *material = new TPZElasticity3D(matId);//criando material que implementa a formulacao fraca do problema modelo
        
        
        // Setting up paremeters
        //  copy this link http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
        REAL Eyoung = 15.3e+9 , ni = 0.24, fbx = 0., fby = 0., fbz = 0.0;//-2500*9.81;
        
        TPZManVector<STATE> f(3,0);
        f[0] = fbx;
        f[1] = fby;
        f[2] = fbz;
        material->SetMaterialDataHook(Eyoung, ni);
        material->SetForce(f);
        
        
        
        
        /******* Calculating Inicial Stresses *******/
        // direction = direction/azimuth
        // inclination = wellbore inclination
        // problem assumption, inclined wellbore state = 1
        // Pwb = pressao da lama em MPa
        REAL Pi = M_PI;
        REAL direction = 0., inclination = 0.; // graus
        REAL directionT   = direction*(Pi/180); // rad
        REAL inclinationT = inclination*(Pi/180); // rad
        int inclinedwellbore = 1;
        REAL Pwb = 30.0e+6; // Pa
        
        // Tensoes in Situ, horizontais e vertical em Pa
        REAL SigmaVV = -50.0e6, Sigmahh = -40.0e6, SigmaHH = -60.0e6;
        
        material->SetPreStress(SigmaHH, Sigmahh, SigmaVV);
        
        ///criar malha computacional
        TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
        cmesh->SetDimModel(dim);//seta dimensao do modelo
        
        // Inserindo material na malha
        cmesh->InsertMaterialObject(material);
        
        int bcw, bce, bcs, bcn, bcb, bct, bcwell;
        
        // Matrial ids for boundaries 3D case
        
        bcw = 2;
        bce = 3;
        bcs = 4;
        bcn = 5;
        bcb = 6;
        bct = 7;
        bcwell = 8;
        
        // Matrial ids for boundaries 3D case
        
        const int stressfield = 4, neumann = 1, fixed_u = 0; // tipo de condicao de contorno
        
        TPZFMatrix<REAL> val1(3,3,0.0), val2(3,1,0.0);
        
        ///Inserir condicao de contorno parede do poco
        val1(0,0) = Pwb;
        val1(1,1) = Pwb;
        val1(2,2) = Pwb;
        TPZMaterial * BCond1 = material->CreateBC(material, bcwell, stressfield, val1, val2);//cria material
        
        val1(0,0) = -1.0*SigmaHH;
        val1(1,1) = -1.0*Sigmahh;
        val1(2,2) = -1.0*SigmaVV;
        TPZMaterial * BCond2 = material->CreateBC(material, bct, stressfield, val1, val2);//cria material
        
        
        val2(0,0) = 0;
        val2(1,0) = 0;
        val2(2,0) = 0;
        TPZMaterial * BCond3 = material->CreateBC(material, bcb, fixed_u, val1, val2);//cria material
        
        val1(0,0) = -1.0*SigmaHH;
        val1(1,1) = -1.0*Sigmahh;
        val1(2,2) = -1.0*SigmaVV;
        TPZMaterial * BCond4 = material->CreateBC(material, bcw, stressfield, val1, val2);
        TPZMaterial * BCond5 = material->CreateBC(material, bce, stressfield, val1, val2);
        TPZMaterial * BCond6 = material->CreateBC(material, bcs, stressfield, val1, val2);
        TPZMaterial * BCond7 = material->CreateBC(material, bcn, stressfield, val1, val2);
        
        cmesh->InsertMaterialObject(BCond1);//insere material na malha
        cmesh->InsertMaterialObject(BCond2);//insere material na malha
        cmesh->InsertMaterialObject(BCond3);//insere material na malha
        cmesh->InsertMaterialObject(BCond4);//insere material na malha
        cmesh->InsertMaterialObject(BCond5);//insere material na malha
        cmesh->InsertMaterialObject(BCond6);//insere material na malha
        cmesh->InsertMaterialObject(BCond7);//insere material na malha
        
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
        cmesh->AutoBuild();
        
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            cmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        return cmesh;
        
    }
    else {
        
        //**************** Criando material  ********************************
        TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
        
        // Setting up paremeters
        //  copy this link http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
        REAL Eyoung = 15.3e+9 , ni = 0.24, fbx = 0., fby = 0., fbz = 0.0;//-2500*9.81;
        
        TPZManVector<STATE> f(3,0);
        f[0] = fbx;
        f[1] = fby;
        f[2] = fbz;
        
        material->SetElasticity(Eyoung, ni, fbx, fby);
        
        //        /******* Calculating Inicial Stresses *******/
        //        // direction = direction/azimuth
        //        // inclination = wellbore inclination
        //        // problem assumption, inclined wellbore state = 1
        //        // Pwb = pressao da lama em MPa
        //        REAL Pi = M_PI;
        //        REAL direction = 0., inclination = 0.; // graus
        //        REAL directionT   = direction*(Pi/180); // rad
        //        REAL inclinationT = inclination*(Pi/180); // rad
        //        int inclinedwellbore = 1;
        //        REAL Pwb = 30.0e+6; // Pa
        //
        // Tensoes in Situ, horizontais e vertical em Pa
        REAL SigmaVV = 50.0e+6, Sigmahh = -40.0e6, SigmaHH = -60.0e6;
        //
        //        material->SetPreStress(SigmaHH, Sigmahh, SigmaVV);
        
        ///criar malha computacional
        TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
        cmesh->SetDimModel(dim);//seta dimensao do modelo
        
        // Inserindo material na malha
        cmesh->InsertMaterialObject(material);
        
        int bcw, bce, bcs, bcn, bcwell;
        
        // Material ids for boundaries 3D case
        bcw = 2;
        bce = 3;
        bcs = 4;
        bcn = 5;
        bcwell = 6;
        
        // Matrial ids for boundaries 2D case
        const int pressure = 6, neumann = 1, fixed_u = 0; // tipo de condicao de contorno
        
        TPZFMatrix<REAL> val1(2,2,0.0), val2(2,1,0.0);
        
        ///Inserir condicao de contorno parede do poco
        val1(0,0) = 0.0;
        val1(1,1) = 0.0;
        TPZMaterial * BCond1 = material->CreateBC(material, bcwell, pressure, val1, val2);//cria material
        
        val1.Zero();
        TPZMaterial * BCond2 = material->CreateBC(material, bce, fixed_u, val1, val2);
        TPZMaterial * BCond3 = material->CreateBC(material, bcw, fixed_u, val1, val2);
        
        val2(0,0) = 0.0;
        val2(1,0) = +1.0*SigmaVV;
        TPZMaterial * BCond4 = material->CreateBC(material, bcn, neumann, val1, val2);
        
        val2(0,0) = 0.0;
        val2(1,0) = -1.0*SigmaVV;
        TPZMaterial * BCond5 = material->CreateBC(material, bcs, neumann, val1, val2);
        
        cmesh->InsertMaterialObject(BCond1);//insere material na malha
        cmesh->InsertMaterialObject(BCond2);//insere material na malha
        cmesh->InsertMaterialObject(BCond3);//insere material na malha
        cmesh->InsertMaterialObject(BCond4);//insere material na malha
        cmesh->InsertMaterialObject(BCond5);//insere material na malha
        
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Cria elementos computacionais que gerenciarao o espaco de aprox. da malha
        cmesh->AutoBuild();
        
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            cmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        return cmesh;
        
    }
}

/******************* Le malha do GID ******************************************/

TPZGeoMesh * ReadGeoMesh(std::string GridFileName, int dim)
{
    TPZReadGIDGrid GeometryInfo;
    REAL s = 1.0;
    GeometryInfo.SetfDimensionlessL(s);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    gmesh->SetDimension(dim);
    return gmesh;
}

