//
//  CircularCMesh.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "ComputationalMesh.hpp"

// Cria malha Computacional para malha 360 graus

TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder, int projection, int inclinedwellbore,
                           int analytic, REAL SigmaV, REAL Sigmah, REAL SigmaH, REAL Pwb, REAL rw,
                           REAL rext, REAL direction, REAL inclination, bool isStochastic,
                           int nSquareElements, TPZFMatrix<STATE> &M, REAL scale, int funcE, int funcnu, int distribE, int distribnu) {
    
    //criando material que implementa a formulacao fraca do problema modelo
    TPZMatElasticity2D *material = new TPZMatElasticity2D(MATERIAL_ID);
    
    // http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
    REAL Eyoung = 29269.00, ni = 0.203, fbx = 0., fby = 0.;
    
    // REAL Eyoung = 29269, ni = 0.203, fbx = 0., fby = 0.; //Dados tese Diogo
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    // transforma direcao e inclinacao de graus em rad
    REAL directionT = direction * (M_PI / 180);
    REAL inclinationT = inclination * (M_PI / 180);
    
    // Seta os parametros do poco (Inclined or not)
    material->SetInclinedWellboreParameters(SigmaH, Sigmah, SigmaV, directionT, inclinationT,
                                            inclinedwellbore, Pwb, rw, analytic, projection);
    
    //Dados tese Diogo
    REAL A = 152.54; //MPa
    REAL B = 0.0015489; // MPae-1
    REAL C = 146.29; //MPa
    material->SetSandlerDiMaggioParameters(A, B, C);
    
    //Dados tese Diogo, angulo de friccao e coesao
    REAL c = 27.38; //MPa
    REAL frict = 11.7; // graus
    REAL frictRad = frict * (M_PI / 180); // rad
    material->SetMogiAndMohrCoulombParameters(c, frictRad);
    
    //Obtem tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
    
    ///criar malha computacional
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(DIMENSION_2D);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    // cond contorno - ids igual da malha geometrica
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4;
    //bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6;
    
    // tipo de condicao de contorno
    const int normalpressure = 6, stressfield = 4, mixed = 2, dirichlet = 0;
    //neumann = 1;
    
    TPZFMatrix<REAL> val1(3,3,0.), val2(2,1,0.);
    
    if (analytic == 0) {
        // Inserir condicao de contorno parede do poco
        val1(0,0) = Pwb;
        val1(1,1) = Pwb;
        val1(2,2) = Pwb;
        TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material
        
        // Inserir condicao de contorno circunferencia externa
        val1(0,0) = SigmaX;
        val1(1,0) = SigmaXY;
        val1(0,1) = SigmaXY;
        val1(1,1) = SigmaY;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
        
        // Inserir condicao de contorno ponto externo bottom
        val1(0,0) = 0.0;
        val1(1,0) = 0.0;
        val1(0,1) = 0.0;
        val1(1,1) = 1.0;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond2 = material->CreateBC(material, bc2, mixed, val1, val2);//cria material
        
        // Inserir condicao de contorno ponto externo lateral direita
        val1(0,0) = 1.0;
        val1(1,0) = 0.0;
        val1(0,1) = 0.0;
        val1(1,1) = 0.0;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond3 = material->CreateBC(material, bc3, mixed, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
        
        cmesh->InsertMaterialObject(BCond0);//insere material na malha
        cmesh->InsertMaterialObject(BCond1);//insere material na malha
        cmesh->InsertMaterialObject(BCond2);//insere material na malha
        cmesh->InsertMaterialObject(BCond3);//insere material na malha
        
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
        cmesh->AutoBuild();
        
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
    else if (analytic == 1 || analytic == 2) {
        
        ///Inserir condicao de contorno circunferencia interna
        val1(0,0) = 0.; //SigmaX;
        val1(1,0) = 0.; //SigmaXY;
        val1(0,1) = 0.; //SigmaXY;
        val1(1,1) = 0.; //SigmaY;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond0 = material->CreateBC(material, bc0, stressfield, val1, val2);
        
        ///Inserir condicao de contorno circunferencia externa
        val1(0,0) = 0.; //SigmaX;
        val1(1,0) = 0.; //SigmaXY;
        val1(0,1) = 0.; //SigmaXY;
        val1(1,1) = 0.; //SigmaY;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);
        
        ///Inserir condicao de contorno ponto externo bottom
        val1(0,0) = 0.0;
        val1(1,0) = 0.0;
        val1(0,1) = 0.0;
        val1(1,1) = 1.0;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond2 = material->CreateBC(material, bc2, mixed, val1, val2);
        
        ///Inserir condicao de contorno ponto externo lateral direita
        val1(0,0) = 1.0;
        val1(1,0) = 0.0;
        val1(0,1) = 0.0;
        val1(1,1) = 0.0;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond3 = material->CreateBC(material, bc3, mixed, val1, val2);
        
        cmesh->InsertMaterialObject(BCond0);//insere material na malha
        cmesh->InsertMaterialObject(BCond1);
        cmesh->InsertMaterialObject(BCond2);
        cmesh->InsertMaterialObject(BCond3);
        
        cmesh->SetAllCreateFunctionsContinuous();
        
        //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
        cmesh->AutoBuild();
        cmesh->AdjustBoundaryElements();
        cmesh->CleanUpUnconnectedNodes();
    }
    else {
        DebugStop();
    }
    
    // Set Forcing Function for Stochastic Analysis
    if(isStochastic == true) {
        TPZAutoPointer<TPZFunction<STATE> > force = new TPZRandomField<STATE>(gmesh, nSquareElements, inclinedwellbore, direction, inclination, rw, rext, M, scale,funcE, funcnu,
                                                                              distribE,distribnu);
        material->SetForcingFunction(force);
        material->GetNSquareElements(nSquareElements);
    }
    
    return cmesh;
}

// *********** Cria malha Computacional para 1/4 do Poco **********************/
TPZCompMesh *QuarterCMesh(TPZGeoMesh *gmesh, int pOrder) {
    //criando material que implementa a formulacao fraca do problema modelo
    TPZMatElasticity2D *material = new TPZMatElasticity2D(MATERIAL_ID);
    
    // Setting up paremeters
    REAL Eyoung = 15300, ni = 0.24, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    REAL Pi = M_PI;
    REAL direction = 0., inclination = 0.; //graus
    direction = 0.;
    inclination = 0.;
    REAL directionT = 0.,inclinationT = 0.;
    directionT = direction*(Pi/180); // rad
    inclinationT = inclination*(Pi/180); // rad
    int inclinedwellbore = 1;
    
    REAL Pwb = -30.0; //MPa
    
    // Tensoes in Situ, horizontais e vertical em Pa
    REAL SigmaVV = 0., Sigmahh = 0., SigmaHH = 0.; // inicializa
    SigmaVV = -50.0, Sigmahh = -40.0, SigmaHH = -60.0; //preenche
    
    REAL rw = 1.0;
    int analytic = 1;
    int projection = 0;
    
    // Seta os parametros do poco
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT,
                                            inclinedwellbore, Pwb, rw, analytic, projection);
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(DIMENSION_2D);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    // cond contorno - ids igual da malha geometrica
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6, bc6 = -7, bc7 = -8;
    
    // tipo de condicao de contorno
    const int stressfield = 4, mixed = 2, neumann = 1, normalpressure = 6;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    ///Inserir condicao de contorno parede do poco
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond0 = material->CreateBC(material, bc0, stressfield, val1, val2);//cria material
    
    //    ///Inserir condicao de contorno parede do poco
    //    val1(0,0) = Pwb;
    //    val1(1,0) = 0.;
    //    val1(0,1) = 0.;
    //    val1(1,1) = Pwb;
    //    val2(0,0) = 0.;
    //    val2(1,0) = 0.;
    //    TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);
    
    ///Inserir condicao de contorno externo bottom
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);
    
    ///Inserir condicao de contorno externo upper
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, stressfield, val1, val2);
    
    ///Inserir condicao de contorno arco externo
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, stressfield, val1, val2);
    
    ///Inserir condicao de contorno ponto parede do poco bottom
    val1(0,0) = 1.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 0.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material, bc4, mixed, val1, val2);
    
    ///Inserir condicao de contorno ponto arco externo do poco bottom
    val1(0,0) = 1.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 0.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material, bc5, mixed, val1, val2);
    
    ///Inserir condicao de contorno ponto parede do poco top
    val1(0,0) = 0.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 1.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond6 = material->CreateBC(material, bc6, mixed, val1, val2);
    
    ///Inserir condicao de contorno ponto arco externo do poco top
    val1(0,0) = 0.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 1.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond7 = material->CreateBC(material, bc7, mixed, val1, val2);
    
    cmesh->InsertMaterialObject(BCond0);//insere material na malha
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
