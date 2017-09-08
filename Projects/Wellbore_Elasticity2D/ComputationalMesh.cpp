//
//  CircularCMesh.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "ComputationalMesh.hpp"

//************************************ Cria malha Computacional para malha 360 graus ***********************************************************//

TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder, TPZFMatrix<REAL> SetKCorr)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    // Setting up paremeters
    //  copy this link http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
    REAL Eyoung = 15300, ni = 0.24, fbx = 0., fby = 0.;
    //      REAL Eyoung = 29269, ni = 0.203, fbx = 0., fby = 0.; //Dados tese Diogo
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    // Pwb = pressao da lama em MPa
    REAL Pi = M_PI;
    
    /************ Define Posicao do Poco **************/
    REAL direction = 0., inclination = 0.; //inicializa angulos
    direction   = 0.; // graus******** 30
    inclination = 0.; // graus******** 50
    
    // transforma graus em rad
    REAL directionT = 0.,inclinationT = 0.; // inicializa
    directionT = direction*(Pi/180); // rad
    inclinationT = inclination*(Pi/180); // rad
    
    // define disposicao do poco
    // inclined == 1
    int inclinedwellbore = 0;
    
    // pressao da lama de perfuracao
    REAL Pwb = -10.5; // MPa
    //  REAL Pwb = -30.; // MPa
    
    
    // Tensoes in Situ, horizontais e vertical em MPa
    REAL SigmaVV = 0., Sigmahh = 0., SigmaHH = 0.; // inicializa
    SigmaVV = -50.0, Sigmahh = -40.0, SigmaHH = -60.0; //preenche
    //    SigmaVV = -48.2, Sigmahh = -62.1, SigmaHH = -45.9; //preenche //Dados Diogo //NANANANA INVERTIDAS AS TENSOES
    //     SigmaVV = -48.2, Sigmahh = -45.9, SigmaHH = -62.1; //preenche //Dados Diogo //NANANANA
    //  SigmaVV = -30.0, Sigmahh = -30.0, SigmaHH = -30.0; //preenche
    
    // REAL rw = 0.10795;
    REAL rw = 0.1;
    
    //analytic=0 nao usa sol analitica como prestress e BC
    //analytic=1 usa sol analitica como prestress e BC (zerar BCond0 e BCond1)
    //analytic=2 nao usa sol analitica como prestress mas usa como BC (zerar BCond0 e BCond1)
    int analytic = 0;
    
    // para projecao horizontal, projection == 1
    int projection = 0;
    
    // Seta os parametros do poco (Inclined or not)
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT, inclinedwellbore, Pwb, rw, analytic, projection);
    
    REAL A = 0., B = 0., C = 0.;
    
    //Dados tese Diogo
    A = 152.54; //MPa
    B = 0.0015489; // MPae-1
    C = 146.29; //MPa
    material->SetSandlerDiMaggioParameters(A, B, C);
    
    REAL c = 0., frict = 0., frictRad = 0.;
    
    //Dados tese Diogo, angulo de friccao e coesao
    c = 27.38; //MPa
    frict = 11.7; // graus
    
    frictRad = frict*(Pi/180); // rad
    material->SetMogiAndMohrCoulombParameters(c, frictRad);
    
    
    //Obtem tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
    //
    //#ifdef PZDEBUG
    //    #ifdef LOG4CXX
    //        if(logger->isDebugEnabled())
    //        {
    //
    //            std::stringstream out;
    //            out << " Stress rotation " << std::endl;
    //            out << " SigmaX     = " << SigmaX << std::endl;
    //            out << " SigmaXY    = " << SigmaXY << std::endl;
    //            out << " SigmaY     = " << SigmaY << std::endl;
    //            out << " SigmaZ     = " << SigmaZ << std::endl;
    //            LOGPZ_DEBUG(logger,out.str())
    //        }
    //    #endif
    //#endif
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    // cond contorno
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4; // ids igual da malha geometrica
    //bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6;
    const int normalpressure = 6, stressfield = 4, mixed = 2, dirichlet = 0; // tipo de condicao de contorno
    //neumann = 1;
    
    TPZFMatrix<REAL> val1(3,3,0.), val2(2,1,0.);
    
    if (analytic==0) {
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
    else if (analytic == 1 || 2) {
        
        ///Inserir condicao de contorno circunferencia interna
        val1(0,0) = 0.; //SigmaX;
        val1(1,0) = 0.; //SigmaXY;
        val1(0,1) = 0.; //SigmaXY;
        val1(1,1) = 0.; //SigmaY;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond0 = material->CreateBC(material, bc0, stressfield, val1, val2);//cria material
        
        ///Inserir condicao de contorno circunferencia externa
        val1(0,0) = 0.; //SigmaX;
        val1(1,0) = 0.; //SigmaXY;
        val1(0,1) = 0.; //SigmaXY;
        val1(1,1) = 0.; //SigmaY;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
        
        ///Inserir condicao de contorno ponto externo bottom
        val1(0,0) = 0.0;
        val1(1,0) = 0.0;
        val1(0,1) = 0.0;
        val1(1,1) = 1.0;
        val2(0,0) = 0.0;
        val2(1,0) = 0.0;
        TPZMaterial * BCond2 = material->CreateBC(material, bc2, mixed, val1, val2);//cria material
        
        ///Inserir condicao de contorno ponto externo lateral direita
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
    else {
        DebugStop();
    }
    
    
    
    //#ifdef LOG4CXX
    //    if (logger->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logger, sout.str())
    //    }
    //#endif
    
    
    // Set Forcing Function for Stochastic Analysis
    TPZAutoPointer<TPZFunction<STATE> > force  = new TPZRandomField<STATE>();
    material->SetForcingFunction(force);
    
    
    //St Correlation Matrix
    material->SetCorrelationMatrix(SetKCorr);
    
    return cmesh;
}

// *********** Cria malha Computacional para 1/4 do Poco **********************/
TPZCompMesh *QuarterCMesh(TPZGeoMesh *gmesh, int pOrder, TPZFMatrix<REAL> SetKCorr)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema
    
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
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
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT, inclinedwellbore, Pwb, rw, analytic, projection);
    
    
    
    //    //Obtem tensor de tensoes iniciais
    //    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    //    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
    //
    //
    //#ifdef PZDEBUG
    //#ifdef LOG4CXX
    //    if(logger->isDebugEnabled())
    //    {
    //
    //        std::stringstream out;
    //        out << " Stress rotation " << std::endl;
    //        out << " SigmaX     = " << SigmaX << std::endl;
    //        out << " SigmaXY    = " << SigmaXY << std::endl;
    //        out << " SigmaY     = " << SigmaY << std::endl;
    //        out << " SigmaZ     = " << SigmaZ << std::endl;
    //        LOGPZ_DEBUG(logger,out.str())
    //    }
    //#endif
    //#endif
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    
    //    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    //    REAL theta = 0.;
    //    theta = 0;
    //    material->AnalyticalWellboreSolution(SigmaX, SigmaY, SigmaXY, SigmaZ, theta, rw);
    
    
    // cond contorno
    const int bc0 = -1, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6, bc6 = -7, bc7 = -8; // ids igual da malha geometrica
    const int stressfield = 4, mixed = 2, neumann = 1, normalpressure = 6; // tipo de condicao de contorno
    
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
    //    TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material
    
    ///Inserir condicao de contorno externo bottom
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno externo upper
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno arco externo
    val1(0,0) = 0.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 0.;
    val2(0,0) = 0.;
    val2(1,0) = 0.;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto parede do poco bottom
    val1(0,0) = 1.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 0.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond4 = material->CreateBC(material, bc4, mixed, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto arco externo do poco bottom
    val1(0,0) = 1.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 0.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond5 = material->CreateBC(material, bc5, mixed, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    ///Inserir condicao de contorno ponto parede do poco top
    val1(0,0) = 0.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 1.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond6 = material->CreateBC(material, bc6, mixed, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto arco externo do poco top
    val1(0,0) = 0.0;
    val1(1,0) = 0.0;
    val1(0,1) = 0.0;
    val1(1,1) = 1.0;
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMaterial * BCond7 = material->CreateBC(material, bc7, mixed, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    
    
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
