//
//  Wellbore_Elasticity_2D_360d.cpp
//  PZ
//
//  Created by Nathalia on 6/7/16.
//
//

#include "Wellbore_Elasticity_2D_360d.hpp"



int Problem2D(){
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    FileName = dirname + "/Projects/Wellbore_Elasticity2D/";
    FileName += "Wellbore_Elasticity2DLog.cfg";
    InitializePZLOG(FileName);
#endif
    
    
    
    //******** Configura malha geometrica ***************/
    // rw = raio do poco (metros)
    // rext = raio externo do contorno (metros)
    // ncircle = nro elementos em 1/4 da parede do poco
    // nradial = nro de elementos da parede do poco ate o raio externo
    // drdcirc = proporcao do primeiro elemento
    REAL rw = 0.1;
    REAL rext = 5.0;
    int ncircle = 20;
    int nradial = 20;
    REAL drdcirc = 1.5;
    
    
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc); //funcao para criar a malha GEOMETRICA de todo o poco
    
    
    const std::string nm("line");
    gmesh->SetName(nm);
    
#ifdef LOG4CXX
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
#endif
    
    
    //******** Configura malha Computacional ***************/
    
    int p = 1;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *cmesh = CircularCMesh(gmesh, p); //funcao para criar a malha COMPUTACIONAL de todo o poco
    
    // Solving linear equations
    // Initial steps
    TPZAnalysis an (cmesh);
    int numthreads = 1;
    
    bool UseIterativeSolverQ = false;
    
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
    else{
        TPZSkylineStructMatrix strskyl(cmesh);
        strskyl.SetNumThreads(numthreads);
        an.SetStructuralMatrix(strskyl);
        TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
        direct->SetDirect(ECholesky);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
    }
    
    
    std::cout << "Entering into Assemble ..." << std::endl;
    std::cout << "number of dof = " << cmesh->NEquations() << std::endl;
    an.Assemble();
    
    //  an.Rhs() ;
    
    //    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    //    TPZFMatrix<STATE> FGlobal;
    //    KGlobal =   an.Solver().Matrix();
    //    FGlobal =   an.Rhs();
    //
    //#ifdef PZDEBUG
    ////    #ifdef LOG4CXX
    ////            if(logger->isDebugEnabled())
    ////            {
    ////                std::stringstream sout;
    ////                KGlobal->Print("k = ", sout,EMathematicaInput);
    ////                FGlobal.Print("r = ", sout,EMathematicaInput);
    ////                LOGPZ_DEBUG(logger,sout.str())
    ////            }
    ////    #endif
    //#endif
    
    std::cout << "Entering into Solve ..." << std::endl;
    an.Solve();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
#ifdef LOG4CXX
    TPZFMatrix<REAL> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
    solucao.Print("Sol = ",cout,EMathematicaInput);//imprime na formatacao do Mathematica
#endif
    
    std::cout << "Entering into Post processing ..." << std::endl;
    // Post processing
    int ndiv = 1;
    TPZManVector<std::string> scalarnames(5), vecnames(1);
    scalarnames[0] = "SigmaX";
    scalarnames[1] = "SigmaY";
    scalarnames[2] = "SigmaZ";
    scalarnames[3] = "TauXY";
    scalarnames[4] = "SolidPressure";
    vecnames[0] = "Displacement";
    //vecnames[1] = "";
    an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions2D.vtk");
    
    an.PostProcess(ndiv);
    //
    //
    //
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
    
    
}


//********************************** Cria malha Geometrica Circular (360 graus) *********************************************************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc) {
    
    
    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = M_PI;
    szmin = (Pi/2)*(rwb/ncirc)*(DrDcirc);
    
    // calcula comprimento radial da parede do poco ate contorno
    REAL radiallength;
    radiallength = re - rwb;
    
    // definindo variacao do angulo theta ao redor do poco
    // em rads!!!!!
    TPZVec<REAL> theta;
    theta.Resize(ncirc+1);
    REAL firsttheta = (2*Pi) / (ncirc);
    for (int k = 0; k<ncirc+1; k++) {
        REAL sumtheta = 0.;
        sumtheta += firsttheta * k;
        theta[k] = sumtheta;
    }
    
    
    //    // *******Imprime variacao dos angulos (em rads)
    //    std::cout<< "elementos de theta: " << endl;
    //    for (int t=0; t<ncirc+1; t++) {
    //        std::cout<< "Theta[" << t << "] :" << theta[t] << endl;
    //    }
    //    std::cout<< "Print theta " << endl;
    //    theta.Print();
    //    std::cout << endl;
    
    
    // nx = number of nodes in x direction
    // ny = number of nodes in y direction
    int nx,ny;
    nx = nrad+1;
    ny = ncirc+1;
    
    
    // Geometric Progression of the elements
    REAL q;
    if(nrad >1)
    {
        q = TPZGenGrid::GeometricProgression(szmin, radiallength, nrad);
    }
    else
    {
        q=radiallength;
    }
    // std::cout<< "valor de q " << q << endl; // imprime razao da PG
    
    
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    long i,j;
    long id, index;
    
    
    //*********************** Malha Circunferencial (360 graus) *********************//
    
    //vector to store a coordinate
    TPZVec <REAL> coord (2,0.);
    
    // aloca valor acumulado dos raios
    REAL rsum = 0.;
    REAL sz = szmin;
    
    //Nodes initialization
    for(i = 1; i < nx+1; i++){
        for(j = 1; j < ny+1; j++){
            
            // aloca coordenadas em cada no
            coord[0] = (rwb + rsum)* cos(theta[j-1]);
            coord[1] = (rwb + rsum)* sin(theta[j-1]);
            
            // id do elemento
            id = (i) * ny + (j);
            
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            
            
            //            // print
            //
            //            std::cout << "*****Iteracao nro: " << j << endl;
            //            std::cout << "rsum: " << rsum << endl;
            //            std::cout << "cos" << "[" << theta[j-1] << "]" <<": " ;
            //            std::cout << cos(theta[j-1]);
            //            std::cout << endl;
            //            std::cout << "sin" << "[" << theta[j-1] << "]" << ": ";
            //            std::cout << sin(theta[j-1]);
            //            std::cout << endl;
            //            std::cout << "Coord x: " << coord[0] << ";" << " Coord y: " << coord[1] << ";" << " Rad: " << theta[j-1] << endl;
            //            std::cout << endl;
            
            
        }
        
        rsum += sz; //valor acumulado dos raios
        sz *= q;
    }
    
    
    
    //vector to store a element connectivities
    TPZVec <long> connect(4,0);
    
    
    //Element connectivities
    for(i = 0; i < (nx -1); i++){
        for(j = 0; j < (ny -1); j++){
            // index = (i)*(ny - 1)+ (j);
            connect[0] = (i) * ny + (j);
            connect[1] = connect[0]+(ny);
            connect[2] = connect[1]+1;
            connect[3] = connect[0]+1;
            //Allocates and define the geometric element
            gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
            //std::cout << "connect: " << connect << endl;
            
            gmesh->ElementVec()[id];
            
            //std::cout << "id: " << id << endl;
            
        }
        
    }
    
    //*******Conecta os nos dos ultimos elementos da circunferencia com os primeiros *****************
    for (int k=0; k < nrad; k++) {
        TPZGeoEl *gel = gmesh->Element(((k+1)*(ncirc-1))+k);
        gel->SetNodeIndex(3, (nrad+1)*k);
        gel->SetNodeIndex(2, (nrad+1)*(k+1));
        
        //        gel->Print(); // verifica se a conexao esta correta
        //        std::cout << endl;
        
    }
    
    
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    
    
    //********** Criando Geo de BC ********//
    
    // bc = -1 -> Normal Pressure condition na parede do poco
    for (int i = 0; i<ncirc; i++ ) {
        // TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);
    }
    
    
    // bc = -2 -> StressField condition circunferencia externa do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -2);
    }
    
    
    // bc -3 -> Mixed, ponto fixo canto externo do poco (farfield) bottom
    TPZGeoElBC gbc1(gmesh->ElementVec()[(ncirc*nrad)-1],2,-3);
    
    // bc -4 -> Mixed, ponto fixo canto externo do poco (farfield) lateral direita
    TPZGeoElBC gbc2(gmesh->ElementVec()[(ncirc*nrad)-(ncirc/4)],1,-4);
    
    return gmesh;
    
}






//************************************ Cria malha Computacional para malha 360 graus ***********************************************************//

TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder)
{
    int matId = 1;
    const int dim = 2; //dimensao do problema
    
    
    //**************** Criando material  ********************************
    TPZMatElasticity2D *material = new TPZMatElasticity2D(matId);//criando material que implementa a formulacao fraca do problema modelo
    
    
    // Setting up paremeters
    //  copy this link http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
    REAL Eyoung = 29269e+6 , ni = 0.203, fbx = 0., fby = 0.;
    material->SetElasticity(Eyoung, ni, fbx, fby);
    
    
    /******* Calculating Inicial Stresses *******/
    // direction = direction/azimuth
    // inclination = wellbore inclination
    // problem assumption, inclined wellbore state = 1
    // Pwb = pressao da lama em Pa
    REAL Pi = M_PI;
    
    REAL direction = 0., inclination = 0.; //inicializa angulos
    direction   = 0.; // graus
    inclination = 0.; // graus
    
    // transforma graus em rad
    REAL directionT = 0.,inclinationT = 0.; // inicializa
    directionT = direction*(Pi/180); // rad
    inclinationT = inclination*(Pi/180); // rad
    
    // define disposicao do poco
    int inclinedwellbore = 1;
    
    // pressao da lama de perfuracao
    REAL Pwb = 29.9e+6; // Pa
    
    // Tensoes in Situ, horizontais e vertical em Pa
    REAL SigmaVV = 0., Sigmahh = 0., SigmaHH = 0.; // inicializa
    SigmaVV = -48.2e+6, Sigmahh = -45.9e+6, SigmaHH = -62.1e+6; //preenche
    // Seta os parametros do poco
    material->SetInclinedWellboreParameters(SigmaHH, Sigmahh, SigmaVV, directionT, inclinationT, inclinedwellbore);
    
    //Eh necessario chamar esse metodo para que sejam calculadas as tensoes iniciais apos a rotacao
    //Obtem tensor de tensoes iniciais
    REAL SigmaX = 0., SigmaXY = 0., SigmaY = 0., SigmaZ = 0.;
    material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ);
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        
        std::stringstream out;
        out << " Stress rotation " << std::endl;
        out << " SigmaX     = " << SigmaX << std::endl;
        out << " SigmaXY    = " << SigmaXY << std::endl;
        out << " SigmaY     = " << SigmaY << std::endl;
        out << " SigmaZ     = " << SigmaZ << std::endl;
        LOGPZ_DEBUG(logger,out.str())
    }
#endif
#endif
    
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
    cmesh->SetDimModel(dim);//seta dimensao do modelo
    
    // Inserindo material na malha
    cmesh->InsertMaterialObject(material);
    
    
    // cond contorno
    const int bc0 = -1, bc1 = -2, bc2 = -1, bc3 = -3, bc4 = -4; // ids igual da malha geometrica
    //bc2 = -3, bc3 = -4, bc4 = -5, bc5 = -6;
    const int normalpressure = 6, stressfield = 4, mixed = 2; // tipo de condicao de contorno
    //neumann = 1;
    //material->GetPreStress(SigmaX, SigmaXY, SigmaY, SigmaZ); // obtem tensoes iniciais
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    //    ///Inserir condicao de contorno parede do poco
    //    val1(0,0) = Pwb;
    //    val1(1,0) = 0;
    //    val1(0,1) = 0;
    //    val1(1,1) = Pwb;
    //    val2(0,0) = 0;
    //    val2(1,0) = 0;
    //    TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material
    
    ///Inserir condicao de contorno circunferencia externa
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond1 = material->CreateBC(material, bc1, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno circunferencia interna
    val1(0,0) = SigmaX;
    val1(1,0) = SigmaXY;
    val1(0,1) = SigmaXY;
    val1(1,1) = SigmaY;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond2 = material->CreateBC(material, bc2, stressfield, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto externo bottom
    val1(0,0) = 1.;
    val1(1,0) = 0.;
    val1(0,1) = 0.;
    val1(1,1) = 1;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond3 = material->CreateBC(material, bc3, mixed, val1, val2);//cria material
    
    ///Inserir condicao de contorno ponto externo lateral direita
    val1(0,0) = 1;
    val1(1,0) = 0;
    val1(0,1) = 0;
    val1(1,1) = 0;
    val2(0,0) = 0;
    val2(1,0) = 0;
    TPZMaterial * BCond4 = material->CreateBC(material, bc4, mixed, val1, val2);//cria material que implementa a condicao de contorno da parede do poco
    
    
    //cmesh->InsertMaterialObject(BCond0);//insere material na malha
    cmesh->InsertMaterialObject(BCond1);//insere material na malha
    cmesh->InsertMaterialObject(BCond2);//insere material na malha
    cmesh->InsertMaterialObject(BCond3);//insere material na malha
    cmesh->InsertMaterialObject(BCond4);//insere material na malha
    
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

