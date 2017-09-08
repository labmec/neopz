//
//  Problem2D.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "Problem2D.hpp"

//******** Configura malha geometrica ***************/
// rw = raio do poco (metros)
// rext = raio externo do contorno (metros)
// ncircle = nro elementos na parede do poco
// nradial = nro de elementos da parede do poco ate o raio externo
// drdcirc = proporcao do primeiro elemento

int Problem2D(REAL rw, REAL rext, int ncircle, int nradial, REAL drdcirc,
              REAL direction, REAL inclination, bool isStochastic) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    std::string dirname = PZSOURCEDIR;
    
    REAL Pi = M_PI;
    
    // transforma graus em rad
    REAL alpha = direction * (Pi/180);
    REAL beta = inclination * (Pi/180);
    
    int nelemtsr = nradial * ncircle;
    
    TPZFMatrix<REAL> GetKCorr(nelemtsr,nelemtsr,0.0);
    
    //funcao para criar a malha GEOMETRICA de todo o poco
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc,
                                         alpha, beta, GetKCorr);
    
    //funcao para criar a malha GEOMETRICA de 1/4 do poco
    //TPZGeoMesh *gmesh = GetMesh(rw, rext, ncircle, nradial, drdcirc);
    
    const std::string nm("line");
    gmesh->SetName(nm);
    
#ifdef LOG4CXX
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
#endif
    
    //******** Configura malha Computacional ***************/
    
    int p = 2;
    TPZCompEl::SetgOrder(p);
    
    TPZFMatrix<REAL> SetKCorr(nelemtsr, nelemtsr, 0.0);
    SetKCorr = GetKCorr;
    
    // Cria a malha COMPUTACIONAL de todo o poco
    TPZCompMesh *cmesh = NULL;
    
    if(isStochastic) {
        cmesh = CircularCMesh(gmesh, p, SetKCorr);
    }
    else {
        // TODO - implement CircularCMesh without correlation matrix
        // cmesh = CircularCMesh(gmesh, p);
    }
    
    // Cria a malha COMPUTACIONAL de 1/4 do poco
    //TPZCompMesh *cmesh = CMesh(gmesh, p);
    
    // Solving linear equations
    // Initial steps
    TPZAnalysis an (cmesh);
    int numthreads = 2;
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
    } else {
        TPZSkylineStructMatrix strskyl(cmesh);
        strskyl.SetNumThreads(0);
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
    
    //      an.Rhs() ;
    
    //        TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    //        TPZFMatrix<STATE> FGlobal;
    //        KGlobal =   an.Solver().Matrix();
    //        FGlobal =   an.Rhs();
    //
    //    #ifdef PZDEBUG
    //        #ifdef LOG4CXX
    //                if(logger->isDebugEnabled())
    //                {
    //                    std::stringstream sout;
    //                    KGlobal->Print("k = ", sout,EMathematicaInput);
    //                    FGlobal.Print("r = ", sout,EMathematicaInput);
    //                    LOGPZ_DEBUG(logger,sout.str())
    //                }
    //        #endif
    //    #endif
    
    //
    //    std::cout << "Rhs ..." << std::endl;
    //
    //#ifdef LOG4CXX
    //    TPZFMatrix<REAL> FGlobal = an.Rhs();
    //    FGlobal.Print("Rhs = ",cout,EMathematicaInput);
    //#endif
    //
    //    std::cout << std::endl;
    //
    //
    std::cout << "Entering into Solve ..." << std::endl;
    an.Solve();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    
    //#ifdef LOG4CXX
    //    TPZFMatrix<REAL> solucao=cmesh->Solution(); //Pegando o vetor de solucao, alphaj
    //
    ////    std::ofstream fileAlpha("alpha.txt");
    ////    an.Solution().Print("Alpha = ", fileAlpha, EMathematicaInput);
    ////
    //    solucao.Print("Sol = ",cout,EMathematicaInput);//imprime na formatacao do Mathematica
    //#endif
    
    std::cout << "Entering into Post processing ..." << std::endl;
    // Post processing
    int ndiv = 2;
    int projection = 0; // define se havera projecao no plano horizontal
    
    if (projection == 1) {
        TPZStack<std::string> scalarnames, vecnames;
        scalarnames.Push("SigmaXProjected");
        scalarnames.Push("SigmaYProjected");
        scalarnames.Push("SigmaZProjected");
        scalarnames.Push("TauXYProjected");
        scalarnames.Push("SolidPressureProjected");
        scalarnames.Push("SigmaXAnalyticProjected");
        scalarnames.Push("SigmaYAnalyticProjected");
        scalarnames.Push("SigmaZAnalyticProjected");
        scalarnames.Push("TauXYAnalyticProjected");
        scalarnames.Push("SolidPressureAnalyticProjected");
        scalarnames.Push("J2_Projected");
        scalarnames.Push("F1_Projected");
        scalarnames.Push("I1_Projected");
        //vecnames[1] = "";
        an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions2D.vtk");
    } else {
        TPZStack<std::string> scalarnames, vecnames;
        scalarnames.Push("SigmaX");
        scalarnames.Push("SigmaY");
        scalarnames.Push("SigmaZ");
        scalarnames.Push("TauXY");
        scalarnames.Push("SolidPressure");
        scalarnames.Push("SigmaXAnalytic");
        scalarnames.Push("SigmaYAnalytic");
        scalarnames.Push("SigmaZAnalytic");
        scalarnames.Push("TauXYAnalytic");
        scalarnames.Push("SolidPressureAnalytic");
        vecnames.Push("Displacement");
        scalarnames.Push("ExxAnalytic");
        scalarnames.Push("EyyAnalytic");
        scalarnames.Push("ExyAnalytic");
        scalarnames.Push("Exx");
        scalarnames.Push("Eyy");
        scalarnames.Push("Exy");
        scalarnames.Push("SigmaXProjected");
        scalarnames.Push("SigmaYProjected");
        scalarnames.Push("SigmaZProjected");
        scalarnames.Push("TauXYProjected");
        scalarnames.Push("SolidPressureProjected");
        scalarnames.Push("J2");
        scalarnames.Push("F1");
        scalarnames.Push("I1");
        scalarnames.Push("Sigma1");
        scalarnames.Push("Sigma2");
        scalarnames.Push("Sigma3");
        scalarnames.Push("CheckingVM1");
        scalarnames.Push("CheckingVM2");
        scalarnames.Push("CheckingVM3");
        scalarnames.Push("F_Mogi-Coulomb");
        //vecnames[1] = "";
        an.DefineGraphMesh(2,scalarnames,vecnames,"ElasticitySolutions2D.vtk");
        
    }
    
    an.PostProcess(ndiv);
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}
