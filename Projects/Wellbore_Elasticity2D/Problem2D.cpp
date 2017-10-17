//
//  Problem2D.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "Problem2D.hpp"

// Configura malha geometrica
// rw = raio do poco (metros)
// rext = raio externo do contorno (metros)
// ncircle = nro elementos na parede do poco
// nradial = nro de elementos da parede do poco ate o raio externo
// drdcirc = proporcao do primeiro elemento

int Problem2D(REAL rw, REAL rext, int ncircle, int nradial, int projection, int inclinedwellbore,
              int analytic, REAL SigmaV, REAL Sigmah, REAL SigmaH, REAL Pwb, REAL drdcirc,
              REAL direction, REAL inclination, bool isStochastic) {

#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    std::string dirname = PZSOURCEDIR;
    
    // Transforma graus em rad
    REAL alpha = direction * (M_PI / 180);
    REAL beta = inclination * (M_PI / 180);
    
    int nSquareElements = nradial * ncircle;
    
    TPZFMatrix<REAL> GetKCorr(nSquareElements, nSquareElements, 0.0);
    
    // Cria a malha GEOMETRICA de todo o poco
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc, alpha, beta, GetKCorr);
    
    // Cria a malha GEOMETRICA de 1/4 do poco
    //TPZGeoMesh *gmesh = GetMesh(rw, rext, ncircle, nradial, drdcirc);
    
    const std::string nm("line");
    gmesh->SetName(nm);
    
#ifdef LOG4CXX
    std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
    gmesh->Print(outtxt);
    std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
#endif
    
    // Configura malha Computacional
    int p = 2;
    TPZCompEl::SetgOrder(p);
    
    // Cria a malha COMPUTACIONAL de todo o poco
    TPZCompMesh *cmesh = CircularCMesh(gmesh, p, projection, inclinedwellbore, analytic, SigmaV,
                                       Sigmah, SigmaH, Pwb, rw, direction, inclination,
                                       isStochastic, nSquareElements);
    
    // Cria a malha COMPUTACIONAL de 1/4 do poco
    //TPZCompMesh *cmesh = CMesh(gmesh, p);
    
    // Solving linear equations
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
    }
    else {
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
    
    std::cout << "Entering into Solve ..." << std::endl;
    
    // Assembla atriz de rigidez global e o vetor de carga e inverte o sist. de equacoes
    an.Solve();
    
    std::cout << "Entering into Post processing ..." << std::endl;
    
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
        an.DefineGraphMesh(2, scalarnames, vecnames, "ElasticitySolutions2D.vtk");
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
        scalarnames.Push("Gaussian_Field_E");
        
        //vecnames[1] = "";
        an.DefineGraphMesh(2, scalarnames, vecnames, "ElasticitySolutions2D.vtk");
    }
    
    an.PostProcess(NDIV);
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}
