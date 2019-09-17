//
//  Problem2D.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "Problem2D.hpp"
//#include <__config>
#include <ios>
#include <streambuf>
#include <istream>
#include <ostream>



// Configura malha geometrica
// rw = raio do poco (metros)
// rext = raio externo do contorno (metros)
// ncircle = nro elementos na parede do poco
// nradial = nro de elementos da parede do poco ate o raio externo
// drdcirc = proporcao do primeiro elemento

void PrintSolution(std::ofstream &solutionfile,int &icase,TPZGeoMesh *gmesh);

int Problem2D(REAL rw, REAL rext, int ncircle, int nradial, int projection, int inclinedwellbore,
              int analytic, REAL SigmaV, REAL Sigmah, REAL SigmaH, REAL Pwb, REAL drdcirc,
              REAL direction, REAL inclination, bool isStochastic,std::ofstream &solutionfile,
              int &icase, TPZFMatrix<STATE> &M, REAL scale, int funcE, int funcnu, int distribE, int distribnu) {
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    std::string dirname = PZSOURCEDIR;
    
    // Transforma graus em rad
    REAL alpha = direction * (M_PI / 180);
    REAL beta = inclination * (M_PI / 180);
    
    int nSquareElements = nradial * ncircle;
    
    // Cria a malha GEOMETRICA de todo o poco
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc, alpha, beta);
    
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
                                       Sigmah, SigmaH, Pwb, rw, rext, direction, inclination,
                                       isStochastic, nSquareElements, M,  scale, funcE, funcnu, distribE, distribnu);
    
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
    
    int stochasticanalysis = 1;
    
    std::string icaseStr = std::to_string(icase + 1);
    std::string elasticitySolFilename("../simulacao_Journal/" + icaseStr + "_Vertical_Stoch_Pw19_5_2m.vtk");
    
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
        an.DefineGraphMesh(2, scalarnames, vecnames, elasticitySolFilename);
    }
    else if (stochasticanalysis == 1){
        TPZStack<std::string> scalarnames, vecnames;
        scalarnames.Push("Plot_F1");
        scalarnames.Push("Stochastic_Field_E");
        scalarnames.Push("Tensile_Fail");
        
        an.DefineGraphMesh(2, scalarnames, vecnames, elasticitySolFilename);
    }
    else {
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
        scalarnames.Push("J2");
        scalarnames.Push("F1");
        scalarnames.Push("I1");
        scalarnames.Push("Sigma1");
        scalarnames.Push("Sigma2");
        scalarnames.Push("Sigma3");
        scalarnames.Push("Stochastic_Field_E");
        //vecnames[1] = "";
        an.DefineGraphMesh(2, scalarnames, vecnames, elasticitySolFilename);
    }
    
    an.PostProcess(NDIV);
    
    //std::ofstream solutionfile("f1_solution.csv");
    PrintSolution(solutionfile,icase,gmesh);
    
    //Cleanup
    gmesh->ResetReference();
    if(cmesh) delete cmesh;
    if(gmesh) delete gmesh;
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}

void PrintSolution(std::ofstream &solutionfile,int &icase,TPZGeoMesh *gmesh) {
    
    // Verify if the file is open
    if(!solutionfile.is_open()) DebugStop();
    
    //Intermediaries
    TPZGeoEl* geoel     = NULL;
    int matid           = MATERIAL_ID;  // COLOQUE SEU MATERIAL ID
    int var             = 42; // COLOQUE O NUMERO DA SOLUCAO (VARIAVEL DE ESTADO)
    REAL tol            = 1e-8; // talvez esse valor possa ser alterado (aumentar?)
    REAL qsivalue       = 0;
    REAL etavalue       = 0;
    REAL deltaqsi       = 0;
    REAL geoel_area     = 0;
    REAL geoelplast_area= 0;//geoel plastified area
    REAL totalplast_area= 0;//total plastified area
    int ndiv            = 20;//> 0
    int side            = 8;//side 8, to compute the area
    int counter         = 0;
    int ntotal          = 0;
    std::set<long> 	nodeindex;
    TPZManVector<REAL,3> x(3);
    TPZVec<REAL> qsi(2,0);
    TPZVec<STATE> sol;
    
    //check tensile failure
    int tensileFail     = 0;     // each integration point
    int var2            = 43;
    TPZVec<STATE> sol2;
    int totalFailure    = 0;    // sum when tensile failure occur
    int fail = 0;               // print failure status
    
    
    TPZFMatrix<REAL> jac;       //it is not being used
    TPZFMatrix<REAL> axes;      //it is not being used
    REAL master_el_area = 4.;   // area of the master element
    REAL detjac;				//this is important
    REAL weight;				//this is important
    TPZFMatrix<REAL> jacinv;    //it is not being used
    
    deltaqsi = 1.0/ndiv;
    ntotal   = ndiv*ndiv;
    totalplast_area = 0;
    
    for(long i = 0; i < gmesh->NElements(); i++){
        
        geoel = gmesh->ElementVec()[i];
        
        if(!geoel) continue;
        if(geoel->HasSubElement()) continue;
        if(geoel->MaterialId() != matid) continue;
        if(!geoel->Reference()) DebugStop(); //Why did this element lost its comp element?
        
        weight = master_el_area*(1./ntotal);//same weight for each qsi-eta
        qsivalue = -1;
        counter  = 0;
        geoelplast_area = 0.;
        geoel_area = 0.;
        
        for(int j=0;j<ndiv;j++){//qsi
            qsivalue += deltaqsi;
            etavalue  = -1;
            
            for(int k=0;k<ndiv;k++){//eta
                etavalue += deltaqsi;
                
                //fill qsi vector
                qsi[0] = qsivalue;
                qsi[1] = etavalue;
                
                //now, compute solution
                sol.clear();
                geoel->Reference()->Solution(qsi,var,sol);
                
                //original
                if (sol[0]>0) counter++;
                etavalue += deltaqsi;
                
                //now, compute detjac for this qsi-eta coordinate
                geoel->Jacobian(qsi,jac,axes,detjac,jacinv);
                geoel_area += detjac*weight;
                if(sol[0]>0) {
                    geoelplast_area += detjac*weight;
                }
                
                //tensile failure
                geoel->Reference()->Solution(qsi,var2,sol2);
                if (sol2[0]==1) {
                    tensileFail += 1;
                }
                
            }
            qsivalue += deltaqsi;
        }
        
        //sum
        totalplast_area += geoelplast_area;
        
        totalFailure    += tensileFail;
        
    } //loop over elements
    
    if (totalFailure>0) {
        fail = 1;
    }
    
    std::cout << "Case " << icase+1 << " total plastified area " << totalplast_area <<
    " Tensile Failure = " << fail << std::endl;
    
    solutionfile << icase+1 <<","<< totalplast_area <<","<< fail << std::endl;
    
    solutionfile.flush();
}

