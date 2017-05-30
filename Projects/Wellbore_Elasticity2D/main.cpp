 #ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZMatElasticity2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"

#include "pzstepsolver.h"
#include "math.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZParFrontStructMatrix.h"

#include <cmath>
#include <set>

#include "pzelast3d.h"
#include <pzgengrid.h>

#include "pzrandomfield.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
//static LoggerPtr loggeradap(Logger::getLogger("pz.adaptivity"));
//static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
//static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;



TPZGeoMesh *GetMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc);

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, TPZFMatrix<REAL> SetKCorr);

TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc, REAL alpha, REAL beta, TPZFMatrix<REAL> GetKCorr);

TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder, TPZFMatrix<REAL> SetKCorr);

TPZCompMesh *CMesh3D(TPZGeoMesh *gmesh, int pOrder, bool Is3DQ);

TPZGeoMesh * ReadGeoMesh(std::string GridFileName, int dim);

//******* executa a simulação para malha 2D ***********//
int Problem2D();

//******* executa a simulação para malha 3D ***********//
int Problem3D();

/** @brief Compute approximation rates for Inclined wellbore analytic solution */
int ApproximationRates();

/** @brief Compute L2 error for a given computational mesh */
void ComputeErrorL2(TPZCompMesh * cmesh, REAL &error);

/** @brief compute uniform refinement */
void UniformRefinement(TPZGeoMesh *gmesh, int nh);

/** @brief compute inner product between tensors */
REAL Inner_Product(TPZFNMatrix<4,REAL> S, TPZFNMatrix<4,REAL> T);


//******* EXECUTA O CODIGO ***********//
int main(int argc, char *argv[])
{

//    Problem3D();
    Problem2D();
//
    //ApproximationRates();

    return 0;
}


/** @brief Compute approximation rates for Inclined wellbore analytic solution */
int ApproximationRates(){
    
    //******** Configura malha geometrica ***************/
    // rw = raio do poco (metros)
    // rext = raio externo do contorno (metros)
    // ncircle = nro elementos na parede do poco
    // nradial = nro de elementos da parede do poco ate o raio externo
    // drdcirc = proporcao do primeiro elemento
    REAL rw = 0.1;
    REAL rext = 4.0;
    int ncircle = 30;
    int nradial = 40;
    REAL drdcirc = 0.5;
    REAL Pi = M_PI;
    /************ Define Posicao do Poco **************/
    REAL direction = 0., inclination = 0.; //inicializa angulos
    direction   = 60.; // Azimuth em graus********
    inclination = 30.; // Polar Inclination em graus********
    
    // transforma graus em rad
    REAL alpha = 0., beta = 0.; // inicializa
    alpha = direction*(Pi/180); // rad
    beta = inclination*(Pi/180); // rad
    
    int numthreads = 2;
    int nh = 2;
    int np = 2;

    
    std::string plotfile = "ElasticitySolutions2D.vtk";
    TPZVec<REAL> errorvec;
    errorvec.Resize(nh);
    
    TPZVec<REAL> rates;
    rates.Resize(nh-1);
    
    int nelemtsr = nradial*ncircle;
    
    TPZFNMatrix<10,REAL> rates_array(np,nh-1,0.0);
    TPZFNMatrix<10,REAL> error_array(np,nh,0.0);
    TPZFMatrix<REAL> GetKCorr(nelemtsr,nelemtsr,0.0);
    
    int current_p = 0;
    
    for (int ip = 1; ip <= np; ip++) {
                current_p = ip;
        
        for (int ih = 0; ih < nh; ih++) {
            
            //******** geometric mesh ***************/
            TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc, alpha, beta, GetKCorr); //funcao para criar a malha GEOMETRICA de todo o poco
            const std::string nm("wellbore");
            gmesh->SetName(nm);
            
//#ifdef LOG4CXX
            std::ofstream outtxt("gmesh.txt"); //define arquivo de saida para impressao dos dados da malha
            gmesh->Print(outtxt);
            std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
//#endif
            
            
            //******** Apply geometric refinement ***************/
            UniformRefinement(gmesh, ih);
            
            //******** Configura malha Computacional ***************/
            
            TPZCompMesh *cmesh = CircularCMesh(gmesh, current_p, GetKCorr); //funcao para criar a malha COMPUTACIONAL de todo o poco
            TPZAnalysis an (cmesh);
            
            TPZSkylineStructMatrix strskyl(cmesh);
            strskyl.SetNumThreads(numthreads);
            an.SetStructuralMatrix(strskyl);
            TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
            direct->SetDirect(ECholesky);
            an.SetSolver(*direct);
            delete direct;
            direct = 0;
            
            std::cout << "number of dof = " << cmesh->NEquations() << std::endl;
            
            an.Assemble();
            an.Solve();
            std::cout << "problem solved. " << std::endl;
            
            cmesh->LoadSolution(an.Solution());
            REAL error = 0.0;
            ComputeErrorL2(cmesh, error);
            std::cout << "error computed. " << std::endl;
            
            std::cout << error << endl;
            
            errorvec[ih] = sqrt(error);
            error_array(ip-1,ih) = errorvec[ih];
            
        }
        
        for (int i = 0; i < nh - 1; i++) {
            if (i+1 > errorvec.size()) {
                continue;
            }
            
            rates[i]    = (log10(errorvec[i+1])-log10(errorvec[i]))/(log10(0.5));
            rates_array(ip-1,i) = rates[i];
            
        }
        
        
        
        std::cout << "P = " << current_p << " Error = " << errorvec << std::endl;
        std::cout << "P = " << current_p << " Rate = " << rates << std::endl;

        

        }
    
    std::ofstream out_error("error_file.txt");
    std::ofstream out_rates("rates_file.txt");
    error_array.Print("Errors = ",out_error,EMathematicaInput);
    rates_array.Print("Rates = ",out_rates,EMathematicaInput);
    
    return 0;
    
}


/** @brief Compute L2 error for a given computational mesh */
void ComputeErrorL2(TPZCompMesh * cmesh, REAL &error){
    
    int nel = cmesh->NElements();
    int dim =  cmesh->Dimension();
    int int_order = 10;
    int int_typ = 0;
    
    int sx_nu = 3;
    int sy_nu = 4;
    int sxy_nu = 6;
    
    int sx_an = 7;
    int sy_an = 8;
    int sxy_an = 10;
    
    TPZVec<REAL> sxnu;
    TPZVec<REAL> synu;
    TPZVec<REAL> sxynu;
    
    TPZVec<REAL> sxan;
    TPZVec<REAL> syan;
    TPZVec<REAL> sxyan;
    
    
    int defx_nu = 25;
    int defy_nu = 26;
    int defxy_nu = 27;
    
    int defx_an = 22;
    int defy_an = 23;
    int defxy_an = 24;
    
    TPZVec<REAL> defxnu;
    TPZVec<REAL> defynu;
    TPZVec<REAL> defxynu;
    
    TPZVec<REAL> defxan;
    TPZVec<REAL> defyan;
    TPZVec<REAL> defxyan;
    
    error = 0.0;
    

    TPZFNMatrix<4,REAL> snu(2,2,0.0);
    TPZFNMatrix<4,REAL> san(2,2,0.0);
       
    TPZFNMatrix<4,REAL> defnu(2,2,0.0);
    TPZFNMatrix<4,REAL> defan(2,2,0.0);
    
    TPZFNMatrix<4,REAL> esig(2,2,0.0);
    TPZFNMatrix<4,REAL> edef(2,2,0.0);
    
    
    for (int icel = 0; icel < nel; icel++) {
        TPZCompEl * cel = cmesh->Element(icel);
        
        
        if (!cel) {
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        
        if (!gel) {
            DebugStop();
        }
        
        
        /** Filtering boundary elements **/
        if (gel->Dimension() != dim) {
            continue;
        }
        
        
        /** Create integration rule **/
        int gel_side = gel->NSides() - 1;
        TPZIntPoints * IntegrationRule = gel->CreateSideIntegrationRule(gel_side, int_order);
        IntegrationRule->SetType(int_typ, int_order);
        int npoints = IntegrationRule->NPoints();
        
        TPZVec<REAL> duplet_xi_eta(2,0.0);
        REAL w;
        
        TPZFMatrix<REAL> jac,axes,jacinv;
        REAL detjac;
        
        REAL element_error = 0.0;
        
        REAL inner;
        
        for (int i = 0; i < npoints; i++) {
            gel->Jacobian(duplet_xi_eta, jac, axes, detjac, jacinv);
            IntegrationRule->Point(i, duplet_xi_eta, w);
            
            /** Computing analytic stress **/
            cel->Solution(duplet_xi_eta, sx_an, sxan);
            cel->Solution(duplet_xi_eta, sy_an, syan);
            cel->Solution(duplet_xi_eta, sxy_an, sxyan);
            
            /** Computing numeric stress **/
            cel->Solution(duplet_xi_eta, sx_nu, sxnu);
            cel->Solution(duplet_xi_eta, sy_nu, synu);
            cel->Solution(duplet_xi_eta, sxy_nu, sxynu);
            
            san(0,0) = sxan[0];
            san(1,1) = syan[0];
            san(1,0) = sxyan[0];
            san(0,1) = sxyan[0];
            
            snu(0,0) = sxnu[0];
            snu(1,1) = synu[0];
            snu(1,0) = sxynu[0];
            snu(0,1) = sxynu[0];
           
            
            /** Computing analytic deformation **/
            cel->Solution(duplet_xi_eta, defx_an, defxan);
            cel->Solution(duplet_xi_eta, defy_an, defyan);
            cel->Solution(duplet_xi_eta, defxy_an, defxyan);
            
            /** Computing numeric deformation **/
            cel->Solution(duplet_xi_eta, defx_nu, defxnu);
            cel->Solution(duplet_xi_eta, defy_nu, defynu);
            cel->Solution(duplet_xi_eta, defxy_nu, defxynu);
            
            defan(0,0) = defxan[0];
            defan(1,1) = defyan[0];
            defan(1,0) = defxyan[0];
            defan(0,1) = defxyan[0];
            
            defnu(0,0) = defxnu[0];
            defnu(1,1) = defynu[0];
            defnu(1,0) = defxynu[0];
            defnu(0,1) = defxynu[0];
        
            
            esig = san - snu;
            edef = defan - defnu;
            
            inner = Inner_Product(edef,esig);
            element_error += w * detjac * inner;
        }
        
        error += element_error;
    }
    
    std::cout << "error = " << sqrt(error) << std::endl;
    
    
}

REAL Inner_Product(TPZFNMatrix<4,REAL> S , TPZFNMatrix<4,REAL> T){
    
    if ((S.Rows() != T.Rows()) && (S.Cols() != T.Cols())) {
        DebugStop();
    }
    
    
    //** Tr[Transpose[s].t] = S00 T00 + S01 T01 + S10 T10 + S11 T11 */
    REAL inner_product = S(0,0) * T(0,0) + S(0,1) * T(0,1) + S(1,0) * T(1,0) + S(1,1) * T(1,1);
    return inner_product;
}




//******************* Refinamento *****************************//
void UniformRefinement(TPZGeoMesh *gmesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
    gmesh->BuildConnectivity();
}





//*********************************************** PROBLEMA 2D ******************************************************//
int Problem2D(){
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    //******** Configura malha geometrica ***************/
    // rw = raio do poco (metros)
    // rext = raio externo do contorno (metros)
    // ncircle = nro elementos na parede do poco
    // nradial = nro de elementos da parede do poco ate o raio externo
    // drdcirc = proporcao do primeiro elemento
    
    //REAL rw = 0.10795;
    REAL rw = 0.1;
    
    REAL rext = 3.0;
    int ncircle = 30;
    int nradial = 25;
    REAL drdcirc = 0.5;
    
    REAL Pi = M_PI;
    /************ Define Posicao do Poco **************/
    REAL direction = 0., inclination = 0.; //inicializa angulos
    direction   = 0.; // Azimuth em graus******** 30
    inclination = 0.; // Polar Inclination em graus******** 50
    
    // transforma graus em rad
    REAL alpha = 0., beta = 0.; // inicializa
    alpha = direction*(Pi/180); // rad
    beta = inclination*(Pi/180); // rad
    
    int nelemtsr = nradial*ncircle;
 
    TPZFMatrix<REAL> GetKCorr(nelemtsr,nelemtsr,0.0);
    
    TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc, alpha, beta, GetKCorr); //funcao para criar a malha GEOMETRICA de todo o poco
    //TPZGeoMesh *gmesh = GetMesh(rw, rext, ncircle, nradial, drdcirc); //funcao para criar a malha GEOMETRICA de 1/4 do poco
    
    
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
    
    TPZFMatrix<REAL> SetKCorr(nelemtsr,nelemtsr,0.0);
    SetKCorr = GetKCorr;
    
    TPZCompMesh *cmesh = CircularCMesh(gmesh, p, SetKCorr); //funcao para criar a malha COMPUTACIONAL de todo o poco
    //TPZCompMesh *cmesh = CMesh(gmesh, p); //funcao para criar a malha COMPUTACIONAL de 1/4 do poco
    
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
    }
    else{
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
    
    if (projection==1) {
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
    //
    //
    //
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
    
    
}



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

/******************************************************** FIM do MAIN *********************************************************/


//********************************** Cria malha Geometrica Circular (360 graus) *********************************************************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *CircularGeoMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc, REAL alpha, REAL beta, TPZFMatrix<REAL> GetKCorr) {
    
    
    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = M_PI;
    szmin = (2*Pi)*(rwb/ncirc)*(DrDcirc); //NNNNNNNNN
    
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
    
  //  q = 1; // CALCULAR TAXAS DE CONVERGENCIA, malha com mesmo tamanho de elementos
//     std::cout<< "valor de q " << q << endl; // imprime razao da PG
    
   
    //Creates the geometric mesh... The nodes and elements
    //will be inserted into mesh object during initilize process
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    long i,j;
    long id, index;
    
    
    //*********************** Malha Circunferencial (360 graus) *********************//
    
    //vector to store a coordinate
    TPZVec <REAL> coord (3,0.);
    TPZVec <REAL> coordT (3,0.);
    
    
    // aloca valor acumulado dos raios
    REAL rsum = 0.;
    REAL sz = 0.;
    
    if (q == 1) {
    sz = radiallength/nrad;
    }
    else {
    
     sz = szmin; // para taxas
    
    }
    
    
    //Nodes initialization
    for(i = 1; i < nx+1; i++){
        for(j = 1; j < ny+1; j++){
            
            // aloca coordenadas em cada no
            coord[0] = (rwb + rsum)* cos(theta[j-1]);
            coord[1] = (rwb + rsum)* sin(theta[j-1]);
            coord[2] = 0.;
            
            
                // id do elemento
                id = (i) * ny + (j);
                
                //Get the index in the mesh nodes vector for the new node
                index = gmesh->NodeVec().AllocateNewElement();
                //Set the value of the node in the mesh nodes vector
                gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
            
            
////               Print
//            std::cout << "*****Iteracao nro: " << j << endl;
//            std::cout << "rsum: " << rsum << endl;
//            std::cout << "cos" << "[" << theta[j-1] << "]" <<": " ;
//            std::cout << cos(theta[j-1]);
//            std::cout << endl;
//            std::cout << "sin" << "[" << theta[j-1] << "]" << ": ";
//            std::cout << sin(theta[j-1]);
//            std::cout << endl;
//            
//            std::cout << "Coord x: " << coord[0] << ";" << " Coord y: " << coord[1] << ";" << " Coord z: " << coord[2] << endl;
//            std::cout << endl;
//            
//            std::cout << "alpha: " << alpha << ";" << " beta: " << beta  << endl;
//            std::cout << "CoordT x: " << coordT[0] << ";" << " CoordT y: " << coordT[1] << ";" << " CoordT z: " << coordT[2] << ";" << endl;
//            std::cout << endl;
            
            
            
        }
        
        

        
        if (q == 1) {
             rsum = sz  * i * q; // para taxas
        }
        else {
            
            rsum += sz; //valor acumulado dos raios
            sz *= q;
            
        }

//        std::cout << "sz: " << sz << endl;
//
////        
//        std::cout << "rsum: " << rsum << endl;
        

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
//            std::cout << "connect: " << connect << endl;
            
            gmesh->ElementVec()[id];
            
//            std::cout << "id: " << id << endl;
            
        }
        
    }
    
    //*******Conecta os nos dos ultimos elementos da circunferencia com os primeiros *****************
        for (int k=0; k < nrad; k++) {
            TPZGeoEl *gel = gmesh->Element(((k+1)*(ncirc-1))+k);
            
            //gel->Print(); // verifica
            
            TPZGeoEl *gelAbove = gmesh->Element((((k+1)*(ncirc-1))+k)-(ncirc-1));
            TPZVec <int> gelAboveIndex(4,0);
            gelAbove->GetNodeIndices(gelAboveIndex);
            
            //gelAbove->Print(); // verifica
            
            gel->SetNodeIndex(3, gelAboveIndex[0]);
            gel->SetNodeIndex(2, gelAboveIndex[1]);
            
            //gel->Print(); // verifica se a conexao esta correta
        
           //std::cout << endl;
        
        }
   
    
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    //gmesh->Print();
    
    
    
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
    
    
    
    
    
    //****** Cria matriz da norma entre os centroides (para a matriz de correlacao) *****//
    
    TPZManVector<REAL> centerpsi(3), center(3);
    // Refinamento de elementos selecionados
    TPZGeoEl *gel;
    TPZVec<TPZGeoEl *> sub;
    
    //int nelem = 0;
    //int ngelem=gmesh->NElements();
    
    REAL nelemtsr = nrad*ncirc; //Nro de elementos - quadrilateros
    
    TPZFMatrix<REAL> CenterNorm(nelemtsr,nelemtsr,0.0);
    
    TPZManVector<REAL,3> CenterPoint1;
    TPZManVector<REAL,3> CenterPoint2;
    
    // Matriz da distancia entre os centroides
    for (i = 0; i < nelemtsr; i++) {
          for (j = 0; j < nelemtsr; j++) {
        
        gel = gmesh->ElementVec()[i];
        
        gel->CenterPoint(8,centerpsi);
        
        gel->X(centerpsi,center);
            
        CenterPoint1 = center;
         
        gel = gmesh->ElementVec()[j];
              
        gel->CenterPoint(8,centerpsi);
              
        gel->X(centerpsi,center);
            
        CenterPoint2 = center;
              
              CenterNorm(i,j) = sqrt(pow((CenterPoint2[0]-CenterPoint1[0]),2)+pow((CenterPoint2[1]-CenterPoint1[1]),2)+pow((CenterPoint2[2]-CenterPoint1[2]),2));
        }
        
    }
    
    
    // Matriz de correlacao
    
    TPZFMatrix<REAL> KCorr(nelemtsr,nelemtsr,0.0);
    REAL r = 0.;        // define distancia r entre cada centroide
    REAL e = M_E;       //Numero de Euler
    REAL scale = 0.;    //Valor de alpha, escala normalizada
    
    for (i = 0; i < nelemtsr; i++) {
        for (j = 0; j < nelemtsr; j++) {
            
            r = CenterNorm(i,j);
            
            KCorr(i,j) = pow(e, (-scale*(pow(r, 2))));
            
        }
    }
    
    GetKCorr = KCorr;
    
    
  //  std::cout << center[0] << ";" << center[1] << ";" << center[2] << ";" << std::endl;
    
   
    return gmesh;


    
}






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
        
     
       ///Inserir condicao de contorno parede do poco
       val1(0,0) = Pwb;
       val1(1,1) = Pwb;
       val1(2,2) = Pwb;
       TPZMaterial * BCond0 = material->CreateBC(material, bc0, normalpressure, val1, val2);//cria material
        
       ///Inserir condicao de contorno circunferencia externa
       val1(0,0) = SigmaX;
       val1(1,0) = SigmaXY;
       val1(0,1) = SigmaXY;
       val1(1,1) = SigmaY;
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
    
    else  {
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
    else{
        
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
    

    

    

    
}


/*************************************************** Le malha do GID ****************************************************************/

TPZGeoMesh * ReadGeoMesh(std::string GridFileName, int dim)
{
    TPZReadGIDGrid GeometryInfo;
    REAL s = 1.0;
    GeometryInfo.SetfDimensionlessL(s);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    gmesh->SetDimension(dim);
    return gmesh;
}





// ******************************************** Cria malha Geometrica 1/4 do Poco *************************************************************/

// Para gerar a malha eh necessario:
// rwb -> raio do poco
// re -> raio do contorno
// ncirc -> nro elem ao longo 1/4 do poco
// nrad -> nro elem da parede do poco ate contorno externo
// DrDcirc -> proporcao dos elementos da parede do poco

TPZGeoMesh *GetMesh (REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc) {


    // calcula comprimento radial do primeiro elemento
    REAL szmin;
    REAL Pi = 3.14159;
    szmin = (Pi/2)*(rwb/ncirc)*(DrDcirc);

    // calcula comprimento radial da parede do poco ate contorno
    REAL radiallength;
    radiallength = re - rwb;

    // definindo variacao do angulo theta ao redor do poco
    // em rads!!!!!
    TPZVec<REAL> theta;
    theta.Resize(ncirc+1);
    
    REAL firsttheta = ((Pi/2)/5) / (ncirc);
    for (int k = 0; k<ncirc+1; k++) {
        REAL sumtheta = 0.;
        sumtheta += firsttheta * k;
        theta[k] = sumtheta;
    }

//
//    // *******Imprime variacao dos angulos (em rads)
//    std::cout<< "elementos de theta: " << endl;
//        for (int t=0; t<ncirc+1; t++) {
//            std::cout<< "Theta[" << t << "] :" << theta[t] << endl;
//            }
//            std::cout<< "Print theta " << endl;
//            theta.Print();
//            std::cout << endl;


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


//*********************** Malha Circunferencial (1/4) *********************//

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
            
            //std::cout << "id: " << id << endl;

            gmesh->ElementVec()[id];
        }
    }
    //Generate neighborhod information
    gmesh->BuildConnectivity();



    //********** Criando Geo de BC,********//

    // bc = -1 -> Normal Pressure condition
    for (int i = 0; i<ncirc; i++ ) {
//        TPZGeoElBC gbc(gmesh->ElementVec()[i],7,-1);
        gmesh->ElementVec()[i]->CreateBCGeoEl(7, -1);

    }

    // bc = -2 -> Neumann condition contorno externo bottom
    for (int i = 0; i<nrad; i++ ) {
        gmesh->ElementVec()[ncirc*i]->CreateBCGeoEl(4, -2);

    }

    // bc = -3 -> Neumann condition contorno externo upper
    for (int i = 1; i<nrad+1; i++ ) {
        gmesh->ElementVec()[(ncirc*i)-1]->CreateBCGeoEl(6, -3);

    }

    // bc = -4 -> Neumann condition arco externo do farfield
    for (int i = 1; i<ncirc+1; i++ ) {
        gmesh->ElementVec()[(ncirc*nrad)-i]->CreateBCGeoEl(5, -4);

    }


    // bc -5 -> Mixed, ponto fixo canto interno parede do poco bottom
    TPZGeoElBC gbc1(gmesh->ElementVec()[0],0,-5);

    // bc -6 -> Mixed, ponto fixo externo ao 1/4 de poco bottom
    TPZGeoElBC gbc2(gmesh->ElementVec()[0],1,-6);
    
    // bc -5 -> Mixed, ponto fixo canto interno parede do poco bottom
    TPZGeoElBC gbc3(gmesh->ElementVec()[0],2,-7);
    
    // bc -6 -> Mixed, ponto fixo externo ao 1/4 de poco bottom
    TPZGeoElBC gbc4(gmesh->ElementVec()[0],3,-8);
    

//    {
//        std::ofstream malha("../malhageo.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malha);
//    }

    return gmesh;

}



// ********************************************* Cria malha Computacional para 1/4 do Poco *****************************************************/


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
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






