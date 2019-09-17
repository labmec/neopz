//
//  AproximationRates.cpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#include "AproximationRates.hpp"

using namespace std;

/** @brief Compute approximation rates for Inclined wellbore analytic solution */
int ApproximationRates(TPZFMatrix<STATE> &M){
    
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
    
    int nSquareElements = nradial * ncircle;
    
    TPZFNMatrix<10,REAL> rates_array(np,nh-1,0.0);
    TPZFNMatrix<10,REAL> error_array(np,nh,0.0);
    
    int current_p = 0;
    
    for (int ip = 1; ip <= np; ip++) {
        current_p = ip;
        
        for (int ih = 0; ih < nh; ih++) {
            
            //******** geometric mesh ***************/
            TPZGeoMesh *gmesh = CircularGeoMesh (rw, rext, ncircle, nradial, drdcirc, alpha, beta); //funcao para criar a malha GEOMETRICA de todo o poco
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
            
            // Cria a malha COMPUTACIONAL de todo o poco
            int projection = 0;
            int inclinedwellbore = 0;
            REAL Pwb = 0;
            int analytic = 0;
            REAL SigmaV = 0, Sigmah, SigmaH;
            bool isStochastic = false;
            REAL scale = 0.5; // exponential scale
            int funcE = 1;
            int funcnu = 1;
            int distribE = 1;
            int distribnu = 1;
            
            TPZCompMesh *cmesh = CircularCMesh(gmesh, current_p, projection, inclinedwellbore,
                                               analytic, SigmaV, Sigmah, SigmaH, Pwb, rw, rext,
                                               direction, inclination, isStochastic, nSquareElements,
                                               M,  scale, funcE, funcnu, distribE, distribnu);
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

