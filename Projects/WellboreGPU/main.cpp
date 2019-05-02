
#include "path.h"
#include <iostream>
#include <string.h>
#include <ctime>
#include <algorithm>
#include <iterator>

// Neopz
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzinterpolationspace.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "tpzintpoints.h"
#include "TPZMatElasticity2D.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "pzpostprocanalysis.h"
#include "pzfstrmatrix.h"

#include "TPZMatElastoPlastic2D.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElastoPlasticMem.h"
#include "TPZElasticCriterion.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZBndCondWithMem.h"
#include "TPZIntPointsFEM.h"

#include "TElastoPlasticData.h"
#include "TRKSolution.h"

/// Gmsh mesh
TPZGeoMesh * ReadGeometry(std::string geometry_file);
void PrintGeometry(TPZGeoMesh * geometry);

/// CompMesh elastoplasticity
TPZCompMesh *CmeshElastoplasticity(TPZGeoMesh *gmesh, int pOrder, TElastoPlasticData & material_data);

/// Material configuration
TElastoPlasticData WellboreConfig();

/// Material configuration for RK verification
TElastoPlasticData WellboreConfigRK();

/// Solve using assemble residual of all intg points at once
void SolutionAllPoints(TPZAnalysis * analysis, int n_iterations, REAL tolerance, TElastoPlasticData & wellbore_material);

///Set Analysis
TPZAnalysis * Analysis(TPZCompMesh * cmesh, int n_threads);

///Solve using Newton method
void Solution(TPZAnalysis *analysis, int n_iterations, REAL tolerance);

/// Accept solution
void AcceptPseudoTimeStepSolution(TPZAnalysis * an, TPZCompMesh * cmesh);

/// Print mesh memory data
void PrintMemory(TPZCompMesh * cmesh);

/// Post process
void PostProcess(TPZCompMesh *cmesh, TElastoPlasticData material, int n_threads, std::string vtk_file);

///RK Approximation
void RKApproximation (REAL u_re, REAL sigma_re, TElastoPlasticData wellbore_material, int npoints, std::ostream &out, bool euler = false);

int main(int argc, char *argv[]) {
    
    int pOrder = 2; // Computational mesh order
    bool render_vtk_Q = false;
    
// Generates the geometry
    std::string source_dir = SOURCE_DIR;
    std::string msh_file = source_dir + "/Projects/WellboreGPU/gmsh/wellbore.msh";
//    std::string msh_file = source_dir + "/gmsh/wellbore-3.msh";
    TPZGeoMesh *gmesh = ReadGeometry(msh_file);
    PrintGeometry(gmesh);

// Creates the computational mesh
//    TElastoPlasticData wellbore_material = WellboreConfig(); /// NVB this one is for recurrent usage
    TElastoPlasticData wellbore_material = WellboreConfigRK(); /// NVB this one is just for verification purposes
    TPZCompMesh *cmesh = CmeshElastoplasticity(gmesh, pOrder, wellbore_material);

// Defines the analysis
    int n_threads = 0;
    TPZAnalysis *analysis = Analysis(cmesh,n_threads);
    
// Calculates the solution using Newton method
    int n_iterations = 2;
    REAL tolerance = 1.e-5;
   Solution(analysis, n_iterations, tolerance);

// Post process
   if (render_vtk_Q) {
       std::string vtk_file = "Approximation.vtk";
       PostProcess(cmesh, wellbore_material, n_threads, vtk_file);
   }

// Creates the computational mesh
    TPZCompMesh *cmesh_npts = CmeshElastoplasticity(gmesh, pOrder, wellbore_material);
    TPZAnalysis *analysis_npts = Analysis(cmesh_npts,n_threads);
    
// Calculates the solution using all intg points at once
    SolutionAllPoints(analysis_npts, n_iterations, tolerance, wellbore_material);
    std::cout << "main:: After SolutionAllPoints calling ... " << std::endl;
    if (render_vtk_Q) { //Post process
        std::string vtk_file = "Approximation_IntPointFEM.vtk";
        PostProcess(cmesh_npts, wellbore_material, n_threads, vtk_file);
    }

    return 0;
}

void Solution(TPZAnalysis *analysis, int n_iterations, REAL tolerance) {
    bool stop_criterion_Q = false;
    REAL norm_res, norm_delta_du;

    int neq = analysis->Solution().Rows();
    std::cout  << "Solving a NLS with DOF = " << neq << std::endl;

    analysis->Solution().Zero();
    TPZFMatrix<REAL> du(analysis->Solution()), delta_du;
    analysis->Assemble();

    for (int i = 0; i < n_iterations; i++) {
        analysis->Solve();
        delta_du = analysis->Solution();
        du += delta_du;
        analysis->LoadSolution(du);
        analysis->AssembleResidual();
        norm_delta_du = Norm(delta_du);
        norm_res = Norm(analysis->Rhs());
        stop_criterion_Q = norm_res < tolerance;
        std::cout << "Nonlinear process :: delta_du norm = " << norm_delta_du << std::endl;
        std::cout << "Nonlinear process :: residue norm = " << norm_res << std::endl;
        if (stop_criterion_Q) {

            AcceptPseudoTimeStepSolution(analysis, analysis->Mesh());
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i + 1 << std::endl;
            break;
        }
//        analysis->Assemble();
    }

    if (stop_criterion_Q == false) {
        AcceptPseudoTimeStepSolution(analysis, analysis->Mesh());
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

TPZAnalysis *Analysis(TPZCompMesh *cmesh, int n_threads) {
    bool optimizeBandwidth = true;
    TPZAnalysis *analysis = new TPZAnalysis(cmesh, optimizeBandwidth);
    TPZSymetricSpStructMatrix strskyl(cmesh);
    strskyl.SetNumThreads(n_threads);
    analysis->SetStructuralMatrix(strskyl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    analysis->SetSolver(step);
    return analysis;
}

void PostProcess(TPZCompMesh *cmesh, TElastoPlasticData wellbore_material, int n_threads, std::string vtk_file) {
    int div = 1;
    TPZPostProcAnalysis * post_processor = new TPZPostProcAnalysis;
    post_processor->SetCompMesh(cmesh, true);

    int n_regions = 1;
    TPZManVector<int,1> post_mat_id(n_regions);
    post_mat_id[0] = wellbore_material.Id();

    TPZStack<std::string,50> names, scalnames,vecnames,tensnames;
    scalnames.push_back("FailureType");
    vecnames.push_back("Displacement");
    tensnames.push_back("Stress");
    tensnames.push_back("Strain");
    tensnames.push_back("StrainPlastic");
    

    for (auto i : scalnames) {
        names.push_back(i);
    }
    for (auto i : vecnames) {
        names.push_back(i);
    }
    for (auto i : tensnames) {
        names.push_back(i);
    }

    post_processor->SetPostProcessVariables(post_mat_id, names);
    TPZFStructMatrix structmatrix(post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    post_processor->SetStructuralMatrix(structmatrix);

    post_processor->TransferSolution(); /// Computes the L2 projection.
    post_processor->DefineGraphMesh(2,scalnames,vecnames,tensnames,vtk_file);
    post_processor->PostProcess(div,2);
}

TElastoPlasticData WellboreConfig(){
    TPZElasticResponse LER;
    
    REAL Ey = 2000.0;
    REAL nu = 0.2;

    LER.SetEngineeringData(Ey, nu);

    REAL mc_cohesion    = 10000000000.0;
//    REAL mc_cohesion    = 5.0;
    REAL mc_phi         = (20*M_PI/180);

    /// NVB it is important to check the correct sign for ef in TPZMatElastoPlastic and TPZMatElastoPlastic2D materials. It is better to avoid problems with tensile state of stress.

    std::vector<TBCData> bc_data;
    TBCData bc_inner, bc_outer, bc_ux_fixed, bc_uy_fixed;
    bc_inner.SetId(2);
    bc_inner.SetType(6);
    bc_inner.SetInitialValue(-50.);
    bc_inner.SetValue({-1.0*(-20.-bc_inner.InitialValue())}); /// tr(sigma)/3

    bc_outer.SetId(3);
    bc_outer.SetType(6);
    bc_outer.SetInitialValue(-50.0); /// tr(sigma)/3
    bc_outer.SetValue({-1.0*(-50.-bc_outer.InitialValue())}); /// tr(sigma)/3

    bc_ux_fixed.SetId(4);
    bc_ux_fixed.SetType(3);
    bc_ux_fixed.SetValue({1,0});
    
    bc_uy_fixed.SetId(5);
    bc_uy_fixed.SetType(3);
    bc_uy_fixed.SetValue({0,1});
    
    bc_data.push_back(bc_inner);
    bc_data.push_back(bc_outer);
    bc_data.push_back(bc_ux_fixed);
    bc_data.push_back(bc_uy_fixed);
    
    TElastoPlasticData rock;
    rock.SetMaterialParameters(LER, mc_phi, mc_cohesion);
    rock.SetId(1);

    rock.SetBoundaryData(bc_data);
    
    return rock;
}

/// Material configuration for RK verification
TElastoPlasticData WellboreConfigRK(){
    
    /// Elastic verification -> true
    /// ElastoPlastic verification -> false
    bool is_elastic_Q = false;
    
    TPZElasticResponse LER;
    REAL Ey = 2000.0;
    REAL nu = 0.2;
    
    LER.SetEngineeringData(Ey, nu);
    REAL mc_cohesion, mc_phi;
    if (is_elastic_Q) {
        mc_cohesion    = 10000000000.0;
        mc_phi         = (20*M_PI/180);
    }else{
        mc_cohesion    = 5.0;
        mc_phi         = (20*M_PI/180);
    }
    
    std::vector<TBCData> bc_data;
    TBCData bc_inner, bc_outer, bc_ux_fixed, bc_uy_fixed;
    bc_inner.SetId(2);
    bc_inner.SetType(6);
    bc_inner.SetInitialValue(-50.);
    bc_inner.SetValue({-1.0*(-40.-bc_inner.InitialValue())}); /// tr(sigma)/3
    
    bc_outer.SetId(3);
    bc_outer.SetType(6);
    bc_outer.SetInitialValue(-50.0); /// tr(sigma)/3
    bc_outer.SetValue({-1.0*(-50.-bc_outer.InitialValue())+0.001}); /// tr(sigma)/3
    
    bc_ux_fixed.SetId(4);
    bc_ux_fixed.SetType(3);
    bc_ux_fixed.SetValue({1,0});
    
    bc_uy_fixed.SetId(5);
    bc_uy_fixed.SetType(3);
    bc_uy_fixed.SetValue({0,1});
    
    bc_data.push_back(bc_inner);
    bc_data.push_back(bc_outer);
    bc_data.push_back(bc_ux_fixed);
    bc_data.push_back(bc_uy_fixed);
    
    TElastoPlasticData rock;
    rock.SetMaterialParameters(LER, mc_phi, mc_cohesion);
    rock.SetId(1);
    
    rock.SetBoundaryData(bc_data);
    
// Runge Kutta approximation
    int np = 2000;
    bool euler = false;
    
    if(is_elastic_Q){

        REAL u_re, simga_re;
        u_re = 0.0000254403;
        simga_re = 0.00124267;
        std::ofstream rkfile("ElasticRKdata.txt");
        RKApproximation(u_re, simga_re, rock, np, rkfile, euler);
    }else{
        REAL u_re, simga_re;
        u_re = 0.0000318827;
        simga_re = 0.00130781;
        std::ofstream rkfile("ElastoPlasticRKdata.txt");
        RKApproximation(u_re, simga_re, rock, np, rkfile, euler);
    }
    return rock;
}

TPZGeoMesh * ReadGeometry(std::string geometry_file) {
    TPZGmshReader Geometry;
    REAL l = 1.0;
    Geometry.SetCharacteristiclength(l);
    Geometry.SetFormatVersion("4.1");
    TPZGeoMesh * geometry = Geometry.GeometricGmshMesh(geometry_file);
    Geometry.PrintPartitionSummary(std::cout);
#ifdef PZDEBUG
    if (!geometry)
    {
        std::cout << "The geometrical mesh was not generated." << std::endl;
        DebugStop();
    }
#endif
    
    return geometry;
    
}

void PrintGeometry(TPZGeoMesh * geometry) {
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name  << "geometry" << ".txt";
    vtk_name   << "geometry"  << ".vtk";
    std::ofstream textfile(text_name.str().c_str());
    geometry->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, vtkfile, true);
    
#ifdef PZDEBUG
    TPZCheckGeom checker(geometry);
    checker.CheckUniqueId();
    if(checker.PerformCheck())
    {
        DebugStop();
    }
#endif
    
}

void AcceptPseudoTimeStepSolution(TPZAnalysis * an, TPZCompMesh * cmesh){
    
    bool update = true;
    {
        std::map<int, TPZMaterial *> & refMatVec = cmesh->MaterialVec();
        std::map<int, TPZMaterial * >::iterator mit;
        TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
        for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
        {
            pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
            if(pMatWithMem != NULL)
            {
                pMatWithMem->SetUpdateMem(update);
            }
        }
    }
    an->AssembleResidual();
    update = false;
    {
        std::map<int, TPZMaterial *> & refMatVec = cmesh->MaterialVec();
        std::map<int, TPZMaterial * >::iterator mit;
        TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
        for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
        {
            pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
            if(pMatWithMem != NULL)
            {
                pMatWithMem->SetUpdateMem(update);
//                for(auto memory: *pMatWithMem->GetMemory()){
//                    memory.Print();
//                }
            }
        }
    }
    
}

void PrintMemory(TPZCompMesh * cmesh){
    
    std::map<int, TPZMaterial *> & refMatVec = cmesh->MaterialVec();
    std::map<int, TPZMaterial * >::iterator mit;
    TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem;
    for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
    {
        pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
        if(pMatWithMem != NULL)
        {
            for(auto memory: *pMatWithMem->GetMemory()){
                memory.Print();
            }
        }
    }
    
}

TPZCompMesh *CmeshElastoplasticity(TPZGeoMesh *gmesh, int p_order, TElastoPlasticData & wellbore_material) {

// Creates the computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(p_order);
    int dim = gmesh->Dimension();
    int matid = wellbore_material.Id();
    
    // Mohr Coulomb data
    REAL mc_cohesion    = wellbore_material.Cohesion();
    REAL mc_phi         = wellbore_material.FrictionAngle();
    REAL mc_psi         = mc_phi;

    // ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER = wellbore_material.ElasticResponse();

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    int PlaneStrain = 1;
    LEMC.fN.m_eps_t.Zero();
    LEMC.fN.m_eps_p.Zero();
    
    TPZElastoPlasticMem default_memory;
    default_memory.m_ER = ER;
    default_memory.m_sigma.Zero();
    default_memory.m_elastoplastic_state = LEMC.fN;

    // Creates elastoplatic material
    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >(matid,PlaneStrain);
    material->SetPlasticityModel(LEMC);
    material->SetDefaultMem(default_memory);
    cmesh->InsertMaterialObject(material);
    
    // Set the boundary conditions
    TPZFNMatrix<3,REAL> val1(dim,dim), val2(dim,dim);
    val1.Zero();
    val2.Zero();
    int n_bc = wellbore_material.BoundaryData().size();
    for (int i = 0; i < n_bc; i++) {
        int bc_id = wellbore_material.BoundaryData()[i].Id();
        int type = wellbore_material.BoundaryData()[i].Type();
        
        int n_values = wellbore_material.BoundaryData()[i].Value().size();

        for (int k = 0; k < n_values; k++) {
            val2(k,0) = wellbore_material.BoundaryData()[i].Value()[k];
        }
        
        TPZBndCondWithMem<TPZElastoPlasticMem> * bc = new  TPZBndCondWithMem<TPZElastoPlasticMem>(material, bc_id, type, val1, val2);
        cmesh->InsertMaterialObject(bc);
        
    }

    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("cmesh.txt");
    cmesh->Print(out);
#endif
    return cmesh;
}

void SolutionAllPoints(TPZAnalysis * analysis, int n_iterations, REAL tolerance, TElastoPlasticData & wellbore_material){
	std::cout << "Solving with IntPointsFEM ..." << std::endl;
    bool stop_criterion_Q = false;
    REAL norm_res, norm_delta_du;
    int neq = analysis->Solution().Rows();
    TPZFMatrix<REAL> du(neq, 1, 0.), delta_du;

    TPZIntPointsFEM solveintpoints(analysis->Mesh(), wellbore_material.Id());
    solveintpoints.SetDataStructure();

    analysis->Solution().Zero();
    analysis->Assemble();
    for (int i = 0; i < n_iterations; i++) {
        analysis->Solve();
        delta_du = analysis->Solution();
        du += delta_du;
        solveintpoints.LoadSolution(du);
        solveintpoints.AssembleResidual();
        norm_delta_du = Norm(delta_du);
        norm_res = Norm(solveintpoints.Rhs());
        stop_criterion_Q = norm_res < tolerance;
        std::cout << "Nonlinear process :: delta_du norm = " << norm_delta_du << std::endl;
        std::cout << "Nonlinear process :: residue norm = " << norm_res << std::endl;
//        PrintMemory(analysis->Mesh());
        if (stop_criterion_Q) {
            AcceptPseudoTimeStepSolution(analysis, analysis->Mesh());
//            PrintMemory(analysis->Mesh());
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i + 1 << std::endl;
            break;
        }
//        analysis->Assemble();
        analysis->Rhs() = solveintpoints.Rhs();

    }

    if (stop_criterion_Q == false) {
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

void RKApproximation (REAL u_re, REAL sigma_re, TElastoPlasticData wellbore_material, int npoints, std::ostream &out, bool euler) {
    REAL rw = 0.1;
    REAL re = 4.0;
    REAL theta = 0.;

    //Initial stress and wellbore pressure
    REAL sigma0 = wellbore_material.BoundaryData()[0].InitialValue();
    REAL pw = wellbore_material.BoundaryData()[0].Value()[0] + sigma0;

    //Outer stress
    TPZTensor<REAL> sigmaXYZ;
    REAL sigma = wellbore_material.BoundaryData()[1].Value()[0] + sigma0;
    sigmaXYZ.Identity();
    sigmaXYZ.Multiply(sigma,1);

    // Mohr Coulomb data
    REAL mc_cohesion    = wellbore_material.Cohesion();
    REAL mc_phi         = wellbore_material.FrictionAngle();
    REAL mc_psi         = mc_phi;

    // ElastoPlastic Material using Mohr Coulomb
    // Elastic predictor
    TPZElasticResponse ER = wellbore_material.ElasticResponse();

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    int PlaneStrain = 1;
    LEMC.fN.m_eps_t.Zero();
    LEMC.fN.m_eps_p.Zero();

    TPZElastoPlasticMem default_memory;
    default_memory.m_ER = ER;
    default_memory.m_sigma.Zero();
    default_memory.m_elastoplastic_state = LEMC.fN;

    // Creates elastoplatic material
    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >(wellbore_material.Id(),PlaneStrain);
    material->SetPlasticityModel(LEMC);
    material->SetDefaultMem(default_memory);

    TRKSolution rkmethod;
    rkmethod.SetWellboreRadius(rw);
    rkmethod.SetExternalRadius(re);
    rkmethod.SetNumberOfPoints(npoints);
    rkmethod.SetInitialStateMemory(default_memory);
    rkmethod.SetElastoPlasticModel(LEMC);
    rkmethod.SetRadialDisplacement(u_re);
    rkmethod.SetRadialStress(sigma_re);
    rkmethod.FillPointsMemory();
    rkmethod.RKProcess(out, euler);
}

