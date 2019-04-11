
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
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"

#include "TPZSolveMatrix.h"
#include "TElastoPlasticData.h"
#include "TRKSolution.h"

/// Gmsh mesh
TPZGeoMesh * ReadGeometry(std::string geometry_file);
void PrintGeometry(TPZGeoMesh * geometry);

/// CompMesh elastoplasticity
TPZCompMesh *CmeshElastoplasticity(TPZGeoMesh *gmesh, int pOrder, TElastoPlasticData & material_data);
TPZCompMesh *CmeshElastoplasticityNoBoundary(TPZGeoMesh *gmesh, int pOrder);

/// Material configuration
TElastoPlasticData WellboreConfig();

/// Residual calculation
void SolutionAllPoints(int n_iterations, REAL tolerance, TPZAnalysis * an, TElastoPlasticData & wellbore_material);
TPZFMatrix<REAL> AssembleResidualAllPoints(TPZFMatrix<REAL> &du, TPZCompMesh *cmesh, TElastoPlasticData & wellbore_material);
TPZFMatrix<REAL> ResidualWithoutBoundary(TPZCompMesh *cmesh, TPZCompMesh *cmesh_noboundary);

///Set Analysis
TPZAnalysis * Analysis(TPZCompMesh * cmesh, int n_threads);

///Solve using Newton method
void Solution(TPZAnalysis *analysis, int n_iterations, REAL tolerance);

/// Accept solution
void AcceptPseudoTimeStepSolution(TPZAnalysis * an, TPZCompMesh * cmesh);

/// Post process
void PostProcess(TPZCompMesh *cmesh, TElastoPlasticData material, int n_threads);

///RK Approximation
void RKApproximation (TElastoPlasticData wellbore_material, int npoints, std::ostream &outbool, bool euler = false);

int main(int argc, char *argv[]) {
    
    int pOrder = 1; // Computational mesh order
    
// Generates the geometry
    std::string file("wellbore.msh");
    TPZGeoMesh *gmesh = ReadGeometry(file);
    PrintGeometry(gmesh);

// Creates the computational mesh
    TElastoPlasticData wellbore_material = WellboreConfig();
    TPZCompMesh *cmesh = CmeshElastoplasticity(gmesh, pOrder, wellbore_material);
//    TPZCompMesh *cmesh_noboundary = CmeshElastoplasticityNoBoundary(gmesh, pOrder);

// Runge Kutta approximation
    int np = 100;
    ofstream rkfile("rkdata.txt");
    bool euler = true;
    RKApproximation(wellbore_material, np, rkfile, euler);

// Defines the analysis
    int n_threads = 0;
    TPZAnalysis *analysis = Analysis(cmesh,n_threads);

//Calculates the solution using Newton method
    int n_iterations = 10;
    REAL tolerance = 1.e-8;
    Solution(analysis, n_iterations, tolerance);

    // Post process
//    PostProcess(cmesh, wellbore_material, n_threads);

// Calculates residual without boundary conditions
//    TPZFMatrix<REAL> residual = Residual(cmesh, cmesh_noboundary);

// Calculates residual using matrix operations and check if the result is ok
    SolutionAllPoints(n_iterations, tolerance, analysis, wellbore_material);

    return 0;
}

void Solution(TPZAnalysis *analysis, int n_iterations, REAL tolerance) {
    bool stop_criterion_Q = false;
    REAL norm_res;
    int neq = analysis->Solution().Rows();
    TPZFMatrix<REAL> du(neq, 1, 0.), delta_du;

    analysis->Assemble();

    for (int i = 0; i < n_iterations; i++) {
        analysis->Solve();
        delta_du = analysis->Solution();
        du += delta_du;
        analysis->LoadSolution(du);
        analysis->AssembleResidual();
        norm_res = Norm(analysis->Rhs());
        analysis->Rhs().Print(std::cout);
        stop_criterion_Q = norm_res < tolerance;

        if (stop_criterion_Q) {
            AcceptPseudoTimeStepSolution(analysis, analysis->Mesh());
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i + 1 << std::endl;
            break;
        }
        analysis->Assemble();
    }

    if (stop_criterion_Q == false) {
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

void PostProcess(TPZCompMesh *cmesh, TElastoPlasticData wellbore_material, int n_threads) {
    int div = 1;
    TPZPostProcAnalysis * post_processor = new TPZPostProcAnalysis;
    post_processor->SetCompMesh(cmesh, true);

    int n_regions = 1;
    TPZManVector<int,1> post_mat_id(n_regions);
    post_mat_id[0] = wellbore_material.Id();

    TPZStack<std::string,50> names, scalnames,vecnames,tensnames;
    vecnames.push_back("Displacement");
    tensnames.push_back("Stress");

    for (auto i : scalnames) {
        names.push_back(i);
    }
    for (auto i : vecnames) {
        names.push_back(i);
    }
    for (auto i : tensnames) {
        names.push_back(i);
    }


    std::string vtk_file("Approximation.vtk");

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
    REAL mc_phi         = (20*M_PI/180);
    

    std::vector<TBCData> bc_data;
    TBCData bc_inner, bc_outer, bc_ux_fixed, bc_uy_fixed;
    bc_inner.SetId(2);
    bc_inner.SetType(6);
    bc_inner.SetInitialValue(-50.);
    bc_inner.SetValue({-10.-bc_inner.InitialValue()}); /// tr(sigma)/3

    bc_outer.SetId(3);
    bc_outer.SetType(6);
    bc_outer.SetInitialValue(-50.0); /// tr(sigma)/3
    bc_outer.SetValue({-50.-bc_outer.InitialValue()}); /// tr(sigma)/3

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
        TPZMaterial *bc = material->CreateBC(material, bc_id, type, val1, val2);
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

TPZCompMesh *CmeshElastoplasticityNoBoundary(TPZGeoMesh * gmesh, int p_order) {

// Creates the computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(p_order);

// Mohr Coulomb data
    REAL mc_cohesion    = 10.0;
    REAL mc_phi         = (20.0*M_PI/180);
    REAL mc_psi         = mc_phi;

// ElastoPlastic Material using Mohr Coulomb
// Elastic predictor
    TPZElasticResponse ER;
    REAL G = 400*mc_cohesion;
    REAL nu = 0.3;
    REAL E = 2.0*G*(1+nu);

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetEngineeringData(E,nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    int PlaneStrain = 1;
    int matid = 1;

// Creates elastoplatic material
    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >(matid,PlaneStrain);
    material->SetPlasticityModel(LEMC);

    cmesh->InsertMaterialObject(material);
    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();

    return cmesh;
}

void SolutionAllPoints(int n_iterations, REAL tolerance, TPZAnalysis * analysis, TElastoPlasticData & wellbore_material){
    TPZCompMesh *cmesh = analysis->Mesh();
    analysis->Assemble();
    analysis->Solve();

    bool stop_criterion_Q = false;
    REAL norm_res;
    int neq = cmesh->NEquations();
    TPZFMatrix<REAL> du(neq, 1, 0.), delta_du;
    TPZFMatrix<REAL> res;
    delta_du = analysis->Solution();

    for (int i = 0; i < n_iterations; i++) {
        du += delta_du;
        res = AssembleResidualAllPoints(du, cmesh, wellbore_material);

        norm_res = Norm(res);
        stop_criterion_Q = norm_res < tolerance;

        if (stop_criterion_Q) {
            std::cout << "Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "Number of iterations = " << i + 1 << std::endl;
            break;
        }
    }

    if (stop_criterion_Q == false) {
        std::cout << "Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

TPZFMatrix<REAL> AssembleResidualAllPoints(TPZFMatrix<REAL> &du, TPZCompMesh *cmesh, TElastoPlasticData & wellbore_material) {
    TPZFMatrix<REAL> res;
    TPZFMatrix<REAL> sigma_trial;
    TPZFMatrix<REAL> eigenvalues;
    TPZFMatrix<REAL> eigenvectors;
    TPZFMatrix<REAL> nodal_forces_vec;
    TPZFMatrix<REAL> gather_solution;
    TPZFMatrix<REAL> delta_strain;
    TPZFMatrix<REAL> delta_strainXYZ;
    TPZFMatrix<REAL> total_strain;
    TPZFMatrix<REAL> plastic_strain;
    TPZFMatrix<REAL> elastic_strain;
    TPZFMatrix<REAL> sigma_projected;
    TPZFMatrix<REAL> sigma;
    TPZFMatrix<REAL> nodal_forces;

    TPZSolveMatrix solmat(cmesh, wellbore_material);

    solmat.GatherSolution(du, gather_solution);
    solmat.DeltaStrain(gather_solution, delta_strain);
    solmat.TotalStrain(delta_strain, total_strain);
    solmat.ElasticStrain(delta_strain, total_strain, plastic_strain, elastic_strain);
    solmat.ComputeStress(elastic_strain, sigma_trial);
    solmat.SpectralDecomposition(sigma_trial, eigenvalues, eigenvectors);
    solmat.ProjectSigma(eigenvalues, sigma_projected, plastic_strain);
    solmat.StressCompleteTensor(sigma_projected, eigenvectors, sigma);
    solmat.ComputeStrain(sigma, elastic_strain);
    sigma.Print(std::cout);
    elastic_strain.Print(std::cout);
//        plastic_strain = total_strain - elastic_strain;
    solmat.NodalForces(sigma, nodal_forces);
    solmat.ColoredAssemble(nodal_forces,res);

    return res;
}

TPZFMatrix<REAL> ResidualWithoutBoundary(TPZCompMesh *cmesh, TPZCompMesh *cmesh_noboundary) {
    bool optimizeBandwidth = true;
    int n_threads = 16;

    TPZAnalysis an_d(cmesh_noboundary, optimizeBandwidth);
    TPZSymetricSpStructMatrix strskyl(cmesh_noboundary);
    strskyl.SetNumThreads(n_threads);
    an_d.SetStructuralMatrix(strskyl);

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an_d.SetSolver(step);
    an_d.Assemble();
    an_d.Solve();

    TPZFMatrix<STATE> res;
    an_d.Solver().Matrix()->Multiply(cmesh->Solution(), res);
    return res;
}

void RKApproximation (TElastoPlasticData wellbore_material, int npoints, std::ostream &out, bool euler) {
    REAL rw = 0.1;
    REAL re = 4.;
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
    rkmethod.SetExternRadius(re);
    rkmethod.SetMaterial(material);
    rkmethod.SetStressXYZ(sigmaXYZ,theta);
    rkmethod.SetInitialStress(sigma0);
    rkmethod.SetWellborePressure(pw);

    rkmethod.RKProcess(npoints,out, euler);
}

