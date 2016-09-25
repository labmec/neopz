
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <iostream>
#include <math.h>

// Geometry
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "tpzhierarquicalgrid.h"

// Computational mesh
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzreducedspace.h"

#include "pzlog.h"

// Materials
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZPoroPermCoupling.h"
#include "TPZNonLinearElliptic.h"

// Analysis
#include "pzanalysis.h"
#include "TPZGeomechanicAnalysis.h"

// Matrix
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

// Simulation data structure
#include "TPZSimulationData.h"

// Methods declarations


// Rectangular geometry
TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n);
void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);

void UniformRefinement(TPZGeoMesh *gmesh, int nh);
void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id);

static void f_xy(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
static void Sigma(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);
static void u_y(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);
static void u_xy(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);

// Create a computational mesh for nonlinear elliptic benchmark
TPZCompMesh * CMesh_Elliptic(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for nonlinear elliptic reduce base
TPZCompMeshReferred * CMesh_Elliptic_RB(TPZCompMesh * cmesh);

// Create a computational mesh for nonlinear elliptic benchmark multiphysisc version
TPZCompMesh * CMesh_Elliptic_M(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

// Create a computational mesh for nonlinear elliptic benchmark multiphysisc version
TPZCompMesh * CMesh_Elliptic_M_RB(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

// Create a computational mesh for deformation
TPZCompMesh * CMesh_Deformation(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_PorePressure(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_PorePermeabilityCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_PorePermeabilityCouplingII(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.permeabilityc"));
#endif


TPZCompMesh * OffLine_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data);

void OnLine_Benchmark(TPZCompMesh * cmesh, TPZSimulationData * sim_data);

void SetParameters(TPZCompMesh * cmesh, TPZVec<REAL> mu_vector);

int NonLinearElliptic();

int Geomchanic();

int main(int argc, char *argv[])
{
    
    NonLinearElliptic();
    
//    Geomchanic();
}

int NonLinearElliptic(){
    
    TPZMaterial::gBigNumber = 1.0e14;
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/PermeabilityCoupling/";
    FileName += "geomechanics_rom_log.cfg";
    InitializePZLOG(FileName);
#endif
    
    TPZSimulationData * sim_data = new TPZSimulationData;
    
    REAL dt = 1.0;
    int n_steps = 30;
    REAL epsilon_res = 1.0e-5;
    REAL epsilon_corr = 1.0e-10;
    int n_corrections = 50;
    
    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    g[1] = -0.0*9.81;
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    
    TPZVec<REAL> dx_dy(2);
    TPZVec<int> n(2);
    
    REAL Lx = 1.0; // meters
    REAL Ly = 1.0; // meters
    
    n[0] = 5; // x - direction
    n[1] = 5; // y - direction
    
    dx_dy[0] = Lx/REAL(n[0]); // x - direction
    dx_dy[1] = Ly/REAL(n[1]); // y - direction
    
    TPZGeoMesh * gmesh = RockBox(dx_dy,n);
    std::cout<< "Geometry done. " << std::endl;
    
    std::cout<< "off line process. " << std::endl;
    
    TPZCompMesh * nonlinear_cmesh = OffLine_Benchmark(gmesh, sim_data);
    
    std::cout<< "off line process done . " << std::endl;
    
    std::cout<< "on line process. " << std::endl;
    
    OnLine_Benchmark(nonlinear_cmesh, sim_data);
    
    std::cout<< "on line process done . " << std::endl;
    
    std::cout << " Execution finished " << std::endl;
    return EXIT_SUCCESS;
    
}

int Geomchanic(){
    
    TPZMaterial::gBigNumber = 1.0e14;
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/PermeabilityCoupling/";
    FileName += "geomechanics_rom_log.cfg";
    InitializePZLOG(FileName);
#endif
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Defining the case ... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    TPZSimulationData * sim_data = new TPZSimulationData;
    
    REAL dt = 1.0;
    int n_steps = 30;
    REAL epsilon_res = 1.0e-3;
    REAL epsilon_corr = 1.0e-10;
    int n_corrections = 50;
    
    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    g[1] = -0.0*9.81;
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Computing Geometry ... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    TPZVec<REAL> dx_dy(2);
    TPZVec<int> n(2);
    
    REAL Lx = 5.0; // meters
    REAL Ly = 10.0; // meters
    
    n[0] = 5; // x - direction
    n[1] = 10; // y - direction
    
    dx_dy[0] = Lx/REAL(n[0]); // x - direction
    dx_dy[1] = Ly/REAL(n[1]); // y - direction
    
    TPZGeoMesh * gmesh = RockBox(dx_dy,n);
    
#ifdef LOG4CXX
    
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Computing Geometry accomplished... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    // Create the approximation space
    int deformation_order = 3;
    int pore_pressure_order = 2;
    
    // Create multiphysisc mesh
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    
    mesh_vector[0] = CMesh_Deformation(gmesh, deformation_order);
    mesh_vector[1] = CMesh_PorePressure(gmesh, pore_pressure_order);
    
    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_PorePermeabilityCoupling(gmesh, mesh_vector, sim_data);
    //    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_PorePermeabilityCouplingII(gmesh, mesh_vector, sim_data);
    
    // Create the static analysis
    
    // Run Static analysis
    // @omar:: the initial condition is set up to zero for displacement and pore pressure excess
    
    // Create the Transient analysis
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 8;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(cmesh_poro_perm_coupling,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    //    TPZSkylineNSymStructMatrix skyl(cmesh_poro_perm_coupling);
    //    TPZSkylineStructMatrix struct_mat(cmesh_poro_perm_coupling);
    
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(cmesh_poro_perm_coupling);
    struct_mat.SetDecomposeType(ELDLt);
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    TPZVec<REAL> x(3);
    x[0] = Lx/2.0;
    x[1] = Ly/2.0;
    x[2] = 0.0;
    std::string file_ss_name("plot.nb");
    std::string file_sp_name("porosity.nb");
    std::string file_sk_name("permeability.nb");
    std::string file_spex_name("porepressure.nb");
    
    // Run Transient analysis
    time_analysis->Run_Evolution(x);
    time_analysis->PlotStrainStress(file_ss_name);
    time_analysis->PlotStrainPorosity(file_sp_name);
    time_analysis->PlotStrainPermeability(file_sk_name);
    time_analysis->PlotStrainPressure(file_spex_name);
    std::cout << " Execution finished" << std::endl;
    return EXIT_SUCCESS;
    
}

TPZCompMesh * OffLine_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data){
    
    std::string plotfile("Nonlinear_Elliptic.vtk");
    
    // Create the approximation space
    int potential_order = 2;
    
    // Create multiphysisc mesh
    TPZVec<TPZCompMesh * > mesh_vector(1);
    mesh_vector[0] = CMesh_Elliptic(gmesh, potential_order);
    
    TPZCompMesh * nonlinear_cmesh = CMesh_Elliptic_M(gmesh, mesh_vector, sim_data);
    
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 0;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(nonlinear_cmesh,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    //    TPZSkylineNSymStructMatrix struct_mat(nonlinear_cmesh);
//    TPZSkylineStructMatrix struct_mat(nonlinear_cmesh);
    
    TPZSymetricSpStructMatrix struct_mat(nonlinear_cmesh);
    struct_mat.SetNumThreads(number_threads);
    
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(nonlinear_cmesh);
//    struct_mat.SetDecomposeType(ELDLt);
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = nonlinear_cmesh->NEquations();
    std::cout<< "ndof = " << ndof << std::endl;
 
    // Set up the empirical interpolation
    
    TPZVec<REAL> mu_vector(3,0.0);
    int n = 2;
    mu_vector[1] = 0.0;
    mu_vector[2] = 1.0;
    
    TPZFMatrix<REAL> sol(ndof,n,0.0);
    
    for (int i = 0; i < n; i++) {
        mu_vector[0] = 2.0 + REAL(i) * 0.1;
        SetParameters(time_analysis->Mesh(), mu_vector);
        time_analysis->ExcecuteOneStep();
        time_analysis->PostNonlinearProcessStep(plotfile);
        sol.AddSub(0, i, time_analysis->X_n());
        time_analysis->X().Zero();
        time_analysis->X_n().Zero();
        time_analysis->Mesh()->Solution().Zero();
    }
    
    mesh_vector[0]->LoadSolution(sol);
    nonlinear_cmesh->LoadSolution(sol);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, nonlinear_cmesh);
//    nonlinear_cmesh->Solution().Print("sol = ");
    return nonlinear_cmesh;
    
}

void SetParameters(TPZCompMesh * cmesh, TPZVec<REAL> mu_vector){
    int mat_id = 1;
    TPZMaterial * material = cmesh->MaterialVec().find(mat_id)->second;
#ifdef PZDEBUG
    if(!material)
    {
        DebugStop();
    }
#endif
    
    TPZNonLinearElliptic * elliptic_material = dynamic_cast<TPZNonLinearElliptic * >(material);
    
#ifdef PZDEBUG
    if(!elliptic_material)
    {
        DebugStop();
    }
#endif
    REAL mu_0 = mu_vector[0];
    REAL mu_1 = mu_vector[1];
    REAL mu_2 = mu_vector[2];
    elliptic_material->SetParameters(mu_0, mu_1, mu_2);
    
    
}

void OnLine_Benchmark(TPZCompMesh * cmesh, TPZSimulationData * sim_data){
    
    // Create the approximation space
    TPZGeoMesh * gmesh = cmesh->Reference();
    
    std::string plotfile("Nonlinear_Elliptic_bc.vtk");
    
    // Create multiphysisc mesh
    TPZVec<TPZCompMesh * > mesh_vector(1);
    mesh_vector[0] = CMesh_Elliptic_RB(cmesh);

    TPZCompMesh * nonlinear_rb_cmesh = CMesh_Elliptic_M_RB(gmesh, mesh_vector, sim_data);
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 0;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(nonlinear_rb_cmesh,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    //    TPZSkylineNSymStructMatrix struct_mat(nonlinear_rb_cmesh);
    TPZSkylineStructMatrix struct_mat(nonlinear_rb_cmesh);
    
    //    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(nonlinear_rb_cmesh);
    //    struct_mat.SetDecomposeType(ELDLt);
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    std::cout<< "ndof = " << nonlinear_rb_cmesh->NEquations() << std::endl;
    
    time_analysis->ExcecuteOneStep();
    time_analysis->PostNonlinearProcessStep(plotfile);
    
    time_analysis->ExcecuteOneStep();
    time_analysis->PostNonlinearProcessStep(plotfile);
}

void UniformRefinement(TPZGeoMesh *gmesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (sons);
        }//for i
    }//ref
}

void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1){
                if (gel->MaterialId()== mat_id){
                    gel->Divide (sons);
                }
            }
        }//for i
    }//ref
}

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_Elliptic(TPZGeoMesh * gmesh, int order){

    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZNonLinearElliptic * material = new TPZNonLinearElliptic(matid,dim);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();

#ifdef PZDEBUG
    std::ofstream out("CmeshElliptic.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

// Create a computational mesh for nonlinear elliptic reduce base
TPZCompMeshReferred * CMesh_Elliptic_RB(TPZCompMesh * cmesh){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh * gmesh = cmesh->Reference();
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMeshReferred * cmesh_rb = new TPZCompMeshReferred(gmesh);
    
    
    // Creating a material object
    TPZNonLinearElliptic * material = new TPZNonLinearElliptic(matid,dim);
    cmesh_rb->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh_rb->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh_rb->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh_rb->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh_rb->InsertMaterialObject(bc_left_mat);
    
    // Setting RB approximation space
    cmesh_rb->SetDimModel(dim);
    int numsol = cmesh->Solution().Cols();
    cmesh_rb->AllocateNewConnect(numsol, 1, 1);
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh_rb);
    cmesh_rb->AutoBuild();
    
    cmesh_rb->AdjustBoundaryElements();
    cmesh_rb->CleanUpUnconnectedNodes();
    cmesh_rb->LoadReferred(cmesh);
#ifdef PZDEBUG
    std::ofstream out("CmeshEllipticRB.txt");
    cmesh_rb->Print(out);
#endif
    
    return cmesh_rb;
    
}

// Create a computational mesh for nonlinear elliptic benchmark multiphysisc version
TPZCompMesh * CMesh_Elliptic_M(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    REAL mu_0 = 2.0;
    REAL mu_1 = 0.0;
    REAL mu_2 = 20.0;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZNonLinearElliptic * material = new TPZNonLinearElliptic(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetParameters(mu_0, mu_1, mu_2);
    
    TPZFunction<REAL> * fun_xy = new TPZDummyFunction<REAL>(f_xy);
    material->SetTimeDependentForcingFunction(fun_xy);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet_u    = 0;
    int neumann_sigma  = 1;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet_u, val1, val2);
//    TPZFunction<REAL> * boundary_data = new TPZDummyFunction<REAL>(Sigma);
//    bc_top_mat->SetTimedependentBCForcingFunction(boundary_data);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    
    long nel = cmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CMeshEllipticMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

// Create a computational mesh for nonlinear elliptic benchmark multiphysisc version
TPZCompMesh * CMesh_Elliptic_M_RB(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    REAL mu_0 = 2.0;
    REAL mu_1 = 0.0;
    REAL mu_2 = 20.0;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);

    // Creating a material object
    TPZNonLinearElliptic * material = new TPZNonLinearElliptic(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetParameters(mu_0, mu_1, mu_2);
    
    TPZFunction<REAL> * fun_xy = new TPZDummyFunction<REAL>(f_xy);
    material->SetTimeDependentForcingFunction(fun_xy);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet_u    = 0;
    int neumann_sigma  = 1;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet_u, val1, val2);
    //    TPZFunction<REAL> * boundary_data = new TPZDummyFunction<REAL>(Sigma);
    //    bc_top_mat->SetTimedependentBCForcingFunction(boundary_data);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_u, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
  
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    long nel = cmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();       
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CMeshEllipticMultiPhysicsRB.txt");
    cmesh->Print(out);
#endif

    return cmesh;
    
}



TPZCompMesh * CMesh_PorePermeabilityCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    // Plane strain assumption
    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    REAL MPa = 1.0e6;
    REAL rad = M_PI/180.0;
    
    // Getting mesh dimension
    int dim = 2;
    
    int kmodel = 3;
    REAL l = 15.3333e8;
    REAL mu = 5.1111e8;
    REAL l_u = 16.3333e8;
    REAL alpha = 0.95;
    REAL Se = 1.0e-7;
    REAL k = 1.0e-14;
    REAL porosity = 0.25;
    REAL eta = 0.001;
    
    REAL c = 27.2*MPa;
    REAL phi_f = 30.0*rad;

    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZPoroPermCoupling * material = new TPZPoroPermCoupling(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPlaneProblem(planestress);
    material->SetPorolasticParameters(l, mu, l_u);
    material->SetBiotParameters(alpha, Se);
    material->SetParameters(k, porosity, eta);
    material->SetKModel(kmodel);
    material->SetDruckerPragerParameters(phi_f, c);
    cmesh->InsertMaterialObject(material);
    
    
    // Inserting boundary conditions
    int dirichlet_x_vn   = 7;
    int dirichlet_y_vn   = 8;
    int neumann_y_p      = 5;
    int dirichlet_y_p    = 2;

    REAL s_n = -10.0*MPa;
//    REAL u_y = -0.000333333;
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0013801;
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_y_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet_x_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = s_n;
    val2(2,0) = 0.0;
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, neumann_y_p, val1, val2);
    TPZFunction<REAL> * boundary_data = new TPZDummyFunction<REAL>(Sigma);
    bc_top_mat->SetTimedependentBCForcingFunction(boundary_data);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_x_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
//    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    
    long nel = cmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CMeshMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * CMesh_PorePermeabilityCouplingII(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    // Plane strain assumption
    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    REAL MPa = 1.0e6;
    REAL rad = M_PI/180.0;
    
    // Getting mesh dimension
    int dim = 2;
    int kmodel = 3;
    REAL l = 15.3333e8;
    REAL mu = 5.1111e8;
    REAL l_u = 16.3333e8;
    REAL alpha = 0.25;
    REAL Se = 1.0e-8;
    REAL k = 1.0e-13;
    REAL porosity = 0.25;
    REAL eta = 0.001;
    
    REAL c = 27.2*MPa;
    REAL phi_f = 10.0*rad;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZPoroPermCoupling * material = new TPZPoroPermCoupling(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPlaneProblem(planestress);
    material->SetPorolasticParameters(l, mu, l_u);
    material->SetBiotParameters(alpha, Se);
    material->SetParameters(k, porosity, eta);
    material->SetKModel(kmodel);
    material->SetDruckerPragerParameters(phi_f, c);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int neumann_xy_vn    = 9;
    int dirichlet_xy_vn  = 6;
    int neumann_y_p      = 5;
    int dirichlet_xy_p   = 0;
    int dirichlet_y_vn   = 8;
    
    REAL s_n = -10.0*MPa;
    //    REAL u_y = -0.000333333;
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = -100.0013801;
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_y_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, neumann_xy_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = s_n;
    val2(2,0) = 0.0;
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet_xy_p, val1, val2);
    TPZFunction<REAL> * boundary_data = new TPZDummyFunction<REAL>(u_xy);
    bc_top_mat->SetTimedependentBCForcingFunction(boundary_data);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, neumann_xy_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    //    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(mesh_vector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, cmesh);
    
    
    long nel = cmesh->NElements();
    TPZVec<long> indices;
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->InitializeIntegrationRule();
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CMeshMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

void f_xy(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL x = pt[0];
    REAL y = pt[1];
    
    REAL f_val = 100.0*sin(2.0*M_PI*x)*sin(2.0*M_PI*y);
    
    f[0] = f_val;
    return;
}

void Sigma(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP)
{
    
    REAL MPa = 1.0e6;
    REAL s_n = (1.0*time+0.25)*MPa;
    
    f[0] = 0.0;
    f[1] = -s_n;
    f[2] = 0.0;
    return;
}

void u_y(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP)
{
    
    REAL uy = -(0.1/100.0)*time*0.2;
    
    f[0] = 0.0;
    f[1] = uy;
    f[2] = 0.0;
    return;
}

void u_xy(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP)
{
    
    REAL u = -(1.0/100.0)*time*0.2;
    
    f[0] = u;
    f[1] = u;
    f[2] = 0.0;
    return;
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_PorePressure(TPZGeoMesh * gmesh, int order){
    
    // Plane strain assumption
    //    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    // Creating a material object
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshPorePressure.txt");
    cmesh->Print(out);
#endif

    return cmesh;
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_Deformation(TPZGeoMesh * gmesh, int order){
        
    // Plane strain assumption
//    int planestress = 0;
    
    // Material identifiers
    int matid =1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    // Getting mesh dimension
    int dim = 2;
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);


    // Creating a material object
    int nstate = 2;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet = 0;

    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    TPZMaterial * bc_right_mat = material->CreateBC(material, bc_right, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    TPZMaterial * bc_top_mat = material->CreateBC(material, bc_top, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshDeformation.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}



TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n){
    
    REAL t=0.0,dx,dy;
    int n_elements;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh1 = new TPZGeoMesh;
    GeoMesh1->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh1->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    int bc_bottom, bc_right, bc_top, bc_left;
    bc_bottom = -1;
    bc_right = -2;
    bc_top = -3;
    bc_left = -4;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh1);
    GeoMesh1->BuildConnectivity();
    GeoMesh1->SetDimension(0);
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew1.txt");
//        GeoMesh1->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew1.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh1,Dummyfile, true);
//    }
    
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh1);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom.SetParametricFunction(ParFunc);
    CreateGridFrom.SetFrontBackMatId(bc_left,bc_right);
    dx = dx_dy[0];
    n_elements = n[0];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh2 = CreateGridFrom.ComputeExtrusion(t, dx, n_elements);
    
//    {
//        //  Print Geometrical Base Mesh
//        std::ofstream argument("GeometicMeshNew2.txt");
//        GeoMesh2->Print(argument);
//        std::ofstream Dummyfile("GeometricMeshNew2.vtk");
//        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh2,Dummyfile, true);
//    }
    
    
    
    TPZHierarquicalGrid CreateGridFrom2(GeoMesh2);
    TPZAutoPointer<TPZFunction<STATE> > ParFunc2 = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(bc_bottom,bc_top);
    dy = dx_dy[1];
    n_elements = n[1];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dy, n_elements);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("Geometry.txt");
        GeoMesh3->Print(argument);
        std::ofstream Dummyfile("Geometry.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh3,Dummyfile, true);
    }
    
    long last_node = GeoMesh3->NNodes() - 1;
    long last_element = GeoMesh3->NElements() - 1;
    long node_id = GeoMesh3->NodeVec()[last_node].Id();
    long element_id = GeoMesh3->Element(last_element)->Id();
    const std::string name("Geomechanic Reservoir box");
    GeoMesh3->SetName(name);
    GeoMesh3->SetMaxNodeId(node_id);
    GeoMesh3->SetMaxElementId(element_id);
    GeoMesh3->SetDimension(2);
    return GeoMesh3;
    
}

void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0]; // x
    X[1] = 0.0; // y
    X[2] = 0.0; // z
}

void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0; // x
    X[1] = par[0]; // y
    X[2] = 0.0; // z
}

