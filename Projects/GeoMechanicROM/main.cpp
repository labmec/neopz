
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
#include "TPZGmshReader.h"

// Computational mesh
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzreducedspace.h"

#include "pzlog.h"

// Materials
#include "pzelast3d.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZPoroPermCoupling.h"
#include "TPZNonLinearElliptic.h"
#include "TPZLinearElliptic.h"
#include "TPZBiotPoroelasticity.h"
#include "TPZPoroelasticModes.h"
#include "TPZElasticBiot.h"
#include "TPZDarcyFlow.h"

// Analysis
#include "pzanalysis.h"
#include "TPZGeomechanicAnalysis.h"
#include "TPZElasticAnalysis.h"
#include "TPZFLuxPressureAnalysis.h"
#include "TPZSegregatedSolver.h"

// Matrix
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZBiotIrregularBlockDiagonal.h"

// Simulation data structure
#include "TPZSimulationData.h"

#include "tpzautopointer.h"
#include "pzfunction.h"

// Transfer object
#include "TPZTransferFunctions.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

// Solutions
#define Solution1

static void Analytic(const TPZVec<REAL> &x, REAL time, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu);

// Methods declarations


// Rectangular geometry
TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n, bool IsTriangleMeshQ);
void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);

/** @brief Create a reservoir-box geometry with cylindrical wells */
TPZGeoMesh * CreateGeometricGmshMesh(std::string &grid);

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

// Create a computational mesh for flux
TPZCompMesh * CMesh_Flux(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for mixed pressure
TPZCompMesh * CMesh_MFPorePressure(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for pore pressure excess
TPZCompMesh * CMesh_GeomechanicCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data, bool IsMixedQ);

#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.permeabilityc"));
#endif

// Nonlinear Methods

TPZCompMesh * Reference_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data);

TPZCompMesh * OffLine_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data, TPZVec<TPZCompMesh * > & mesh_vector, int n, TPZFMatrix<REAL> &G_projts);

TPZCompMesh * OnLine_Benchmark(TPZCompMesh * cmesh, TPZSimulationData * sim_data);

REAL Error_Benchmark(TPZCompMesh * cmesh_rb, TPZCompMesh * cmesh_ref);


// Geomechanic Methods H1-H1 case


// Create a computational mesh for define unit pressures
TPZCompMesh * CMesh_Pressures(TPZGeoMesh * gmesh);

// Create a computational mesh for integrate the projections of the unit pressures
TPZCompMesh * CMesh_Elasticity(TPZGeoMesh * gmesh, int order);

// Create a computational mesh for basis generation multiphysisc version
TPZCompMesh * CMesh_GeoModes_M(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

// New segregated scheme
TPZCompMesh * CMesh_Elliptic(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

TPZCompMesh * CMesh_Parabolic(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data, bool IsMixedQ);

// Compute Galerkin projections of unit pressures
TPZCompMesh * Galerkin_Projections(TPZGeoMesh * gmesh, TPZSimulationData * sim_data, int order, int level);

int DrawUnitPressuresBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures, int level);

int DrawingPressureBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures, TPZSimulationData * sim_data);

bool DrawingGeometryOutline(TPZGeoMesh * gmesh, TPZStack< REAL > & min_x, TPZStack< REAL > & max_x);

void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes);

// Create a computational mesh for reduced deformation
TPZCompMesh * CMesh_Deformation_rb(TPZCompMesh * cmesh, int order);


void SetParameters(TPZCompMesh * cmesh, TPZVec<REAL> mu_vector);

int NonLinearElliptic();

int Geomechanic();

int Segregated_Geomechanic();

// Material identifiers
const int matid = 1;
const int bc_bottom = 2;
const int bc_right = 3;
const int bc_left = 3;
const int bc_top = 4;
const int bc_top_null = 5;

int main(int argc, char *argv[])
{
    
//    NonLinearElliptic();
    

    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Running whole process
    Geomechanic();
    
//    Segregated_Geomechanic();
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    std::cout  << " Overal execution time = " << (t2-t1) << std::endl;
#endif
}

int NonLinearElliptic(){
    
    TPZMaterial::gBigNumber = 1.0e12;
    
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
    REAL epsilon_res = 1.0e-2;
    REAL epsilon_corr = 1.0e-6;
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
    
    n[0] = 1; // x - direction
    n[1] = 1; // y - direction
    
    dx_dy[0] = Lx/REAL(n[0]); // x - direction
    dx_dy[1] = Ly/REAL(n[1]); // y - direction

    TPZGeoMesh * gmesh = RockBox(dx_dy,n,false);
    UniformRefinement(gmesh, 4);
    std::cout<< "Geometry done. " << std::endl;

    int n_data = 3; // @omar:: ??
    TPZFMatrix<REAL> projections, sel_projections;
    TPZVec<TPZCompMesh * > mesh_vector;
    std::cout<< "off line process. " << std::endl;
    TPZCompMesh * nonlinear_cmesh = OffLine_Benchmark(gmesh, sim_data, mesh_vector, n_data, projections);
    std::cout<< "off line process done . " << std::endl;
    
    std::cout<< "Reference solution. " << std::endl;
    TPZCompMesh * nonlinear_ref_cmesh = Reference_Benchmark(gmesh, sim_data);
    std::cout<< "Reference solution. " << std::endl;
    
    int n_rb = 3;
    TPZFMatrix<REAL> errors(n_rb,2,0.0);
    for (int irb = 1; irb <= n_rb; irb ++) {
        projections.GetSub(0, 0, projections.Rows(), irb*irb, sel_projections);
        mesh_vector[0]->LoadSolution(sel_projections);
        nonlinear_cmesh->LoadSolution(sel_projections);
        
        std::cout<< "on line process. " << std::endl;
        TPZCompMesh * nonlinear_rb_cmesh = OnLine_Benchmark(nonlinear_cmesh, sim_data);
        std::cout<< "on line process done . " << std::endl;
        
        REAL error = Error_Benchmark(nonlinear_rb_cmesh,nonlinear_ref_cmesh);
        
        errors(irb-1,0) = REAL(irb-1);
        errors(irb-1,1) = sqrt(error);
    }

    errors.Print("data = ",std::cout,EMathematicaInput);
    
    std::cout << " Execution finished " << std::endl;
    return EXIT_SUCCESS;
    
}

int Geomechanic(){
    
    HDivPiola = 1;
    TPZMaterial::gBigNumber = 1.0e14;
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/PermeabilityCoupling/";
    FileName += "geomechanics_rom_log.cfg";
    InitializePZLOG(FileName);
#endif
    
    TPZSimulationData * sim_data = new TPZSimulationData;
    
    REAL dt = 1.0e-5;
    int n_steps = 1e6;
    REAL epsilon_res = 1.0e-8;
    REAL epsilon_corr = 1.0e-8;
    int n_corrections = 10;
    bool IsMixedQ = false;
    bool IsRBQ    = false;
    

    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    sim_data->SetRBApproxQ(IsRBQ);
    
    std::string dirname = PZSOURCEDIR;
    std::string file;
    file = dirname + "/Projects/GeoMechanicROM/mesh/Column_Problem.msh";
//    file = dirname + "/Projects/GeoMechanicROM/mesh/Footing_Problem.msh";
    TPZGeoMesh * gmesh = CreateGeometricGmshMesh(file);

    int order = 2;
    int level = 0; // deprecated
    int hlevel = 3;
    
    UniformRefinement(gmesh, hlevel);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("Geometry.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("Geometry.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    std::cout<< "Geometry done. " << std::endl;

    TPZTransferFunctions * transfer = new TPZTransferFunctions;
    TPZCompMesh * cmesh_gp = new TPZCompMesh;
    if (IsRBQ) {
        cmesh_gp = Galerkin_Projections(gmesh, sim_data, order,level);
        transfer->SetCmeshGalerkingProjections(cmesh_gp);
        transfer->SetSimulationData(sim_data);
    }
    
    int n_meshes = 2;
    if (IsMixedQ) {
        n_meshes = 3;
    }

    TPZVec<TPZCompMesh * > mesh_vector(n_meshes);

    if (IsRBQ) {
        mesh_vector[0] = CMesh_Deformation_rb(cmesh_gp, order); // RB mesh
    }
    else{
        mesh_vector[0] = CMesh_Deformation(gmesh, order); // Full order mesh
    }
    
    if (IsMixedQ) {
        mesh_vector[1] = CMesh_Flux(gmesh, order-1);
        mesh_vector[2] = CMesh_MFPorePressure(gmesh, order-1);
    }
    else{
        mesh_vector[1] = CMesh_PorePressure(gmesh, order-1);
    }

    TPZCompMesh * geomechanic = CMesh_GeomechanicCoupling(gmesh, mesh_vector, sim_data,IsMixedQ);
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 16;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(geomechanic,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    if (IsRBQ) {

#ifdef USING_BOOST
        boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
        transfer->FillGeomechanicElPairs(geomechanic);
        transfer->RB_basis_To_Geomechanic_Memory(geomechanic);
        time_analysis->SetTransfer_object(transfer);
        
#ifdef USING_BOOST
        boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
        
#ifdef USING_BOOST
        std::cout  << "RB:: construction of transfer object execution time = " << (t2-t1) << std::endl;
#endif
    }
    
//    TPZSkylineNSymStructMatrix struct_mat(geomechanic);
//    TPZSkylineStructMatrix struct_mat(geomechanic);

    TPZSymetricSpStructMatrix struct_mat(geomechanic);
    struct_mat.SetNumThreads(number_threads);
    
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(geomechanic);
//    struct_mat.SetDecomposeType(ELDLt);

    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = geomechanic->NEquations();
    
    if (IsRBQ) {
        std::cout << " RB order model ndof =  " << ndof << std::endl;
    }
    else{
        std::cout << " Full order model ndof =  " << ndof << std::endl;
    }
    
    TPZVec<REAL> x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    
    
    std::string plotfile("geomechanic_0.vtk");

    if (IsRBQ) {
        plotfile = "geomechanic_rb_0.vtk";
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Running whole process
    // Run Transient analysis
    time_analysis->Run_Evolution(x,plotfile);
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    if (IsRBQ) {
        std::cout  << "RB:: online stage execution time = " << (t2-t1) << std::endl;
    }
    else{
        std::cout  << "Full order:: execution time = " << (t2-t1) << std::endl;
    }
#endif

    std::cout << " Full coupled DoF = " << time_analysis->Solution().Rows() << std::endl;
    
    std::cout << " Execution finished " << std::endl;
    return EXIT_SUCCESS;
    
}

int Segregated_Geomechanic(){
    
    HDivPiola = 1;
    TPZMaterial::gBigNumber = 1.0e14;
    
#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + "/Projects/PermeabilityCoupling/";
    FileName += "geomechanics_rom_log.cfg";
    InitializePZLOG(FileName);
#endif
    
    TPZSimulationData * sim_data = new TPZSimulationData;
    
    REAL dt = 0.01;
    int n_steps = 1000;
    REAL epsilon_res  = 1.0e-3;
    REAL epsilon_corr = 1.0e-3;
    int n_corrections = 150;
    bool IsMixedQ = false;
    bool IsRBQ    = false;
    
    int order = 2;
    int hlevel = 0;
    
    TPZStack< int > blocks;
    blocks.Push(20);
    blocks.Push(20);
    blocks.Push(1);
    
    std::string elliptic_file   = "elliptic.vtk";
    std::string parabolic_file  = "parabolic.vtk";
    
    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    sim_data->SetMixedApproxQ(IsMixedQ);
    sim_data->SetRBApproxQ(IsRBQ);
    sim_data->SetBlocks(blocks);
    
    std::string dirname = PZSOURCEDIR;
    std::string file;
    file = dirname + "/Projects/GeoMechanicROM/mesh/Column_Problem.msh";
//    file = dirname + "/Projects/GeoMechanicROM/mesh/Footing_Problem.msh";
    TPZGeoMesh * gmesh = CreateGeometricGmshMesh(file);
    
    UniformRefinement(gmesh, hlevel);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("Geometry.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("Geometry.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    std::cout<< "Geometry done. " << std::endl;
    
    
    TPZTransferFunctions * transfer = new TPZTransferFunctions;
    
    TPZCompMesh * cmesh_gp = new TPZCompMesh;
    if (IsRBQ) {
        cmesh_gp = Galerkin_Projections(gmesh, sim_data, order,0);
        transfer->SetCmeshGalerkingProjections(cmesh_gp);
    }
    
    int n_meshes = 1;
    if (IsMixedQ) {
        n_meshes = 2;
    }
    
    TPZVec<TPZCompMesh * > elliptic_ini_mesh_vec(1);
    TPZVec<TPZCompMesh * > elliptic_mesh_vec(1);
    TPZVec<TPZCompMesh * > parabolic_mesh_vec(n_meshes);
    
    //    Space Creation
    if (IsRBQ) {
        elliptic_mesh_vec[0] = CMesh_Deformation_rb(cmesh_gp, order); // RB mesh
    }
    else{
        elliptic_mesh_vec[0]  = CMesh_Deformation(gmesh, order); // Full order mesh
    }
    
    elliptic_ini_mesh_vec[0]  = CMesh_Deformation(gmesh, order); // Full order mesh

    if (IsMixedQ) {
        order--;        
        parabolic_mesh_vec[0] = CMesh_Flux(gmesh, order);
        parabolic_mesh_vec[1] = CMesh_MFPorePressure(gmesh, order);
    }
    else{
        order--;
        parabolic_mesh_vec[0] = CMesh_PorePressure(gmesh, order);
    }
    
    // Filling the transfer object
    transfer->SetSimulationData(sim_data);
    int number_threads = 16;
    
    // Elliptic problem
    TPZCompMesh * cmesh_elliptic_ini = CMesh_Elliptic(gmesh, elliptic_ini_mesh_vec, sim_data);
    TPZCompMesh * cmesh_elliptic     = CMesh_Elliptic(gmesh, elliptic_mesh_vec, sim_data);
    
    bool OptimizeBand_e = true;
    TPZElasticAnalysis * elliptic_ini = new TPZElasticAnalysis;
    elliptic_ini->SetCompMesh(cmesh_elliptic_ini,OptimizeBand_e);
    elliptic_ini->SetSimulationData(sim_data);
    elliptic_ini->SetMeshvec(elliptic_ini_mesh_vec);
    elliptic_ini->AdjustVectors();
    elliptic_ini->SetTransfer_object(transfer);
    
    TPZSkylineStructMatrix struct_mat_e_ini(cmesh_elliptic_ini);
    TPZStepSolver<STATE> step_e_ini;
    struct_mat_e_ini.SetNumThreads(number_threads);
    step_e_ini.SetDirect(ECholesky);
    elliptic_ini->SetSolver(step_e_ini);
    elliptic_ini->SetStructuralMatrix(struct_mat_e_ini);
    
    TPZElasticAnalysis * elliptic = new TPZElasticAnalysis;
    elliptic->SetCompMesh(cmesh_elliptic,OptimizeBand_e);
    elliptic->SetSimulationData(sim_data);
    elliptic->SetMeshvec(elliptic_mesh_vec);
    elliptic->AdjustVectors();
    elliptic->SetTransfer_object(transfer);
    
    TPZSkylineStructMatrix struct_mat_e(cmesh_elliptic);
    TPZStepSolver<STATE> step_e;
    struct_mat_e.SetNumThreads(number_threads);
    step_e.SetDirect(ELDLt);
    elliptic->SetSolver(step_e);
    elliptic->SetStructuralMatrix(struct_mat_e);
    
    
    // Parabolic problem
    TPZCompMesh * cmesh_parabolic = CMesh_Parabolic(gmesh, parabolic_mesh_vec, sim_data, IsMixedQ);
    TPZFLuxPressureAnalysis * parabolic = new TPZFLuxPressureAnalysis;
    bool OptimizeBand_p = true;
    parabolic->SetCompMesh(cmesh_parabolic,OptimizeBand_p);
    parabolic->SetSimulationData(sim_data);
    parabolic->SetMeshvec(parabolic_mesh_vec);
    parabolic->AdjustVectors();
    parabolic->SetTransfer_object(transfer);
    
    TPZSkylineStructMatrix struct_mat_p(cmesh_parabolic);
    TPZStepSolver<STATE> step_p;
    struct_mat_p.SetNumThreads(number_threads);
    step_p.SetDirect(ELDLt);
    parabolic->SetSolver(step_p);
    parabolic->SetStructuralMatrix(struct_mat_p);
    
    std::cout << " elliptic DoF = " << elliptic->Solution().Rows() << std::endl;
    std::cout << " parabolic DoF = " << parabolic->Solution().Rows() << std::endl;
    
    // Transfer object
    
    if (IsRBQ) {
        
        // Build linear tranformations initial
        transfer->Fill_elliptic_To_elliptic(cmesh_elliptic_ini);
        transfer->Fill_elliptic_To_parabolic(cmesh_elliptic_ini, cmesh_parabolic);
        
        // Build linear tranformations
        transfer->Fill_gp_elliptic_To_rb_elliptic(cmesh_gp, cmesh_elliptic);
        transfer->Fill_gp_elliptic_To_parabolic(cmesh_gp, cmesh_parabolic);
        transfer->rb_elliptic_To_rb_elliptic(cmesh_elliptic);
        transfer->rb_elliptic_To_parabolic(cmesh_elliptic, cmesh_parabolic);
        
        transfer->Fill_parabolic_To_parabolic(cmesh_parabolic);
        transfer->Fill_parabolic_To_elliptic(cmesh_parabolic, cmesh_elliptic);
        
        transfer->parabolic_To_elliptic(cmesh_parabolic, cmesh_elliptic);
        transfer->parabolic_To_parabolic(cmesh_parabolic);
        
        // transfer approximation space to integration points
        transfer->space_To_elliptic(cmesh_elliptic_ini);        
        transfer->rb_space_To_rb_elliptic(cmesh_elliptic);
        transfer->space_To_parabolic(cmesh_parabolic);
    }
    else{

        // Build linear tranformations initial
        transfer->Fill_elliptic_To_elliptic(cmesh_elliptic_ini);
        transfer->Fill_elliptic_To_parabolic(cmesh_elliptic_ini, cmesh_parabolic);
        
        // Build linear tranformations
        transfer->Fill_elliptic_To_elliptic(cmesh_elliptic);
        transfer->Fill_elliptic_To_parabolic(cmesh_elliptic, cmesh_parabolic);
        transfer->Fill_parabolic_To_parabolic(cmesh_parabolic);
        transfer->Fill_parabolic_To_elliptic(cmesh_parabolic, cmesh_elliptic);
        
        // transfer approximation space to integration points
        transfer->space_To_elliptic(cmesh_elliptic_ini);
        transfer->space_To_elliptic(cmesh_elliptic);
        transfer->space_To_parabolic(cmesh_parabolic);
    }
    
    // Run segregated solution
    TPZSegregatedSolver * segregated = new TPZSegregatedSolver;
    segregated->Set_elliptic_ini(elliptic_ini);
    segregated->Set_elliptic(elliptic);
    segregated->Set_parabolic(parabolic);
    segregated->SetTransfer_object(transfer);
    segregated->SetSimulationData(sim_data);
    
//    if(IsMixedQ){
//        elliptic_file = "elliptic_mf.vtk";
//        parabolic_file = "parabolic_mf.vtk";
//        
//        if (IsRBQ) {
//            elliptic_file = "elliptic_mf_rb.vtk";
//            parabolic_file = "parabolic_mf_rb.vtk";
//        }
//    }
//    else{
//        if (IsRBQ) {
//            elliptic_file = "elliptic_rb.vtk";
//            parabolic_file = "parabolic_rb.vtk";
//        }
//    }
    
    if (IsRBQ) {
        elliptic_file = "elliptic_rb_" + std::to_string(elliptic->Solution().Rows()) + ".vtk";
        parabolic_file = "parabolic_rb_" + std::to_string(elliptic->Solution().Rows()) + ".vtk";
    }

    
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    segregated->Run_Evolution(elliptic_file, parabolic_file);
    
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    if (IsRBQ) {
        std::cout  << "RB Fixed-stress split :: online stage execution time = " << (t2-t1) << std::endl;
    }
    else{
        std::cout  << "Fixed-stress split with full order:: execution time = " << (t2-t1) << std::endl;
    }
#endif
    
    std::cout << " elliptic DoF = " << elliptic->Solution().Rows() << std::endl;
    std::cout << " parabolic DoF = " << parabolic->Solution().Rows() << std::endl;
    
    return 0;
    
}


TPZCompMesh * Galerkin_Projections(TPZGeoMesh * gmesh, TPZSimulationData * sim_data, int order, int level){
    
    TPZVec<TPZCompMesh * > mesh_vector(2);
    mesh_vector[0] = CMesh_Elasticity(gmesh, order);
    mesh_vector[1] = CMesh_Pressures(gmesh);
    
    TPZCompMesh * geo_modes = CMesh_GeoModes_M(gmesh, mesh_vector, sim_data);
    
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 16;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(geo_modes,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    TPZSymetricSpStructMatrix struct_mat(geo_modes);
    struct_mat.SetNumThreads(number_threads);
    
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(geo_modes);
//    struct_mat.SetDecomposeType(ELDLt);
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = geo_modes->NEquations();
    
    std::cout<< "RB:: offline strategy ndof = " << ndof << std::endl;
    
    time_analysis->SimulationData()->SetCurrentStateQ(true);
    time_analysis->Assemble();
    
    // Setting up the empirical interpolation based on unitary pressures
    std::string plotfile("Geo_Modes_rb_0.vtk");
    
    REAL unit_p = 1.0e2;
    TPZStack<TPZVec<long> > cts_pressures;
    
//    int n_blocks = DrawUnitPressuresBlocks(mesh_vector[1],cts_pressures,level);
    int n_blocks = DrawingPressureBlocks(mesh_vector[1], cts_pressures, sim_data);
    
    int ndof_elastic = mesh_vector[0]->NEquations();
    TPZFMatrix<REAL> galerkin_projts(ndof_elastic,n_blocks);
    galerkin_projts.Zero();
    
    int progress = (n_blocks/10) + 1;
    REAL percent = -10.0;
    std::cout<< "RB:: number of geomodes = " << n_blocks << std::endl;
    for (int ip = 0; ip < n_blocks; ip++) {

        
        mesh_vector[0]->Solution().Zero();
        mesh_vector[1]->Solution().Zero();
        
        for (int jp = 0; jp < cts_pressures[ip].size(); jp++) {
            mesh_vector[1]->Solution()(cts_pressures[ip][jp],0) = unit_p;
        }

        TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, time_analysis->Mesh());
        time_analysis->X_n() = time_analysis->Mesh()->Solution();
        time_analysis->AssembleResidual();
        time_analysis->Rhs() *= -1.0;
        time_analysis->Solve();
        time_analysis->Solution() += time_analysis->X_n();
        time_analysis->LoadSolution();
//        time_analysis->PostProcessStep(plotfile);
        
        if(ip%progress == 0){
            percent += 10.0;
            std::cout << " Progress on offline stage " << percent  << " % " <<std::endl;
//            std::cout << " Progress on offline stage " << setw(3) << percent << setw(2)  << " % " <<std::endl;
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vector, time_analysis->Mesh());
        galerkin_projts.AddSub(0, ip, mesh_vector[0]->Solution());
    }
    
    mesh_vector[0]->LoadSolution(galerkin_projts);
    
#ifdef PZDEBUG
    std::ofstream out("CmeshElasticity.txt");
    mesh_vector[0]->Print(out);
#endif
    
    return mesh_vector[0];
    
    
}

int DrawUnitPressuresBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures, int level){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif
    
//    TPZVec<long> drained(0,0);
//    constant_pressures.Push(drained);

    TPZGeoMesh * gmesh = cmesh->Reference();
    int dim = gmesh->Dimension();
    cmesh->LoadReferences();
    
    TPZVec<long> dof_indexes;
    
#ifdef PZDEBUG
    if(!gmesh){
        DebugStop();
    }
#endif
    
    int nel = gmesh->NElements();
    
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = gmesh->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        
        if(gel->Level() != level || gel->Dimension() != dim){
            continue;
        }
        
        // For a selectec group of elements
        
        TPZVec<TPZGeoEl *> unrefined_sons;
        gel->GetHigherSubElements(unrefined_sons);
        int nsub = unrefined_sons.size();
        TPZStack<long> dof_stack;
        for (int isub = 0; isub < nsub; isub++) {
            TPZGeoEl * subgel = unrefined_sons[isub];
            
            if (subgel->Dimension() != dim) {
                continue;
            }
            
            TPZCompEl *cel = subgel->Reference();
#ifdef PZDEBUG
            if(!cel){
                DebugStop();
            }
#endif
            
            TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
            if(!intel){
                DebugStop();
            }
#endif
                
            ElementDofIndexes(intel, dof_indexes);
            for (int i = 0;  i < dof_indexes.size(); i++) {
                dof_stack.Push(dof_indexes[i]);
            }

        }
        
        if(nsub == 0){
            TPZCompEl *cel = gel->Reference();
#ifdef PZDEBUG
            if(!cel){
                DebugStop();
            }
#endif
            
            TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
            if(!intel){
                DebugStop();
            }
#endif
            
            ElementDofIndexes(intel, dof_indexes);
            for (int i = 0;  i < dof_indexes.size(); i++) {
                dof_stack.Push(dof_indexes[i]);
            }
        }

        TPZVec<long> dofs(dof_stack);
        constant_pressures.Push(dofs);
        
    }
    
    int n_modes = constant_pressures.size();
    return n_modes;
}

int DrawingPressureBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<long> > & constant_pressures, TPZSimulationData * sim_data){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif

    TPZGeoMesh * geometry = cmesh->Reference();
    
#ifdef PZDEBUG
    if(!geometry){
        DebugStop();
    }
#endif
    
    int dim = geometry->Dimension();
    cmesh->LoadReferences();
    
    TPZStack<REAL> min_x;
    TPZStack<REAL> max_x;
    bool Outline_is_rightQ = DrawingGeometryOutline(geometry, min_x, max_x);
    
    if(!Outline_is_rightQ){
        DebugStop();
    }
    
    int ni = sim_data->Blocks()[0];
    int nj = sim_data->Blocks()[1];
    int nk = sim_data->Blocks()[2];
    
    TPZManVector<REAL,3> x0(3,0.0);
    
    TPZStack<REAL> rule_x;
    REAL dx = (max_x[0]-min_x[0])/ni;
    REAL xv;
    for (int i = 0; i < ni + 1; i++) {
        xv = REAL(i)*dx + min_x[0];
        rule_x.Push(xv);
    }
    
    TPZStack<REAL> rule_y;
    REAL dy = (max_x[1]-min_x[1])/nj;
    REAL yv;
    for (int j = 0; j < nj + 1; j++) {
        yv = REAL(j)*dy + min_x[1];
        rule_y.Push(yv);
    }
    
    TPZStack<REAL> rule_z;
    REAL dz = (max_x[2]-min_x[2])/nk;
    REAL zv;
    for (int k = 0; k < nk + 1; k++) {
        zv = REAL(k)*dz + min_x[2];
        rule_z.Push(zv);
    }
    
    if (rule_z[0] == rule_z[1] && dim == 2) {
        std::cout << "RB:: Drawing Pressure Blocks for 2D geometry " <<std::endl;
    }
    
    if (dim == 3) {
        std::cout << "RB:: Drawing Pressure Blocks for 3D geometry " <<std::endl;
        DebugStop();
    }
    
    // counting volumetric elements
    int nel = geometry->NElements();
    int n_volumes = 0;
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        
        if(gel->HasSubElement() || gel->Dimension() != dim){
            continue;
        }
        n_volumes++;
    }
    
    std::cout << "RB:: Number of volumetric elements =  " << n_volumes << std::endl;
    
    // goup elements by Cartesian Grid
    
    TPZStack< TPZStack<long> > geo_groups;
    // Check if the element belong to the cartesian block
    for (int i = 0; i < rule_x.size() - 1; i++) {
        for (int j = 0; j < rule_y.size() - 1; j++) {
            
            // for each i,j,k box
            TPZStack<long> box_group;
            
            TPZManVector<REAL,3> x_c(3,0.0);
            TPZManVector<REAL,3> par_c(dim,0.0);
            int nel = geometry->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl * gel = geometry->Element(iel);
                
#ifdef PZDEBUG
                if(!gel){
                    DebugStop();
                }
#endif
                
                if(gel->Level() != 0 || gel->Dimension() != dim){
                    continue;
                }
                
                gel->CenterPoint(gel->NSides()-1, par_c);
                gel->X(par_c, x_c);
                
                // It is inside x
                if(rule_x[i] <= x_c[0] && x_c[0] < rule_x[i+1]){
                    if(rule_y[j] <= x_c[1] && x_c[1] < rule_y[j+1]){
                        
//                        std::cout << " box x dimension " << rule_x[i] << " x " << rule_x[i+1] << std::endl;
//                        std::cout << " box y dimension " << rule_y[j] << " x " << rule_y[j+1] << std::endl;
//                        std::cout << " x_c " << x_c << std::endl;
                        
                        box_group.Push(gel->Index());
                    }
                }
                
            }
            
            if (box_group.size() == 0) {
                continue;
            }
            
            geo_groups.Push(box_group);
//            std::cout << " group of geo elements with indexes =  " << box_group << std::endl;
        }
    }
    
//    // divide groups by max number of elements, adapted case
//    int div = 2;
//    int cut_off = 2;
////    TPZStack< TPZStack<long> > geo_groups_adapted;
//    TPZStack<long> adapted_group;
//    int groups = geo_groups.size();
//    for (int ig = 0; ig < groups; ig++) {
//        adapted_group.Resize(0);
//        int n_elements = geo_groups[ig].size();
//        
//        if (n_elements >= cut_off) {
//            int n_newgroups = int(n_elements/div);
//            for (int iel = n_newgroups - 1 ; iel >= 0; iel--) {
//                adapted_group.Push(geo_groups[ig][iel]);
//                geo_groups[ig].Pop();
//            }
//            geo_groups.Push(adapted_group);
//        }
//
//    
//    }
    
    
#ifdef PZDEBUG
    if(geo_groups.size() == 0){
        DebugStop();
    }
#endif
    
    // drained response mode
    TPZVec<long> dofs(0);
    constant_pressures.Push(dofs);
    
    // Pick blocks dofs
    TPZVec<long> dof_indexes;
    TPZVec<long> igroup;
    int n_groups = geo_groups.size();
    
    // For all groups
    int group_n_vol = 0;
    for (int ig = 0; ig < n_groups; ig++) {
        igroup = geo_groups[ig];

        // For a selecteced group of elements
        TPZStack<long> dof_stack;
        for(int igel = 0; igel < igroup.size(); igel++){

            TPZGeoEl * gel = geometry->Element(igroup[igel]);
            
#ifdef PZDEBUG
            if(!gel || gel->Dimension() != dim){
                DebugStop();
            }
#endif
            
            TPZVec<TPZGeoEl *> unrefined_sons;
            gel->GetHigherSubElements(unrefined_sons);
            int nsub = unrefined_sons.size();
            for (int isub = 0; isub < nsub; isub++) {
                TPZGeoEl * subgel = unrefined_sons[isub];
                
                if (subgel->Dimension() != dim) {
                    continue;
                }
                
                TPZCompEl *cel = subgel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
                
            }
            
            if(nsub == 0){
                TPZCompEl *cel = gel->Reference();
#ifdef PZDEBUG
                if(!cel){
                    DebugStop();
                }
#endif
                
                TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace * >(cel);
#ifdef PZDEBUG
                if(!intel){
                    DebugStop();
                }
#endif
                
                group_n_vol++;
                ElementDofIndexes(intel, dof_indexes);
                for (int i = 0;  i < dof_indexes.size(); i++) {
                    dof_stack.Push(dof_indexes[i]);
                }
            }
        }
        
        TPZVec<long> dofs(dof_stack);
        constant_pressures.Push(dofs);
        
    }
    
    
    if(group_n_vol != n_volumes){
        std::cout << "RB:: Drawing Pressure Blocks left some elements out! " <<std::endl;
//        DebugStop();
    }
    
    int n_pressure_blocks = constant_pressures.size();
    if(n_pressure_blocks == 0){
        DebugStop();
    }
    
    return n_pressure_blocks;
}

bool DrawingGeometryOutline(TPZGeoMesh * gmesh, TPZStack< REAL > & min_x, TPZStack< REAL > & max_x){
    
#ifdef PZDEBUG
    if(!gmesh){
        DebugStop();
    }
#endif
    
    long n_nodes = gmesh->NNodes();
    REAL min_xc = +1.0e12;
    REAL max_xc = -1.0e12;
    REAL min_yc = +1.0e12;
    REAL max_yc = -1.0e12;
    REAL min_zc = +1.0e12;
    REAL max_zc = -1.0e12;
    
    TPZManVector<REAL,3> co(3,0.0);
    for (long inode = 0; inode < n_nodes; inode++) {
        TPZGeoNode node = gmesh->NodeVec()[inode];
        node.GetCoordinates(co);
        
        // x limits
        if (min_xc > co[0]) {
            min_xc = co[0];
        }
        
        if (max_xc < co[0]) {
            max_xc = co[0];
        }
        
        // y limits
        if (min_yc > co[1]) {
            min_yc = co[1];
        }
        
        if (max_yc < co[1]) {
            max_yc = co[1];
        }
        
        // z limits
        if (min_zc > co[2]) {
            min_zc = co[2];
        }
        
        if (max_zc < co[2]) {
            max_zc = co[2];
        }
        
    }
    
    min_x.Push(min_xc);
    min_x.Push(min_yc);
    min_x.Push(min_zc);
    
    max_x.Push(max_xc);
    max_x.Push(max_yc);
    max_x.Push(max_zc);

    return true;
}


void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

TPZCompMesh * CMesh_Pressures(TPZGeoMesh * gmesh){
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Setting piecewise constant approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(0);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshUnitPressures.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * CMesh_Elasticity(TPZGeoMesh * gmesh, int order){
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
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
    
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshElasticity.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

// Create a computational mesh for basis generation multiphysisc version
TPZCompMesh * CMesh_GeoModes_M(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    int dim = gmesh->Dimension();
    
    REAL MPa = 1.0e6;
    
    // soil parameters
    // http://www.sciencedirect.com/science/article/pii/S0045782505001532
    
    REAL l          = 8.333e3;
    REAL mu         = 12.50e3;
    REAL l_u        = 8.333e3;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL k          = 1.0e-10;
    REAL porosity   = 0.25;
    REAL eta        = 0.001;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
   
    TPZPoroelasticModes * material = new TPZPoroelasticModes(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPorolasticParameters(l, mu, l_u);
    material->SetBiotParameters(alpha, Se);
    material->SetFlowParameters(k, porosity, eta);
    cmesh->InsertMaterialObject(material);
    
    
    // Inserting boundary conditions
    int dirichlet_x_vn   = 7;
    int dirichlet_xy_vn  = 6;
    int neumann_y_p      = 5;
    
    REAL s_n = -(1.0e-3)*MPa;
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_bottom_mat = material->CreateBC(material, bc_bottom, dirichlet_xy_vn, val1, val2);
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
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_x_vn, val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, neumann_y_p, val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
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
    std::ofstream out("CMeshGeoModesMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
    
}

REAL Error_Benchmark(TPZCompMesh * cmesh_rb, TPZCompMesh * cmesh_ref){
    
#ifdef PZDEBUG
    
    if(!cmesh_rb){
        DebugStop();
    }
    
    if(!cmesh_ref){
        DebugStop();
    }
    
#endif
    
    REAL error  = 0.0;
    
    int nel_rb = cmesh_rb->NElements();
    int nel_re = cmesh_ref->NElements();
    
#ifdef PZDEBUG
    
    if(nel_rb != nel_re){
        DebugStop();
    }
    
#endif
    int side_gel = 0;
    int order = 4;
    
    int var_index = 1; // u
    TPZManVector<REAL,1> u_rb;
    TPZManVector<REAL,1> u_re;
    
    TPZVec<REAL> pos(3,0.0);
    REAL w, el_error;
    
    TPZFMatrix<REAL> jac,axes,jacinv;
    REAL detjac;
    for (int icel = 0 ; icel < nel_rb; icel++) {
        TPZCompEl * cel_rb = cmesh_rb->Element(icel);
        TPZCompEl * cel_re = cmesh_ref->Element(icel);
        
        TPZGeoEl *gel  = cel_re->Reference();
        if (gel->Dimension()  != 2) {
            continue;
        }
        side_gel = gel->NSides() - 1;
        TPZIntPoints * int_points = gel->CreateSideIntegrationRule(side_gel, order);
        el_error = 0.0;
        int npoints = int_points->NPoints();
        for (int i = 0; i < npoints; i++) {
            int_points->Point(i, pos, w);
            gel->Jacobian(pos, jac, axes, detjac, jacinv);
            cel_rb->Solution(pos, var_index, u_rb);
            cel_re->Solution(pos, var_index, u_re);
            el_error += w * detjac * (u_rb[0] - u_re[0])*(u_rb[0] - u_re[0]);
        }
        
        error += el_error;
        
    }
    
    return error;
    
}

TPZCompMesh * Reference_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data){
    
    // Create the approximation space
    int potential_order = 2;
    
    // Create multiphysisc mesh
    TPZVec<TPZCompMesh * > mesh_vector(1);
    mesh_vector[0] = CMesh_Elliptic(gmesh, potential_order);
    
    TPZCompMesh * nonlinear_cmesh = CMesh_Elliptic_M(gmesh, mesh_vector, sim_data);
    
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 16;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(nonlinear_cmesh,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    //    TPZSkylineNSymStructMatrix struct_mat(nonlinear_cmesh);
    TPZSkylineStructMatrix struct_mat(nonlinear_cmesh);
    
    //    TPZSymetricSpStructMatrix struct_mat(nonlinear_cmesh);
    //    struct_mat.SetNumThreads(number_threads);
    
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(nonlinear_cmesh);
//    struct_mat.SetDecomposeType(ELDLt);
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = nonlinear_cmesh->NEquations();
    std::cout<< "ndof = " << ndof << std::endl;
    
    // Reference Solution
    std::string plotfile_ref("Nonlinear_Elliptic_ref.vtk");
    
    time_analysis->ExcecuteOneStep();
    time_analysis->PostNonlinearProcessStep(plotfile_ref);
    time_analysis->PostNonlinearProcessStep(plotfile_ref);

    return nonlinear_cmesh;
    
}

TPZCompMesh * OffLine_Benchmark(TPZGeoMesh * gmesh, TPZSimulationData * sim_data, TPZVec<TPZCompMesh * > & mesh_vector, int n, TPZFMatrix<REAL> &G_projts){
    
    // Create the approximation space
    int potential_order = 2;
    
    // Create multiphysisc mesh
    mesh_vector.Resize(1);
    mesh_vector[0] = CMesh_Elliptic(gmesh, potential_order);
    
    TPZCompMesh * nonlinear_cmesh = CMesh_Elliptic_M(gmesh, mesh_vector, sim_data);
    
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 16;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(nonlinear_cmesh,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
    //    TPZSkylineNSymStructMatrix struct_mat(nonlinear_cmesh);
    TPZSkylineStructMatrix struct_mat(nonlinear_cmesh);
    
//    TPZSymetricSpStructMatrix struct_mat(nonlinear_cmesh);
//    struct_mat.SetNumThreads(number_threads);
    
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
    std::string plotfile("Nonlinear_Elliptic.vtk");
    TPZManVector<REAL,3> mu_vector(3,0.0);
    
    TPZManVector<REAL,3> mu_1(3,0.0),mu_2(3,0.0),mu_3(3,0.0);
    mu_1[0] = 0.1;
    mu_1[1] = 1.0;
    mu_1[2] = M_PI;
    
    mu_2[0] = 0.1;
    mu_2[1] = 1.0;
    mu_2[2] = M_PI;
    
    mu_2[0] = 0.1;
    mu_2[1] = 1.0;
    mu_2[2] = M_PI;
    
    int n0 = mu_1.size();
    int n1 = mu_2.size();
    int n2 = mu_3.size();
    
    G_projts.Resize(ndof,n0*n1*n2);
    G_projts.Zero();
    
    int count = 0;
    for (int i = 0; i < n0; i++) {
        mu_vector[0] = mu_1[i];
        for (int j = 0; j < n1; j++) {
            mu_vector[1] = mu_2[j];
            for (int k = 0; k < n2; k++) {
                mu_vector[2] = mu_3[k];
                SetParameters(time_analysis->Mesh(), mu_vector);
                time_analysis->ExcecuteOneStep();
                time_analysis->PostNonlinearProcessStep(plotfile);
                G_projts.AddSub(0, count, time_analysis->X_n());
                count++;
            }
        }
    }

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

TPZCompMesh * OnLine_Benchmark(TPZCompMesh * cmesh, TPZSimulationData * sim_data){
    
    // Create the approximation space
    TPZGeoMesh * gmesh = cmesh->Reference();
    
    std::string plotfile("Nonlinear_Elliptic_rb.vtk");
    
    // Create multiphysisc mesh
    TPZVec<TPZCompMesh * > mesh_vector(1);
    mesh_vector[0] = CMesh_Elliptic_RB(cmesh);

    TPZCompMesh * nonlinear_rb_cmesh = CMesh_Elliptic_M_RB(gmesh, mesh_vector, sim_data);
    
    bool mustOptimizeBandwidth = false;
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
    time_analysis->PostNonlinearProcessStep(plotfile);
//    time_analysis->X_n().Print("Selected modes = ");
    
    return nonlinear_rb_cmesh;
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
    
    REAL mu_0 = 1.0;
    REAL mu_1 = M_PI/10000.0;
    REAL mu_2 = 4.0*M_PI;
    
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
    
    REAL mu_0 = 1.0;
    REAL mu_1 = M_PI/10000000.0;
    REAL mu_2 = 4.0*M_PI;
    
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
        
//        if(mfcel->Dimension()){
//            continue;
//        }
        
        mfcel->InitializeIntegrationRule();       
        mfcel->PrepareIntPtIndices();
    }
    
#ifdef PZDEBUG
    std::ofstream out("CMeshEllipticMultiPhysicsRB.txt");
    cmesh->Print(out);
#endif

    return cmesh;
    
}

TPZCompMesh * CMesh_Elliptic(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
    
    REAL MPa = 1.0e6;
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    // soil parameters
    // http://www.sciencedirect.com/science/article/pii/S0045782505001532
    
    REAL l          = 8.333e3;
    REAL mu         = 12.50e3;
    REAL l_u        = 8.333e3;
    REAL l_qin      = 5.0e11;
    REAL mu_qin     = 10000.0;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL phi        = 0.25;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    
    TPZElasticBiot * material = new TPZElasticBiot(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPorolasticParameters(l, mu, l_u, l_qin, mu_qin);
    material->SetBiotParameters(alpha, Se, phi);
    
    TPZAutoPointer<TPZFunction<STATE> > f_analytic = new TPZDummyFunction<STATE>(Analytic);
    material->SetTimeDependentForcingFunction(f_analytic);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet_v   = 0;
    int dirichlet_x   = 1;
    int neumann_y     = 5;
    
    REAL s_n = -(1.0e-3)*MPa;
    
    TPZFMatrix<STATE> val1(dim,dim,0.), val2(dim,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * bc_bottom_mat = new TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>;
    bc_bottom_mat->SetNumLoadCases(1);
    bc_bottom_mat->SetMaterial(material);
    bc_bottom_mat->SetId(bc_bottom);
    bc_bottom_mat->SetType(dirichlet_v);
    bc_bottom_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * bc_right_mat = new TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>;
    bc_right_mat->SetNumLoadCases(1);
    bc_right_mat->SetMaterial(material);
    bc_right_mat->SetId(bc_right);
    bc_right_mat->SetType(dirichlet_x);
    bc_right_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = s_n;
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * bc_top_mat = new TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>;
    bc_top_mat->SetNumLoadCases(1);
    bc_top_mat->SetMaterial(material);
    bc_top_mat->SetId(bc_top);
    bc_top_mat->SetType(neumann_y);
    bc_top_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * bc_left_mat = new TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>;
    bc_left_mat->SetNumLoadCases(1);
    bc_left_mat->SetMaterial(material);
    bc_left_mat->SetId(bc_left);
    bc_left_mat->SetType(dirichlet_x);
    bc_left_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond> * bc_top_null_mat = new TPZMatWithMem<TPZElasticBiotMemory,TPZBndCond>;
    bc_top_null_mat->SetNumLoadCases(1);
    bc_top_null_mat->SetMaterial(material);
    bc_top_null_mat->SetId(bc_top_null);
    bc_top_null_mat->SetType(neumann_y);
    bc_top_null_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);
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

TPZCompMesh * CMesh_Parabolic(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data, bool IsMixedQ){
    
//    REAL MPa = 1.0e6;
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    // soil parameters
    // http://www.sciencedirect.com/science/article/pii/S0045782505001532
    
    REAL k          = 1.0e-10;
    REAL eta        = 0.001;
    
    REAL l          = 8.333e3;
    REAL mu         = 12.50e3;
    REAL l_u        = 8.333e3;
    REAL l_qin      = 5.0e11;
    REAL mu_qin     = 10000.0;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL phi        = 0.25;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    
    TPZDarcyFlow * material = new TPZDarcyFlow(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetFlowParameters(k, eta);
    material->SetMixedFormulation(IsMixedQ);
    material->SetPorolasticParameters(l, mu, l_u, l_qin, mu_qin);
    material->SetBiotParameters(alpha, Se, phi);
    
    TPZAutoPointer<TPZFunction<STATE> > f_analytic = new TPZDummyFunction<STATE>(Analytic);
    material->SetTimeDependentForcingFunction(f_analytic);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet   = 0;
    int neumann     = 1;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    
    val2(0,0) = 0.0;
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * bc_bottom_mat = new TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>;
    bc_bottom_mat->SetNumLoadCases(1);
    bc_bottom_mat->SetMaterial(material);
    bc_bottom_mat->SetId(bc_bottom);
    bc_bottom_mat->SetType(neumann);
    bc_bottom_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    val2(0,0) = 0.0;
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * bc_right_mat = new TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>;
    bc_right_mat->SetNumLoadCases(1);
    bc_right_mat->SetMaterial(material);
    bc_right_mat->SetId(bc_right);
    bc_right_mat->SetType(neumann);
    bc_right_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * bc_top_mat = new TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>;
    bc_top_mat->SetNumLoadCases(1);
    bc_top_mat->SetMaterial(material);
    bc_top_mat->SetId(bc_top);
    bc_top_mat->SetType(dirichlet);
    bc_top_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * bc_left_mat = new TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>;
    bc_left_mat->SetNumLoadCases(1);
    bc_left_mat->SetMaterial(material);
    bc_left_mat->SetId(bc_left);
    bc_left_mat->SetType(neumann);
    bc_left_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    val2(0,0) = 0.0;
    TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond> * bc_top_null_mat = new TPZMatWithMem<TPZDarcyFlowMemory,TPZBndCond>;
    bc_top_null_mat->SetNumLoadCases(1);
    bc_top_null_mat->SetMaterial(material);
    bc_top_null_mat->SetId(bc_top_null);
    bc_top_null_mat->SetType(neumann);
    bc_top_null_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);
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
    std::ofstream out("CMeshParabolicMultiPhysics.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}



TPZCompMesh * CMesh_GeomechanicCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data, bool IsMixedQ){
    
    
    REAL MPa = 1.0e6;
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    // soil parameters
    // http://www.sciencedirect.com/science/article/pii/S0045782505001532
    
    REAL l          = 8.333e3;
    REAL mu         = 12.50e3;
    REAL l_u        = 8.333e3;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL k          = 1.0e-10;
    REAL porosity   = 0.25;
    REAL eta        = 0.001;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZBiotPoroelasticity * material = new TPZBiotPoroelasticity(matid,dim);
    material->SetSimulationData(sim_data);
    material->SetPorolasticParameters(l, mu, l_u);
    material->SetBiotParameters(alpha, Se);
    material->SetFlowParameters(k, porosity, eta);
    material->SetMixedFormulation(IsMixedQ);
    material->SetSymmetricFormulation(true);
    
    TPZAutoPointer<TPZFunction<STATE> > f_analytic = new TPZDummyFunction<STATE>(Analytic);
    material->SetTimeDependentForcingFunction(f_analytic);
    cmesh->InsertMaterialObject(material);
    
    // Inserting boundary conditions
    int dirichlet_x_vn   = 7;
    int dirichlet_xy_vn   = 6;
    int neumann_y_p      = 5;
    int neumann_y_vn     = 11;

    REAL s_n = -(1.0e-3)*MPa;
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * bc_bottom_mat = new TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>;
    bc_bottom_mat->SetNumLoadCases(1);
    bc_bottom_mat->SetMaterial(material);
    bc_bottom_mat->SetId(bc_bottom);
    bc_bottom_mat->SetType(dirichlet_xy_vn);
    bc_bottom_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_bottom_mat);
    
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * bc_right_mat = new TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>;
    bc_right_mat->SetNumLoadCases(1);
    bc_right_mat->SetMaterial(material);
    bc_right_mat->SetId(bc_right);
    bc_right_mat->SetType(dirichlet_x_vn);
    bc_right_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_right_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = s_n;
    val2(2,0) = 0.0;
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * bc_top_mat = new TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>;
    bc_top_mat->SetNumLoadCases(1);
    bc_top_mat->SetMaterial(material);
    bc_top_mat->SetId(bc_top);
    bc_top_mat->SetType(neumann_y_p);
    bc_top_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * bc_left_mat = new TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>;
    bc_left_mat->SetNumLoadCases(1);
    bc_left_mat->SetMaterial(material);
    bc_left_mat->SetId(bc_left);
    bc_left_mat->SetType(dirichlet_x_vn);
    bc_left_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_left_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * bc_top_null_mat = new TPZMatWithMem<TPZPoroPermMemory,TPZBndCond>;
    bc_top_null_mat->SetNumLoadCases(1);
    bc_top_null_mat->SetMaterial(material);
    bc_top_null_mat->SetId(bc_top_null);
    bc_top_null_mat->SetType(neumann_y_vn);
    bc_top_null_mat->SetValues(val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->ApproxSpace().CreateWithMemory(true);
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

// Create a computational mesh for H1 pore pressure;
TPZCompMesh * CMesh_PorePressure(TPZGeoMesh * gmesh, int order){
    
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
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
    
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
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
TPZCompMesh * CMesh_Flux(TPZGeoMesh * gmesh, int order){
    

    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
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
    
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    // Setting Hdiv approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();

    
#ifdef PZDEBUG
    std::ofstream out("CmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

// Create a computational mesh for mixed pore pressure;
TPZCompMesh * CMesh_MFPorePressure(TPZGeoMesh * gmesh, int order){
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    // Aproximation Space of order -> pOrder
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    
    // Creating a material object
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
    cmesh->InsertMaterialObject(material);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(order);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshMFPorePressure.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_Deformation_rb(TPZCompMesh * cmesh, int order){
        
    // Plane strain assumption
//    int planestress = 0;
    TPZGeoMesh * gmesh = cmesh->Reference();
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
    TPZCompMeshReferred * cmesh_rb = new TPZCompMeshReferred(gmesh);


    // Creating a material object
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
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
    
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, dirichlet, val1, val2);
    cmesh_rb->InsertMaterialObject(bc_top_null_mat);
    
    
    // Setting RB approximation space
    cmesh_rb->SetDimModel(dim);
    int numsol = cmesh->Solution().Cols();
    cmesh_rb->AllocateNewConnect(numsol, 1, order);
    TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmesh_rb);
    cmesh_rb->AutoBuild();
    
    cmesh_rb->AdjustBoundaryElements();
    cmesh_rb->CleanUpUnconnectedNodes();
    cmesh_rb->LoadReferred(cmesh);
    
#ifdef PZDEBUG
    std::ofstream out("CmeshDeformationRB.txt");
    cmesh_rb->Print(out);
#endif
    
    return cmesh_rb;
    
}

// Create a computational mesh for deformation;
TPZCompMesh * CMesh_Deformation(TPZGeoMesh * gmesh, int order){
    
    
    // Getting mesh dimension
    int dim = gmesh->Dimension();
    
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
    
    TPZMaterial * bc_top_null_mat = material->CreateBC(material, bc_top_null, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bc_top_null_mat);
    
    
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

/** @brief Create a reservoir-box geometry with cylindrical wells */
TPZGeoMesh * CreateGeometricGmshMesh(std::string &grid){
    
    TPZGeoMesh * geometry = new TPZGeoMesh;
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    geometry = Geometry.GeometricGmshMesh(grid);
    const std::string name("Reduced base of geomechanic coupling");
    geometry->SetName(name);
    return geometry;
}


TPZGeoMesh * RockBox(TPZVec<REAL> dx_dy, TPZVec<int> n, bool IsTriangleMeshQ){
    
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
    if (IsTriangleMeshQ) {
        CreateGridFrom2.SetTriangleExtrusion();
    }

    dy = dx_dy[1];
    n_elements = n[1];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dy, n_elements);
    
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

void Analytic(const TPZVec<REAL> &x, REAL time, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu){
    
    u.resize(5);
    gradu.Resize(3, 1);
    gradu.Zero();
    
    REAL y_c = x[1];
    
    REAL kappa = 1.0e-10;
    REAL eta = 0.001;
    REAL l = 8.333e3;
    REAL mu = 12.50e3;
    REAL h = 1.0;
    REAL p0 = 1.0e3;

    REAL yD = (h-(y_c+h/2))/h;
    REAL tD = (l+2.0*mu)*kappa*time/(eta*h*h);
    
    int n = 1000;
    REAL sump = 0.0;
    REAL sumu = 0.0;
    REAL sumv = 0.0;
    REAL M;
    for (int k = 0; k < n; k++) {
        M = 0.5*M_PI*REAL(2*k+1);
        sump += (2.0/M)*sin(M*yD)*exp(-M*M*tD);
        sumu += (2.0/(M*M))*cos(M*yD)*exp(-M*M*tD);
        sumv += (-2.0*cos(M*yD)*exp(-M*M*tD));
    }

    u[0] = 0.0;
    u[1] = -((p0*h)/(l + 2.0*mu))*(1.0 - yD - sumu);
    u[2] = p0*sump;
    u[3] = 0.0;
    u[4] = -1.0*(kappa/eta)*(p0/h)*sumv;
}

