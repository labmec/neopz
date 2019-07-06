
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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
#include "pzelast3d.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZPoroPermCoupling.h"
#include "TPZNonLinearElliptic.h"
#include "TPZLinearElliptic.h"

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
void ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X);

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
TPZCompMesh * CMesh_GeomechanicCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data);

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

// Compute Galerkin projections of unit pressures
TPZCompMesh * Galerkin_Projections(TPZGeoMesh * gmesh, TPZSimulationData * sim_data, int order, int level);

int DrawUnitPressuresBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<int64_t> > & constant_pressures, int level);

void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes);

// Create a computational mesh for reduced deformation
TPZCompMesh * CMesh_Deformation_rb(TPZCompMesh * cmesh);


void SetParameters(TPZCompMesh * cmesh, TPZVec<REAL> mu_vector);

int NonLinearElliptic();

int Geomechanic();

int main(int argc, char *argv[])
{
    
//    NonLinearElliptic();
    
    Geomechanic();
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
    REAL epsilon_res = 1.0e-4;
    REAL epsilon_corr = 1.0e-8;
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
    
    TPZGeoMesh * gmesh = RockBox(dx_dy,n);
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
    int n_steps = 10;
    REAL epsilon_res = 1.0e-4;
    REAL epsilon_corr = 1.0e-7;
    int n_corrections = 10;
    
    /** @brief Definition gravity field */
    TPZVec<REAL> g(2,0.0);
    g[1] = -0.0*9.81;
    
    sim_data->SetGravity(g);
    sim_data->SetTimeControls(n_steps, dt);
    sim_data->SetNumericControls(n_corrections, epsilon_res, epsilon_corr);
    
    TPZVec<REAL> dx_dy(2);
    TPZVec<int> n(2);
    
    REAL Lx = 10.0; // meters
    REAL Ly = 10.0; // meters
    
    n[0] = 5; // x - direction
    n[1] = 5; // y - direction
    
    int order = 2;
    int level = 0;
    int hlevel = 2;
    
    dx_dy[0] = Lx/REAL(n[0]); // x - direction
    dx_dy[1] = Ly/REAL(n[1]); // y - direction
    
    TPZGeoMesh * gmesh = RockBox(dx_dy,n);
    UniformRefinement(gmesh, hlevel);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("Geometry.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("Geometry.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    std::cout<< "Geometry done. " << std::endl;

    TPZCompMesh * cmesh_gp = Galerkin_Projections(gmesh, sim_data, order,level);
    
//    order = 2;
    // Computing reference solution
    TPZVec<TPZCompMesh * > mesh_vector(2);
//    mesh_vector[0] = CMesh_Deformation(gmesh, order);
    mesh_vector[0] = CMesh_Deformation_rb(cmesh_gp);
    mesh_vector[1] = CMesh_PorePressure(gmesh, order-1);
    TPZCompMesh * geomechanic = CMesh_GeomechanicCoupling(gmesh, mesh_vector, sim_data);
    
    bool mustOptimizeBandwidth = false;
    int number_threads = 16;
    TPZGeomechanicAnalysis * time_analysis = new TPZGeomechanicAnalysis;
    time_analysis->SetCompMesh(geomechanic,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
//    TPZSkylineNSymStructMatrix struct_mat(geomechanic);
//    TPZSkylineStructMatrix struct_mat(geomechanic);

//    TPZSymetricSpStructMatrix struct_mat(geomechanic);
//    struct_mat.SetNumThreads(number_threads);
    
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(geomechanic);
    struct_mat.SetDecomposeType(ELDLt);

    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = geomechanic->NEquations();
//    std::cout << " Number of modes =  " << cmesh_gp->Solution().Cols() << std::endl;
    std::cout << " Full order model ndof =  " << ndof << std::endl;
    
    TPZVec<REAL> x(3);
    x[0] = Lx/2.0;
    x[1] = Ly/2.0;
    x[2] = 0.0;
    std::string plotfile("geomechanic_rb_0.vtk");
    // Run Transient analysis
    time_analysis->Run_Evolution(x,plotfile);
    
    std::cout << " Execution finished " << std::endl;
    return EXIT_SUCCESS;
    
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
    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    int ndof = geo_modes->NEquations();
    
    std::cout<< "ROM:: strategy ndof = " << ndof << std::endl;
    
    time_analysis->SimulationData()->SetCurrentStateQ(true);
    time_analysis->Assemble();
    
    // Setting up the empirical interpolation based on unitary pressures
    std::string plotfile("Geo_Modes_rb_0.vtk");
    
    REAL unit_p = 1.0e6;
    TPZStack<TPZVec<int64_t> > cts_pressures;
    
    int n_blocks = DrawUnitPressuresBlocks(mesh_vector[1],cts_pressures,level);
    int ndof_elastic = mesh_vector[0]->NEquations();
    TPZFMatrix<REAL> galerkin_projts(ndof_elastic,n_blocks);
    galerkin_projts.Zero();
    
    std::cout<< "ROM:: number of geomodes = " << n_blocks << std::endl;
    for (int ip = 0; ip < n_blocks; ip++) {
        
        mesh_vector[0]->Solution().Zero();
        mesh_vector[1]->Solution().Zero();
        
        for (int jp = 0; jp < cts_pressures[ip].size(); jp++) {
            mesh_vector[1]->Solution()(cts_pressures[ip][jp],0) = unit_p;
        }

        TPZBuildMultiphysicsMesh::TransferFromMeshes(mesh_vector, time_analysis->Mesh());
        time_analysis->X_n() = time_analysis->Mesh()->Solution();
        time_analysis->AssembleResidual();
        time_analysis->Solve();
        time_analysis->Solution() += time_analysis->X_n();
        time_analysis->LoadSolution();
//#ifdef PZDEBUG
        time_analysis->PostProcessStep(plotfile);
//#endif
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vector, time_analysis->Mesh());
        galerkin_projts.AddSub(0, ip, mesh_vector[0]->Solution());
    }
    
    mesh_vector[0]->LoadSolution(galerkin_projts);
    return mesh_vector[0];
    
}

int DrawUnitPressuresBlocks(TPZCompMesh * cmesh, TPZStack<TPZVec<int64_t> > & constant_pressures, int level){
    
#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif

    TPZGeoMesh * gmesh = cmesh->Reference();
    int dim = gmesh->Dimension();
    cmesh->LoadReferences();
    
    TPZVec<int64_t> dof_indexes;
    
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
        
        TPZVec<TPZGeoEl *> unrefined_sons;
        gel->GetHigherSubElements(unrefined_sons);
        int nsub = unrefined_sons.size();
        TPZStack<int64_t> dof_stack;
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

        TPZVec<int64_t> dofs(dof_stack);
        constant_pressures.Push(dofs);
        
    }
    
    int n_modes = constant_pressures.size();
    return n_modes;
}

void ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<int64_t> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<int64_t> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        int64_t seqnumber = con.SequenceNumber();
        int64_t position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

TPZCompMesh * CMesh_Pressures(TPZGeoMesh * gmesh){
    
    // Material identifiers
    int matid =1;
    // Getting mesh dimension
    int dim = 2;
    
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
    std::ofstream out("CmeshElasticity.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

// Create a computational mesh for basis generation multiphysisc version
TPZCompMesh * CMesh_GeoModes_M(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
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
    
    int kmodel      = 0;
    REAL l          = 40.38e9;
    REAL mu         = 26.92e9;
    REAL l_u        = 40.38e9;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL k          = 1.0e-14;
    REAL porosity   = 1.0;
    REAL eta        = 0.001;
    
    REAL c = 0.0;
    REAL phi_f = 0.0;
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    
    // Creating a material object
    TPZLinearElliptic * material = new TPZLinearElliptic(matid,dim);
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
    val2(2,0) = 0.0;
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
    
    
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
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
        int64_t n = gmesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (sons);
        }//for i
    }//ref
}

void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> sons;
        int64_t n = gmesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
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
    
    
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
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
    
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
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



TPZCompMesh * CMesh_GeomechanicCoupling(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > mesh_vector, TPZSimulationData * sim_data){
    
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
    
    int kmodel      = 0;
    REAL l          = 40.38e9;
    REAL mu         = 26.92e9;
    REAL l_u        = 40.38e9;
    REAL alpha      = 1.0;
    REAL Se         = 0.0;
    REAL k          = 1.0e-14;
    REAL porosity   = 1.0;
    REAL eta        = 0.001;
    
    REAL c = 0.0;
    REAL phi_f = 0.0;

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
    
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
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
    cmesh->InsertMaterialObject(bc_top_mat);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZMaterial * bc_left_mat = material->CreateBC(material, bc_left, dirichlet_x_vn, val1, val2);
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
    
    
    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> indices;
    for (int64_t el = 0; el<nel; el++) {
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
TPZCompMesh * CMesh_Deformation_rb(TPZCompMesh * cmesh){
        
    // Plane strain assumption
//    int planestress = 0;
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
    
    
    // Setting RB approximation space
    cmesh_rb->SetDimModel(dim);
    int numsol = cmesh->Solution().Cols();
    cmesh_rb->AllocateNewConnect(numsol, 2, 1);
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
    
    TPZVec<int64_t> Topology(1,0);
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
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(ParametricfunctionX);
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
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(ParametricfunctionY);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(bc_bottom,bc_top);
    dy = dx_dy[1];
    n_elements = n[1];
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3 = CreateGridFrom2.ComputeExtrusion(t, dy, n_elements);
    
    int64_t last_node = GeoMesh3->NNodes() - 1;
    int64_t last_element = GeoMesh3->NElements() - 1;
    int64_t node_id = GeoMesh3->NodeVec()[last_node].Id();
    int64_t element_id = GeoMesh3->Element(last_element)->Id();
    const std::string name("Geomechanic Reservoir box");
    GeoMesh3->SetName(name);
    GeoMesh3->SetMaxNodeId(node_id);
    GeoMesh3->SetMaxElementId(element_id);
    GeoMesh3->SetDimension(2);
    return GeoMesh3;
    
}

void ParametricfunctionX(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0]; // x
    X[1] = 0.0; // y
    X[2] = 0.0; // z
}

void ParametricfunctionY(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0; // x
    X[1] = par[0]; // y
    X[2] = 0.0; // z
}

