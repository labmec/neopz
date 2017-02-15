//
//  Hdiv3DCurved.cpp
//  Publication about 3D Curved piola mapping for mixed formulation
//
//  Created by omar on 01/07/2016.
//
//

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"

#include "tpzhierarquicalgrid.h"
#include "TPZReadGIDGrid.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "tpzquadraticquad.h"
#include "tpzarc3d.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzfunction.h"
#include "tpzchangeel.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "TPZPrimalPoisson.h"
#include "TPZDualPoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompMeshTools.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckmesh.h"
#include "TPZGmshReader.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


struct SimulationCase {
    bool            IsHdivQ;
    bool            IsMHMQ;
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    bool            UseGmshMeshQ;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             n_acc_terms;
    int             int_order;
    int             n_threads;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   gamma_ids;
    
    SimulationCase() : IsHdivQ(false), IsMHMQ(false), UsePardisoQ(true), UseFrontalQ(false), UseGmshMeshQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0), mesh_type(""),
    domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),gamma_ids()
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) : IsHdivQ(copy.IsHdivQ), IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
        UseGmshMeshQ(copy.UseGmshMeshQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),
        n_threads(copy.n_threads), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
        dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), gamma_ids(copy.gamma_ids)
    {
        
    }
    
    SimulationCase &operator=(const SimulationCase &copy)
    {
        IsHdivQ = copy.IsHdivQ;
        IsMHMQ = copy.IsMHMQ;
        UsePardisoQ = copy.UsePardisoQ;
        UseFrontalQ = copy.UseFrontalQ;
        UseGmshMeshQ = copy.UseGmshMeshQ;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        n_acc_terms = copy.n_acc_terms;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        mesh_type = copy.mesh_type;
        domain_type = copy.domain_type;
        conv_summary = copy.conv_summary;
        dump_folder = copy.dump_folder;
        omega_ids = copy.omega_ids;
        gamma_ids = copy.gamma_ids;
        return *this;
    }
};

//#define Solution1
//#define Solution2
//#define Solution3
#define Solution6
//#define Thiem

// MHM rates subtructuring level
static int level_mhm = 0;

static void Analytic(const TPZVec<REAL> &x, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu);
static void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &f);
static void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf);

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase  & sim_data);
void PrintGeometry(TPZGeoMesh * gmesh, SimulationCase & sim_data);
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref);
void UniformRefineTetrahedrons(TPZGeoMesh * gmesh, int n_ref);


TPZGeoMesh * MakeSphereFromLinearQuadrilateralFaces(int ndiv, SimulationCase  & sim_data);
TPZGeoMesh * MakeSphereFromQuadrilateralFaces(int ndiv, SimulationCase  & sim_data);

TPZGeoMesh * MakeSphereFromQuadrilateralFacesR(int ndiv, SimulationCase  & sim_data);

// Cylinder
TPZGeoMesh * MakeCylinderFromLinearFaces(int ndiv, SimulationCase & sim_data);
TPZGeoMesh * ExtrudedGIDMesh(TPZGeoMesh * gmesh, SimulationCase sim_data, TPZManVector<REAL,2> dz);
void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);

TPZManVector<STATE,3> ParametricSphere(REAL radius,REAL theta,REAL phi);

TPZManVector<STATE,3> ParametricCylinder(REAL radius,REAL theta,REAL z);

void TransformToQuadratic(TPZGeoMesh *gmesh);
void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);

TPZCompMesh *DualMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof, TPZVec<TPZCompMesh *> &meshvec);

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof); //  Primal approximation
TPZCompMesh * qMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); // Hdiv space
TPZCompMesh * pMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data); // L2 space

/// uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh);

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus);

/// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh);


TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
void PosProcess(TPZAnalysis * an, std::string file, SimulationCase & sim_data);

void ComputeCases(TPZStack<SimulationCase> cases);
void ComputeApproximation(SimulationCase & sim_data);
void ComputeConvergenceRates(TPZVec<STATE> &error, TPZVec<STATE> &convergence);

STATE IntegrateVolume(TPZGeoMesh * geometry, SimulationCase sim_data);
STATE IntegrateSolution(TPZCompMesh * cmesh,  SimulationCase sim_data);

// MHM utilities
/** @brief Sparated connects by hdiv connect neighborhood */
void SeparateConnectsByNeighborhood(TPZCompMesh * mixed_cmesh);

/** @brief Insert low dimentional elements defining the skeleton */
void InsertSkeletonInterfaces(TPZGeoMesh * gmesh);

/** @brief Build mhm macro elements following the mixed sense (space constrains) */
void BuildMacroElements(TPZCompMesh * mixed_cmesh);

void ErrorH1(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_h1);
void ErrorHdiv(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_hdiv);

int main()
{

  HDivPiola = 1;
    
  gRefDBase.InitializeAllUniformRefPatterns();
    //	gRefDBase.InitializeRefPatterns();
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    TPZStack<SimulationCase> simulations;
    
    // Formulations over the sphere 
    struct SimulationCase common;
    common.UsePardisoQ = true;
    common.UseFrontalQ = false;
    common.IsMHMQ      = false;
    common.UseGmshMeshQ = true;
    common.n_h_levels = 4;
    common.n_p_levels = 2;
    common.int_order  = 8;
    common.n_threads  = 16;
    common.domain_type = "sphere";
    common.conv_summary = "convergence_summary";
    common.omega_ids.Push(1);     // Domain
    common.gamma_ids.Push(-1);    // Gamma_D outer surface
    common.gamma_ids.Push(-2);    // Gamma_D inner surface
    
    //    Case_XXX.elemen_type = {0,1,2} as follows:
    //    Case_XXX.elemen_type = 0; -> Tetra
    //    Case_XXX.elemen_type = 1; -> Prism
    //    Case_XX.elemen_type  = 2; -> Hexa
    
    
//    // Primal Formulation over the solid sphere
//    struct SimulationCase H1Case_1 = common;
//    H1Case_1.IsHdivQ = false;
//    H1Case_1.mesh_type = "linear";
//    H1Case_1.elemen_type = 2;
//    H1Case_1.dump_folder = "H1_sphere";
//    simulations.Push(H1Case_1);

//    // Primal Formulation over the solid sphere
//    struct SimulationCase H1Case_2 = common;
//    H1Case_2.IsHdivQ = false;
//    H1Case_2.mesh_type = "quadratic";
//    H1Case_2.elemen_type = 0;
//    H1Case_2.dump_folder = "H1_sphere_T";
//    simulations.Push(H1Case_2);
//    H1Case_2.elemen_type = 1;
//    H1Case_2.dump_folder = "H1_sphere_P";
//    simulations.Push(H1Case_2);
//    H1Case_2.elemen_type = 2;
//    H1Case_2.dump_folder = "H1_sphere_H";
//    simulations.Push(H1Case_2);
    
    
//    // Dual Formulation over the solid sphere
//    struct SimulationCase HdivCase_1 = common;
//    HdivCase_1.IsHdivQ = true;
//    HdivCase_1.mesh_type = "linear";
//    HdivCase_1.elemen_type = 2;
//    HdivCase_1.dump_folder = "Hdiv_sphere";
//    simulations.Push(HdivCase_1);
    
    
    // Dual Formulation over the solid sphere
    struct SimulationCase HdivCase_2 = common;
    HdivCase_2.IsHdivQ = true;
    HdivCase_2.mesh_type = "quadratic";
//    HdivCase_2.elemen_type = 0;
//    HdivCase_2.dump_folder = "Hdiv_sphere_T";
//    simulations.Push(HdivCase_2);
//    HdivCase_2.elemen_type = 1;
//    HdivCase_2.dump_folder = "Hdiv_sphere_P";
//    simulations.Push(HdivCase_2);
    HdivCase_2.elemen_type = 2;
    HdivCase_2.dump_folder = "Hdiv_sphere_H";
    simulations.Push(HdivCase_2);
    

//    // Dual Formulation over the solid sphere
//    struct SimulationCase HdivplusCase_1 = common;
//    HdivplusCase_1.IsHdivQ = true;
//    HdivplusCase_1.n_acc_terms = 1;
//    HdivplusCase_1.mesh_type = "linear";
//    HdivplusCase_1.elemen_type = 2;
//    HdivplusCase_1.dump_folder = "Hdivplus_sphere";
//    simulations.Push(HdivplusCase_1);
    
    
//    // Dual Formulation over the solid sphere
//    struct SimulationCase HdivplusCase_2 = common;
//    HdivplusCase_2.IsHdivQ = true;
//    HdivplusCase_2.n_acc_terms = 1;
//    HdivplusCase_2.mesh_type = "quadratic";
//    HdivplusCase_2.elemen_type = 0;
//    HdivplusCase_2.dump_folder = "Hdivplus_sphere_T";
//    simulations.Push(HdivplusCase_2);
//    HdivplusCase_2.elemen_type = 1;
//    HdivplusCase_2.dump_folder = "Hdivplus_sphere_P";
//    simulations.Push(HdivplusCase_2);
//    HdivplusCase_2.elemen_type = 2;
//    HdivplusCase_2.dump_folder = "Hdivplus_sphere_H";
//    simulations.Push(HdivplusCase_2);
    
    
// Cylinder Darcy formulation for Thiem solution in a vertical well
    
//    // Primal Formulation over the solid cylinder
//    struct SimulationCase H1Case_1_cyl = common;
//    H1Case_1_cyl.IsHdivQ = false;
//    H1Case_1_cyl.domain_type = "cylinder";
//    H1Case_1_cyl.mesh_type = "linear";
//    H1Case_1_cyl.dump_folder = "H1_cylinder";
//    simulations.Push(H1Case_1_cyl);
    
//    // Dual Formulation over the solid cylinder
//    struct SimulationCase HdivCase_1_cyl = common;
//    HdivCase_1_cyl.IsHdivQ = true;
//    HdivCase_1_cyl.domain_type = "cylinder";
//    HdivCase_1_cyl.mesh_type = "linear";
//    HdivCase_1_cyl.dump_folder = "Hdiv_cylinder";
//    simulations.Push(HdivCase_1_cyl);
    
    
//    // MHM Dual Formulation over the solid cylinder
//    struct SimulationCase HdivMHMCase_1_cyl = common;
//    HdivMHMCase_1_cyl.IsHdivQ = true;
//    HdivMHMCase_1_cyl.IsMHMQ = true;
//    HdivMHMCase_1_cyl.domain_type = "cylinder";
//    HdivMHMCase_1_cyl.mesh_type = "linear";
//    HdivMHMCase_1_cyl.dump_folder = "HdivMHM_cylinder";
//    simulations.Push(HdivMHMCase_1_cyl);
    
    ComputeCases(simulations);
    
    return 0;
}

void ComputeCases(TPZStack<SimulationCase> cases){
    

    
#ifdef USING_BOOST
    boost::posix_time::ptime int_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    int n_cases = cases.size();
    for (int i = 0; i < n_cases; i++) {
        ComputeApproximation(cases[i]);
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
    std::cout << "End:: Overal time = " << int_t2-int_t1 << std::endl;
    
}

void ComputeApproximation(SimulationCase & sim_data){
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Creating the directory
    std::string command = "mkdir " + sim_data.dump_folder;
    system(command.c_str());
    
    std::stringstream summary;
    summary   << sim_data.dump_folder << "/" "conv" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << ".txt";
    std::ofstream convergence(summary.str(),std::ios::app);
    
    TPZManVector<STATE,10> p_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<STATE,10> d_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<STATE,10> h_error(sim_data.n_h_levels+1,1.0);
    
    TPZManVector<STATE,10> p_conv(sim_data.n_h_levels,0.0);
    TPZManVector<STATE,10> d_conv(sim_data.n_h_levels,0.0);
    TPZManVector<STATE,10> h_conv(sim_data.n_h_levels,0.0);
    
    int n_h_levels = sim_data.n_h_levels;
    int n_p_levels = sim_data.n_p_levels;
    

    using namespace std;
    
    for (int p = 1; p <= n_p_levels; p++) {
        
        convergence << std::endl;        
        convergence << " Polynomial order  =  " << p << std::endl;
        convergence << setw(5)  << " h" << setw(25) << " ndof" << setw(25) << " ndof_cond" << setw(25) << " assemble_time (msec)" << setw(25) << " solving_time (msec)" << setw(25) << " error_time (msec)" << setw(25) << " Primal l2 error" << setw(25) << " Dual l2 error"  << setw(25) << " H error (H1 or Hdiv)" << endl;
        
        int h_base = 0;
        if (sim_data.IsMHMQ) {
            h_base = n_h_levels;
        }
        for (int h = 0; h <= n_h_levels; h++) {
            
            // Compute the geometry
            TPZGeoMesh * gmesh = GeomtricMesh(h_base, sim_data);
    
#ifdef PZDEBUG
            TPZCheckGeom check(gmesh);
            int checkQ = check.PerformCheck();
            if (checkQ) {
                DebugStop();
            }
#endif
            if (!sim_data.IsMHMQ) {
                if (sim_data.elemen_type == 0) {
                    UniformRefineTetrahedrons(gmesh, h);
                }
                else{
                    UniformRefinement(gmesh, h);
                }
                
            }
            else{
                level_mhm = h;
            }


//#ifdef PZDEBUG
//            
//#ifdef USING_BOOST
//            boost::posix_time::ptime int_t1 = boost::posix_time::microsec_clock::local_time();
//#endif
//            
//            STATE volume = IntegrateVolume(gmesh, sim_data);
//            
//#ifdef USING_BOOST
//            boost::posix_time::ptime int_t2 = boost::posix_time::microsec_clock::local_time();
//#endif
//            
//            std::cout << "Domain volume = " << volume << "; Time for integration = " << int_t2-int_t1 <<std::endl;
//            
//#endif
            std::stringstream text_name;
            std::stringstream vtk_name;
            text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".txt";
            vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
//            ofstream textfile(text_name.str());
//            gmesh->Print(textfile);
            
//            std::ofstream vtkfile(vtk_name.str());
//            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

            
            // Compute the geometry
            long ndof, ndof_cond;
            TPZManVector<TPZCompMesh *,5> meshvec;
            TPZCompMesh * cmesh = ComputationalMesh(gmesh, p, sim_data, ndof, meshvec);

            // Create Analysis
            TPZAnalysis * analysis = CreateAnalysis(cmesh,sim_data);
            
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            analysis->Assemble();
            
            ndof_cond = analysis->Rhs().Rows();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            
            analysis->Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
            STATE assemble_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
            STATE solving_time  = boost::numeric_cast<double>((t3-t2).total_milliseconds());
#endif
            
#ifdef USING_BOOST
                boost::posix_time::ptime int_unwrap_t1 = boost::posix_time::microsec_clock::local_time();
#endif
                UnwrapMesh(cmesh);
                analysis->LoadSolution();
                cmesh->Solution() *= -1.0; /* consequence of newton correction */
                analysis->LoadSolution(cmesh->Solution());
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh);
                
#ifdef USING_BOOST
                boost::posix_time::ptime int_unwrap_t2 = boost::posix_time::microsec_clock::local_time();
#endif
                
            std::cout << "StaticCondensation::Time for uncondense equations = " << int_unwrap_t2-int_unwrap_t1 <<std::endl;

            
            // PostProccessing
            std::stringstream sol_vtk_name;
            sol_vtk_name    << sim_data.dump_folder << "/" "sol" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
//            std::string file(sol_vtk_name.str());
//            PosProcess(analysis, file, sim_data);
            
            
            // compute the error
#ifdef USING_BOOST
            boost::posix_time::ptime error_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            if (sim_data.IsHdivQ) {
                
                bool par_errorQ = true;
                if(par_errorQ){
                    int nthreadsForError = sim_data.n_threads;
                    analysis->SetExact(Analytic);
                    TPZManVector<REAL,3> errors(3,1);
                    analysis->SetThreadsForError(nthreadsForError);
                    analysis->PostProcessError(errors);
                    p_error[h] = errors[0];
                    d_error[h] = errors[1];
                    h_error[h] = errors[2];
                }
                else{
                    ErrorHdiv(cmesh, p_error[h], d_error[h], h_error[h]);
                }
            }
            else{
                ErrorH1(cmesh, p_error[h], d_error[h], h_error[h]);
            }
            
#ifdef USING_BOOST
            boost::posix_time::ptime error_t2 = boost::posix_time::microsec_clock::local_time();
#endif
            STATE error_time = boost::numeric_cast<double>((error_t2 - error_t1).total_milliseconds());
            
            // current summary
            convergence << setw(5) << h << setw(25) << ndof << setw(25) << ndof_cond << setw(25) << assemble_time << setw(25) << solving_time << setw(25) << error_time << setw(25) << p_error[h] << setw(25) << d_error[h]  << setw(25) << h_error[h] << endl;

            // Release memory
            analysis->Solver().Matrix()->Zero();
            delete cmesh;
            for (int i = 0; i < meshvec.size(); i++) {
                delete meshvec[i];
            }
            analysis->CleanUp();
            delete gmesh;
            
        }
        
        // compute rates
        ComputeConvergenceRates(p_error,p_conv);
        ComputeConvergenceRates(d_error,d_conv);
        ComputeConvergenceRates(h_error,h_conv);
        
        
        // print convergence summary
        convergence << std::endl;
        convergence << " Convergence rates summary " << std::endl;
        convergence << " Polynomial order  =  " << p << std::endl;
        convergence << " Primal convergence rates = " << setw(5) << p_conv << std::endl;
        convergence << " Dual convergence rates = " << setw(5) << d_conv << std::endl;
        convergence << " H1 or Hdiv convergence rates = " << setw(5) << h_conv << std::endl;
        convergence << std::endl;
        convergence << " ------------------------------------------------------------------ " << std::endl;
        convergence.flush();
        
    }
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
#ifdef USING_BOOST
    STATE case_solving_time = boost::numeric_cast<double>((int_case_t2-int_case_t1).total_milliseconds());
#endif
    
    convergence <<  "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
    
}

void ComputeConvergenceRates(TPZVec<STATE> &error, TPZVec<STATE> &convergence){
    
    int ndata = error.size();
    STATE log0p5 = log(0.5);
    for (int i = 1; i < ndata; i++) {
        STATE logerror = log(error[i-1]);
        STATE logerrori = log(error[i]);
        convergence[i-1] = (logerrori - logerror)/log0p5;
    }
}


void Analytic(const TPZVec<REAL> &p, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu){
    
    gradu.Resize(4,1);
    
    STATE x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
    STATE r = sqrt(x*x+y*y+z*z);
    STATE theta = acos(z/r);
    STATE phi = atan2(y,x);
    
    STATE costheta = cos(theta);
    STATE cos2theta = cos(2.0*theta);
    STATE cosphi = cos(phi);
    STATE cos2phi = cos(2.0*phi);
    STATE sintheta = sin(theta);
    STATE sin2theta = sin(2.0*theta);
    STATE sinphi = sin(phi);
    STATE sin2phi = sin(2.0*phi);
    
    // Gradient computations
    
    STATE Radialunitx = sintheta*cosphi;
    STATE Radialunity = sintheta*sinphi;
    STATE Radialunitz = costheta;
    
    STATE Thetaunitx = cosphi*costheta;
    STATE Thetaunity = costheta*sinphi;
    STATE Thetaunitz = -sintheta;
    
    STATE Phiunitx = -sinphi;
    STATE Phiunity = cosphi;
    STATE Phiunitz = 0.0;
    
#ifdef Solution1
    
    STATE dfdr       = 2.0*r;
    STATE dfdTheta   = 0.0;
    STATE dfdPhi     = 0.0;
    
    u[0] = r*r;
    
    gradu(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    gradu(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
    gradu(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    
    gradu(3,0) = -6.0;
    
#endif
    
#ifdef Solution2
    
    STATE dfdr       = 4.0*r*r*r*sintheta*sintheta*sinphi*sinphi;
    STATE dfdTheta   = r*r*r*sin2theta*sinphi*sinphi;
    STATE dfdPhi     = r*r*r*sintheta*sin2phi;
    
    u[0] = r*r*r*r*sintheta*sintheta*(1.0-cosphi*cosphi);
    
    gradu(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    gradu(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
    gradu(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    
    gradu(3,0) = -r*r*(2.0*cosphi*cosphi+(9.0-7.0*cos2theta)*sinphi*sinphi);
    
#endif
    
#ifdef Solution3
    
    STATE dfdr       = -4.0*M_PI*r*sin(2.0*M_PI*r*r);
    STATE dfdTheta   = 0.0;
    STATE dfdPhi     = 0.0;
    
    u[0] = cos(2.0*M_PI*r*r);
    
    gradu(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx + dfdPhi * Phiunitx);
    gradu(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity + dfdPhi * Phiunity);
    gradu(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz + dfdPhi * Phiunitz);
    
    gradu(3,0) = 4.0*M_PI*(4.0*M_PI*r*r*cos(2.0*M_PI*r*r) + 3.0*sin(2.0*M_PI*r*r));
    
#endif
    
#ifdef Solution6
    
    REAL a = +5.0/4.0;
    REAL b = -1.0/4.0;
    REAL c = -1.0/4.0;
    
    REAL d = 5.0;
    
    REAL xma = x-a;
    REAL ymb = y-b;
    REAL zmc = z-c;
    REAL piover2 = M_PI/2.0;
    REAL piover3 = M_PI/3.0;
    REAL rad = xma*xma + ymb*ymb+ zmc*zmc;
    REAL sqrt_rad = sqrt(rad);
    REAL artan_arg = atan(d*(sqrt_rad - piover3));
    
    
    REAL denomfactor1 = -9.0 + d*d*(-9.0*rad+M_PI*(-M_PI+6.0*sqrt_rad));
    REAL numfactro1   = 18.0*d*(-9.0+d*d*M_PI*(-M_PI+3.0*sqrt_rad));
    
    u[0] = piover2 - artan_arg;
    
    gradu(0,0) = -1.0*((-d*xma)/( (1.0+d*d*(sqrt_rad - piover3)*(sqrt_rad - piover3)) * (sqrt_rad) ));
    gradu(1,0) = -1.0*((-d*ymb)/( (1.0+d*d*(sqrt_rad - piover3)*(sqrt_rad - piover3)) * (sqrt_rad) ));
    gradu(2,0) = -1.0*((-d*zmc)/( (1.0+d*d*(sqrt_rad - piover3)*(sqrt_rad - piover3)) * (sqrt_rad) ));
    
    gradu(3,0) = -1.0*(numfactro1)/(denomfactor1*denomfactor1*sqrt_rad);
    
#endif
    
#ifdef Thiem
    
    gradu.Resize(4,1);
   
    REAL r0 = 100.0;
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    
    costheta = cos(theta);
    sintheta = sin(theta);
    
    // Gradient computations
    
    Radialunitx = costheta;
    Radialunity = sintheta;
    Radialunitz = 0.0;
    
    Thetaunitx = -sintheta;
    Thetaunity = costheta;
    Thetaunitz = 0.0;
    
    u[0] = log(r/r0);
    
    REAL dfdr = 1.0/r;
    REAL dfdTheta = 0.0;
    
    gradu(0,0) = -1.0*(dfdr * Radialunitx + dfdTheta * Thetaunitx);
    gradu(1,0) = -1.0*(dfdr * Radialunity + dfdTheta * Thetaunity);
    gradu(2,0) = -1.0*(dfdr * Radialunitz + dfdTheta * Thetaunitz);
    
    gradu(3,0) = 0.0;//-4.0;
    
#endif
    
}

void Solution(const TPZVec<REAL> &p, TPZVec<STATE> &f){

    REAL x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
#ifdef Solution1
    
    f[0] = r*r;
    
#endif
   
#ifdef Solution2
    
    STATE cosphi = cos(phi);;
    STATE sintheta = sin(theta);
    
    f[0] = r*r*r*r*sintheta*sintheta*(1.0-cosphi*cosphi);

#endif
    
    
#ifdef Solution3
    
    f[0] = cos(2.0*M_PI*r*r);
    
#endif
    
#ifdef Solution6
    
    REAL a = +5.0/4.0;
    REAL b = -1.0/4.0;
    REAL c = -1.0/4.0;
    
    REAL d = 5.0;
    
    REAL xma = x-a;
    REAL ymb = y-b;
    REAL zmc = z-c;
    REAL piover2 = M_PI/2.0;
    REAL piover3 = M_PI/3.0;
    REAL rad = xma*xma + ymb*ymb+ zmc*zmc;
    REAL sqrt_rad = sqrt(rad);
    REAL artan_arg = atan(d*(sqrt_rad - piover3));
    
    
    REAL denomfactor1 = -9.0 + d*d*(-9.0*rad+M_PI*(-M_PI+6.0*sqrt_rad));
    REAL numfactro1   = 18.0*d*(-9.0+d*d*M_PI*(-M_PI+3.0*sqrt_rad));
    
    f[0] = piover2 - artan_arg;
    
#endif
    
#ifdef Thiem
    
    REAL r0 = 100.0;
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    
    f[0] = log(r/r0);
    
    
#endif
    
    
}

void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
    
    REAL x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
    REAL r = sqrt(x*x+y*y+z*z);
    REAL theta = acos(z/r);
    REAL phi = atan2(y,x);
    
#ifdef Solution1
    
    f[0] = -6.0;
    
#endif
    
#ifdef Solution2
    
    STATE cos2theta = cos(2.0*theta);
    STATE cosphi = cos(phi);
    STATE sinphi = sin(phi);

    f[0] = -r*r*(2.0*cosphi*cosphi+(9.0-7.0*cos2theta)*sinphi*sinphi);
    
#endif
    
#ifdef Solution3
    
    f[0] = 4.0*M_PI*(4.0*M_PI*r*r*cos(2.0*M_PI*r*r) + 3.0*sin(2.0*M_PI*r*r));
    
#endif
    
    
#ifdef Solution6
    
    REAL a = +5.0/4.0;
    REAL b = -1.0/4.0;
    REAL c = -1.0/4.0;
    
    REAL d = 5.0;
    
    REAL xma = x-a;
    REAL ymb = y-b;
    REAL zmc = z-c;
    REAL piover2 = M_PI/2.0;
    REAL piover3 = M_PI/3.0;
    REAL rad = xma*xma + ymb*ymb+ zmc*zmc;
    REAL sqrt_rad = sqrt(rad);
    REAL artan_arg = atan(d*(sqrt_rad - piover3));
    
    
    REAL denomfactor1 = -9.0 + d*d*(-9.0*rad+M_PI*(-M_PI+6.0*sqrt_rad));
    REAL numfactro1   = 18.0*d*(-9.0+d*d*M_PI*(-M_PI+3.0*sqrt_rad));
    
    f[0] = -1.0*(numfactro1)/(denomfactor1*denomfactor1*sqrt_rad);
    
#endif
    
#ifdef Thiem
 
    f[0] = 0.0;//-4.0;
    
#endif
    
}

STATE IntegrateVolume(TPZGeoMesh * geometry, SimulationCase sim_data){

    int order = sim_data.int_order;
    int nel = geometry->NElements();
    STATE domain_volume = 0.0;
    
    for(int iel  = 0; iel < nel; iel++)
    {
    
        TPZGeoEl * gel = geometry->Element(iel);

#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() !=3) {
            continue;
        }
        
        if (gel->HasSubElement()){
            continue;
        }
        
        int gel_volume_side = gel->NSides() - 1;
        TPZIntPoints * int_rule = gel->CreateSideIntegrationRule(gel_volume_side, order);
        int npoints = int_rule->NPoints();
        
        TPZManVector<STATE,3> triplet(3,0.0);
        STATE w;
        
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        REAL detjac;
        
        STATE el_volume = 0.0;
        for (int i = 0; i < npoints ; i++) {
            int_rule->Point(i, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            el_volume += w * detjac;
        }
        
        domain_volume += el_volume;
    }
    return domain_volume;
}

STATE IntegrateSolution(TPZCompMesh * cmesh, SimulationCase sim_data){
    
    int order = sim_data.int_order;
    STATE p_integral = 0.0;
    
    long nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerror(3,0.   );
    for (long iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        
#ifdef PZDEBUG
        if (!cel) {
            DebugStop();
        }
#endif
        
        if(cel->Reference()->Dimension()!=dim) {
            continue;
        }
        
        TPZGeoEl * gel = cel->Reference();
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        int gel_volume_side = gel->NSides() - 1;
        TPZIntPoints * int_rule = gel->CreateSideIntegrationRule(gel_volume_side, order);
        int npoints = int_rule->NPoints();
        
        TPZManVector<STATE,3> triplet(3,0.0);
        STATE w;
        
        TPZFMatrix<REAL> jac;
        TPZFMatrix<REAL> axes;
        TPZFMatrix<REAL> jacinv;
        REAL detjac;
        TPZManVector<STATE,1> p;
        STATE el_interal = 0.0;
        for (int i = 0; i < npoints ; i++) {
            int_rule->Point(i, triplet, w);
            gel->Jacobian(triplet, jac, axes, detjac, jacinv);
            cel->Solution(triplet, 4, p);
            el_interal += w * detjac * p[0];
        }
        
        p_integral += el_interal;
    }
    
    return p_integral;
    
}

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase  & sim_data){
    
    TPZGeoMesh * geometry = NULL;
    
    if (sim_data.domain_type == "sphere") {
        
        if (sim_data.mesh_type == "linear") {
            geometry = MakeSphereFromLinearQuadrilateralFaces(ndiv, sim_data);
            return geometry;
        }
        
        if (sim_data.mesh_type == "quadratic") {
            geometry = MakeSphereFromQuadrilateralFacesR(ndiv, sim_data);
            if(!sim_data.UseGmshMeshQ){
                TransformToQuadratic(geometry);
            }
            return geometry;
        }
        
        if (sim_data.mesh_type == "blended") {
            geometry = MakeSphereFromQuadrilateralFacesR(ndiv, sim_data);
            return geometry;
        }
        
        std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
        DebugStop();
        
    }
    
    if (sim_data.domain_type == "cylinder") {
        
        if (sim_data.mesh_type == "linear") {
            geometry = MakeCylinderFromLinearFaces(ndiv, sim_data);
            return geometry;
        }
        
        std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
        DebugStop();
        
    }
    
    std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
    DebugStop();
    return geometry;
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof, TPZVec<TPZCompMesh *> &meshvec){
    
    TPZCompMesh * mesh = NULL;
    
    if (sim_data.IsHdivQ) {
        mesh = DualMesh(geometry, p, sim_data, meshvec);
        ndof = meshvec[0]->NEquations() + meshvec[1]->NEquations();
        return mesh;
    }
    else
    {
        mesh = PrimalMesh(geometry, p, sim_data, ndof);
        return mesh;        
    }
    
    std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
    DebugStop();
    return mesh;
}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);

    if (sim_data.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix matrix(cmesh);
//        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    
    if (sim_data.UseFrontalQ) {

        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        matrix.SetDecomposeType(ELDLt);
        matrix.SetNumThreads(sim_data.n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);        
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    

    
   return analysis;
    
}

void PosProcess(TPZAnalysis  * an, std::string file, SimulationCase & sim_data)
{
    int dim = 3;
    TPZStack<std::string,10> scalnames, vecnames;
    
    if (sim_data.IsHdivQ) {
        vecnames.Push("q");
        vecnames.Push("q_exact");
        scalnames.Push("p");
        scalnames.Push("p_exact");
        scalnames.Push("f_exact");
        scalnames.Push("div_q");
    }
    else{
        vecnames.Push("q");
        vecnames.Push("q_exact");
        scalnames.Push("p");
        scalnames.Push("p_exact");
        scalnames.Push("f_exact");
    }

    int div = 1;
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
    
}

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, long &ndof){
    
    int dimension = 3;
    int dirichlet = 0;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
    std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
        
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        TPZMaterial * volume = new TPZPrimalPoisson(sim_data.omega_ids[iv]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution);
            analytic_bc->SetPolynomialOrder(sim_data.int_order);
            TPZAutoPointer< TPZFunction<STATE> > solution = analytic_bc;
            
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            face->SetBCForcingFunction(solution);
            cmesh->InsertMaterialObject(face);
        }
        
    }
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    ndof = cmesh->NEquations();
    
    TPZCompMeshTools::GroupElements(cmesh);
    TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
    
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Primal_cmesh" << ".txt";
    std::ofstream sout(file_name.str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}


TPZCompMesh *DualMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec)
{
    
    int dimension = 3;
    int dirichlet = 0;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();
    
    if (sim_data.IsMHMQ) {
        InsertSkeletonInterfaces(geometry);
        PrintGeometry(geometry, sim_data);
    }
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
        std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        TPZMaterial * volume = new TPZDualPoisson(sim_data.omega_ids[iv]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            
            TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution);
            analytic_bc->SetPolynomialOrder(sim_data.int_order);
            TPZAutoPointer< TPZFunction<STATE> > solution = analytic_bc;
            
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            face->SetBCForcingFunction(solution);
            cmesh->InsertMaterialObject(face);
        }
        
    }
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = qMesh(geometry, p, sim_data);
    meshvector[1] = pMesh(geometry, p, sim_data);
    
    AdjustFluxPolynomialOrders(meshvector[0], sim_data.n_acc_terms);
    SetPressureOrders(meshvector[0], meshvector[1]);
    
    if (sim_data.IsMHMQ) {
        SeparateConnectsByNeighborhood(meshvector[0]);
    }
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    
    if (sim_data.IsMHMQ) {
//        BuildMacroElements(cmesh);
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    else{
        
        TPZCompMeshTools::GroupElements(cmesh);
        TPZCompMeshTools::CreatedCondensedElements(cmesh, true);
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str());
    cmesh->Print(sout);
#endif
    
    meshvec = meshvector;
    
    return cmesh;
    

}

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus)
{
    int dim = fluxmesh->Dimension();
    /// loop over all the elements
    long nel = fluxmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        // compute the maxorder
        int maxorder = -1;
        int ncon = intel->NConnects();
        for (int i=0; i<ncon-1; i++) {
            int conorder = intel->Connect(i).Order();
            maxorder = maxorder < conorder ? conorder : maxorder;
        }
        int nsides = gel->NSides();
        int nconside = intel->NSideConnects(nsides-1);
        // tive que tirar para rodar H1
        //        if (nconside != 1 || maxorder == -1) {
        //            DebugStop();
        //        }
        long cindex = intel->SideConnectIndex(nconside-1, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        if (c.NElConnected() != 1) {
            DebugStop();
        }
        if (c.Order()+hdivplusplus != maxorder) {
            //            std::cout << "Changing the order of the central connect " << cindex << " from " << c.Order() << " to " << maxorder+hdivplusplus << std::endl;
            // change the internal connect order to be equal do maxorder
            intel->SetSideOrder(nsides-1, maxorder+hdivplusplus);
        }
    }
    fluxmesh->ExpandSolution();
}

void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh)
{
    // build a vector with the required order of each element in the pressuremesh
    // if an element of the mesh dimension of the fluxmesh does not have a corresponding element in the pressuremesh DebugStop is called
    int meshdim = fluxmesh->Dimension();
    pressuremesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    TPZManVector<long> pressorder(pressuremesh->NElements(),-1);
    long nel = fluxmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (gel->Dimension() != meshdim) {
            continue;
        }
        int nsides = gel->NSides();
        long cindex = intel->SideConnectIndex(0, nsides-1);
        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
        int order = c.Order();
        TPZCompEl *pressureel = gel->Reference();
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(pressureel);
        if (!pintel) {
            DebugStop();
        }
        pressorder[pintel->Index()] = order;
    }
    pressuremesh->Reference()->ResetReference();
    nel = pressorder.size();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        TPZInterpolatedElement *pintel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!pintel) {
            continue;
        }
        if (pressorder[el] == -1) {
            continue;
        }
        pintel->PRefine(pressorder[el]);
    }
    
    pressuremesh->ExpandSolution();
}

TPZCompMesh * qMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data){
    
    int dimension = 3;
    int dirichlet = 0;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
        std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    int Skeleton_material_Id = 100;
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    std::set<int> set_vol, set_skeleton;
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        TPZMaterial * volume = new TPZDualPoisson(sim_data.omega_ids[iv]);
        set_vol.insert(sim_data.omega_ids[iv]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution);
        analytic_bc->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > solution = analytic_bc;
        volume->SetBCForcingFunction(solution);
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            set_vol.insert(sim_data.gamma_ids[ib]);
            cmesh->InsertMaterialObject(face);
        }
        
    }
    

    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild(set_vol);
    
    if (sim_data.IsMHMQ) {
        TPZDualPoisson * mat_skeleton = new TPZDualPoisson(Skeleton_material_Id);
        set_skeleton.insert(Skeleton_material_Id);
        cmesh->InsertMaterialObject(mat_skeleton); // @omar::  skeleton material inserted
        cmesh->SetDefaultOrder(p);
        if (p > 1) {
//            cmesh->SetDefaultOrder(p-1);
            cmesh->SetDefaultOrder(p);
        }
        
        cmesh->AutoBuild(set_skeleton);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "q_cmesh" << ".txt";
    std::ofstream sout(file_name.str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
TPZCompMesh * pMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data){
    
    int dimension = 3;
    int nvolumes = sim_data.omega_ids.size();
    
#ifdef PZDEBUG
    if (nvolumes != 1) {
        std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
        DebugStop();
    }
#endif
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    for (int iv = 0; iv < nvolumes ; iv++) {
        
        TPZMaterial * volume = new TPZMatPoisson3d(sim_data.omega_ids[iv]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution);
        analytic_bc->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > solution = analytic_bc;
        volume->SetBCForcingFunction(solution);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(p);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild();
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "p_cmesh" << ".txt";
    std::ofstream sout(file_name.str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}

TPZGeoMesh * MakeSphereFromLinearQuadrilateralFaces(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Sphere:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    if (sim_data.UseGmshMeshQ) {
        
        std::string dirname = PZSOURCEDIR;
        std::string grid;
        
        switch (sim_data.elemen_type) {
            case 0:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryT.msh";
            }
                break;
            case 1:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryP.msh";
            }
                break;
            case 2:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryH.msh";
            }
                break;
            default:
            {
                DebugStop();
            }
                break;
        }

        
        TPZGmshReader Geometry;
        REAL s = 1.0;
        Geometry.SetfDimensionlessL(s);
        geomesh = Geometry.GeometricGmshMesh(grid);
        const std::string name("Sphere from gmsh script");
        geomesh->SetName(name);
        
        // changin id internally
        sim_data.gamma_ids[0] = 2;
        sim_data.gamma_ids[1] = 3;
    }
    else
    {
        geomesh->SetDimension(3);
        int nl = 2;// Let it fixed
        int basenodes = 8;
        int nodes =  basenodes * (nl);
        REAL radius_o = 1.0;
        REAL radius_i = 0.25;
        
        REAL dr = (radius_o- radius_i)/REAL(nl-1);
        geomesh->SetMaxNodeId(nodes-1);
        geomesh->NodeVec().Resize(nodes);
        TPZManVector<TPZGeoNode,4> Node(nodes);
        
        TPZManVector<long,4> TopolQuad(4);
        TPZManVector<long,8> TopolCube(8);
        TPZManVector<REAL,3> coord(3,0.);
        TPZVec<REAL> xc(3,0.);
        REAL cphi = atan(sqrt(2.0));
        
        int nodeindex = 0;
        long id = 0;
        int matid = sim_data.omega_ids[0];
        
        TPZManVector< TPZManVector<REAL,3> , 8 > points(basenodes,0.);
        for (int il = 0; il < nl; il++) {
            
            if (il==0) {
                matid = sim_data.gamma_ids[0];
            }
            
            if (il==nl-1) {
                matid = sim_data.gamma_ids[1];
            }
            
            REAL radius = radius_o - REAL(il)*dr;
            points[0].Resize(3, 0.0);
            points[0][0]=radius;
            points[0][1]=cphi;
            points[0][2]=M_PI/4.0;
            
            points[1].Resize(3, 0.0);
            points[1][0]=radius;
            points[1][1]=cphi;
            points[1][2]=-M_PI/4.0;
            
            points[2].Resize(3, 0.0);
            points[2][0]=radius;
            points[2][1]=cphi;
            points[2][2]=-3.0*M_PI/4.0;
            
            points[3].Resize(3, 0.0);
            points[3][0]=radius;
            points[3][1]=cphi;
            points[3][2]=3.0*M_PI/4.0;
            
            points[4].Resize(3, 0.0);
            points[4][0]=radius;
            points[4][1]=M_PI-cphi;
            points[4][2]=M_PI/4.0;
            
            points[5].Resize(3, 0.0);
            points[5][0]=radius;
            points[5][1]=M_PI-cphi;
            points[5][2]=-M_PI/4.0;
            
            points[6].Resize(3, 0.0);
            points[6][0]=radius;
            points[6][1]=M_PI-cphi;
            points[6][2]=-3.0*M_PI/4.0;
            
            points[7].Resize(3, 0.0);
            points[7][0]=radius;
            points[7][1]=M_PI-cphi;
            points[7][2]=3.0*M_PI/4.0;
            
            for (int i = 0; i < basenodes; i++) {
                coord = ParametricSphere(points[i][0],points[i][1],points[i][2]);
                geomesh->NodeVec()[nodeindex].SetCoord(coord);
                geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
                nodeindex++;
            }
            
            if (il==0 || il==nl-1) {
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 1+il*basenodes;
                TopolQuad[2] = 2+il*basenodes;
                TopolQuad[3] = 3+il*basenodes;
                new TPZGeoElRefPattern<  pzgeom::TPZGeoQuad  > (id,TopolQuad,matid,*geomesh);
                id++;
                
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 3+il*basenodes;
                TopolQuad[2] = 7+il*basenodes;
                TopolQuad[3] = 4+il*basenodes;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
                id++;
                
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 1+il*basenodes;
                TopolQuad[2] = 5+il*basenodes;
                TopolQuad[3] = 4+il*basenodes;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
                id++;
                
                TopolQuad[0] = 3+il*basenodes;
                TopolQuad[1] = 2+il*basenodes;
                TopolQuad[2] = 6+il*basenodes;
                TopolQuad[3] = 7+il*basenodes;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad  > (id,TopolQuad,matid,*geomesh);
                id++;
                
                TopolQuad[0] = 1+il*basenodes;
                TopolQuad[1] = 2+il*basenodes;
                TopolQuad[2] = 6+il*basenodes;
                TopolQuad[3] = 5+il*basenodes;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
                id++;
                
                TopolQuad[0] = 4+il*basenodes;
                TopolQuad[1] = 5+il*basenodes;
                TopolQuad[2] = 6+il*basenodes;
                TopolQuad[3] = 7+il*basenodes;
                new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*geomesh);
                id++;
            }

        }
        
        matid = sim_data.omega_ids[0];
        
        for (int il = 0; il < nl - 1 ; il++) {
            //      Inserting blend elements
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 1+il*basenodes;
            TopolCube[2] = 2+il*basenodes;
            TopolCube[3] = 3+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 1+(il+1)*basenodes;
            TopolCube[6] = 2+(il+1)*basenodes;
            TopolCube[7] = 3+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 3+il*basenodes;
            TopolCube[2] = 7+il*basenodes;
            TopolCube[3] = 4+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 3+(il+1)*basenodes;
            TopolCube[6] = 7+(il+1)*basenodes;
            TopolCube[7] = 4+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 1+il*basenodes;
            TopolCube[2] = 5+il*basenodes;
            TopolCube[3] = 4+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 1+(il+1)*basenodes;
            TopolCube[6] = 5+(il+1)*basenodes;
            TopolCube[7] = 4+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
            TopolCube[0] = 3+il*basenodes;
            TopolCube[1] = 2+il*basenodes;
            TopolCube[2] = 6+il*basenodes;
            TopolCube[3] = 7+il*basenodes;
            TopolCube[4] = 3+(il+1)*basenodes;
            TopolCube[5] = 2+(il+1)*basenodes;
            TopolCube[6] = 6+(il+1)*basenodes;
            TopolCube[7] = 7+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
            TopolCube[0] = 1+il*basenodes;
            TopolCube[1] = 2+il*basenodes;
            TopolCube[2] = 6+il*basenodes;
            TopolCube[3] = 5+il*basenodes;
            TopolCube[4] = 1+(il+1)*basenodes;
            TopolCube[5] = 2+(il+1)*basenodes;
            TopolCube[6] = 6+(il+1)*basenodes;
            TopolCube[7] = 5+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
            TopolCube[0] = 4+il*basenodes;
            TopolCube[1] = 5+il*basenodes;
            TopolCube[2] = 6+il*basenodes;
            TopolCube[3] = 7+il*basenodes;
            TopolCube[4] = 4+(il+1)*basenodes;
            TopolCube[5] = 5+(il+1)*basenodes;
            TopolCube[6] = 6+(il+1)*basenodes;
            TopolCube[7] = 7+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (id,TopolCube,matid,*geomesh);
            id++;
            
        }
    }
    geomesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    int axis = 1;
    REAL angle = 0.0;//-45.0;
    RotateGeomesh(geomesh, angle, axis);
    return geomesh;
}

TPZGeoMesh * MakeSphereFromQuadrilateralFaces(int ndiv, SimulationCase  & sim_data)
{
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Sphere:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(3);
    int nl = 2;// Let it fixed
    int basenodes = 8;
    int nodes =  basenodes * (nl);
    REAL radius_o = 1.0;
    REAL radius_i = 0.25;
    
    REAL dr = (radius_o - radius_i)/REAL(nl-1);
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<long,4> TopolQuad(4);
    TPZManVector<long,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    long id = 0;
    int matid = sim_data.omega_ids[0];
    
    TPZManVector< TPZManVector<REAL,3> , 8 > points(basenodes,0.);
    for (int il = 0; il < nl; il++) {
        
        if (il==0) {
            matid = sim_data.gamma_ids[0];
        }
        
        if (il==nl-1) {
            matid = sim_data.gamma_ids[1];
        }
        
        REAL radius = radius_o - REAL(il)*dr;
        points[0].Resize(3, 0.0);
        points[0][0]=radius;
        points[0][1]=cphi;
        points[0][2]=M_PI/4.0;
        
        points[1].Resize(3, 0.0);
        points[1][0]=radius;
        points[1][1]=cphi;
        points[1][2]=-M_PI/4.0;
        
        points[2].Resize(3, 0.0);
        points[2][0]=radius;
        points[2][1]=cphi;
        points[2][2]=-3.0*M_PI/4.0;
        
        points[3].Resize(3, 0.0);
        points[3][0]=radius;
        points[3][1]=cphi;
        points[3][2]=3.0*M_PI/4.0;
        
        points[4].Resize(3, 0.0);
        points[4][0]=radius;
        points[4][1]=M_PI-cphi;
        points[4][2]=M_PI/4.0;
        
        points[5].Resize(3, 0.0);
        points[5][0]=radius;
        points[5][1]=M_PI-cphi;
        points[5][2]=-M_PI/4.0;
        
        points[6].Resize(3, 0.0);
        points[6][0]=radius;
        points[6][1]=M_PI-cphi;
        points[6][2]=-3.0*M_PI/4.0;
        
        points[7].Resize(3, 0.0);
        points[7][0]=radius;
        points[7][1]=M_PI-cphi;
        points[7][2]=3.0*M_PI/4.0;
        
        
        for (int i = 0; i < basenodes; i++) {
            coord = ParametricSphere(points[i][0],points[i][1],points[i][2]);
            geomesh->NodeVec()[nodeindex].SetCoord(coord);
            geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
        
        if (il==0 || il==nl-1)
        {
            TopolQuad[0] = 0+il*basenodes;
            TopolQuad[1] = 1+il*basenodes;
            TopolQuad[2] = 2+il*basenodes;
            TopolQuad[3] = 3+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad1->Geom().SetData(radius, xc);
            id++;
            
            TopolQuad[0] = 0+il*basenodes;
            TopolQuad[1] = 3+il*basenodes;
            TopolQuad[2] = 7+il*basenodes;
            TopolQuad[3] = 4+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad2->Geom().SetData(radius, xc);
            id++;
            
            TopolQuad[0] = 0+il*basenodes;
            TopolQuad[1] = 1+il*basenodes;
            TopolQuad[2] = 5+il*basenodes;
            TopolQuad[3] = 4+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad3->Geom().SetData(radius, xc);
            id++;
            
            TopolQuad[0] = 3+il*basenodes;
            TopolQuad[1] = 2+il*basenodes;
            TopolQuad[2] = 6+il*basenodes;
            TopolQuad[3] = 7+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad4->Geom().SetData(radius, xc);
            id++;
            
            TopolQuad[0] = 1+il*basenodes;
            TopolQuad[1] = 2+il*basenodes;
            TopolQuad[2] = 6+il*basenodes;
            TopolQuad[3] = 5+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad5->Geom().SetData(radius, xc);
            id++;
            
            TopolQuad[0] = 4+il*basenodes;
            TopolQuad[1] = 5+il*basenodes;
            TopolQuad[2] = 6+il*basenodes;
            TopolQuad[3] = 7+il*basenodes;
            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
            quad6->Geom().SetData(radius, xc);
            id++;
        }

    }
    
    matid = sim_data.omega_ids[0];
    
    for (int il = 0; il < nl - 1 ; il++) {
        //      Inserting blend elements
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 1+il*basenodes;
        TopolCube[2] = 2+il*basenodes;
        TopolCube[3] = 3+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 1+(il+1)*basenodes;
        TopolCube[6] = 2+(il+1)*basenodes;
        TopolCube[7] = 3+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 3+il*basenodes;
        TopolCube[2] = 7+il*basenodes;
        TopolCube[3] = 4+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 3+(il+1)*basenodes;
        TopolCube[6] = 7+(il+1)*basenodes;
        TopolCube[7] = 4+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 0+il*basenodes;
        TopolCube[1] = 1+il*basenodes;
        TopolCube[2] = 5+il*basenodes;
        TopolCube[3] = 4+il*basenodes;
        TopolCube[4] = 0+(il+1)*basenodes;
        TopolCube[5] = 1+(il+1)*basenodes;
        TopolCube[6] = 5+(il+1)*basenodes;
        TopolCube[7] = 4+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 3+il*basenodes;
        TopolCube[1] = 2+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 7+il*basenodes;
        TopolCube[4] = 3+(il+1)*basenodes;
        TopolCube[5] = 2+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 7+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 1+il*basenodes;
        TopolCube[1] = 2+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 5+il*basenodes;
        TopolCube[4] = 1+(il+1)*basenodes;
        TopolCube[5] = 2+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 5+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
        TopolCube[0] = 4+il*basenodes;
        TopolCube[1] = 5+il*basenodes;
        TopolCube[2] = 6+il*basenodes;
        TopolCube[3] = 7+il*basenodes;
        TopolCube[4] = 4+(il+1)*basenodes;
        TopolCube[5] = 5+(il+1)*basenodes;
        TopolCube[6] = 6+(il+1)*basenodes;
        TopolCube[7] = 7+(il+1)*basenodes;
        geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
        id++;
        
    }
    
    geomesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    int axis = 1;
    REAL angle = 0.0;//-45.0;
    RotateGeomesh(geomesh, angle, axis);

    return geomesh;
}

/** This mesh exlude one sixth of the original mesh, just in the region that the solution 6 is almost constant */
TPZGeoMesh * MakeSphereFromQuadrilateralFacesR(int ndiv, SimulationCase  & sim_data)
{

#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Sphere:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    
    if (sim_data.UseGmshMeshQ) {
        
        std::string dirname = PZSOURCEDIR;
        std::string grid;
        
        switch (sim_data.elemen_type) {
            case 0:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryTQ.msh";
            }
                break;
            case 1:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryPQ.msh";
            }
                break;
            case 2:
            {
                grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryHQ.msh";
            }
                break;
            default:
            {
                DebugStop();
            }
                break;
        }
        
        
        TPZGmshReader Geometry;
        REAL s = 1.0;
        Geometry.SetfDimensionlessL(s);
        geomesh = Geometry.GeometricGmshMesh(grid);
        const std::string name("Sphere from gmsh script");
        geomesh->SetName(name);
        
        // changin id internally
        sim_data.gamma_ids[0] = 2;
        sim_data.gamma_ids[1] = 3;
    }
    else
    {
    
        geomesh->SetDimension(3);
        const int nl = 2;// for this mesh it is fixed
        int basenodes = 8;
        int nodes =  basenodes * (nl) + 12;
        REAL radius_o = 1.0;
        REAL radius_i = 0.25;
        
        REAL dr = (radius_o - radius_i)/REAL(nl-1);
        geomesh->SetMaxNodeId(nodes-1);
        geomesh->NodeVec().Resize(nodes);
        TPZManVector<TPZGeoNode,4> Node(nodes);
        
        TPZManVector<long,1> TopolPoint(1);
        TPZManVector<long,4> TopolQuad(4),TopolArc(3);
        TPZManVector<long,8> TopolCube(8),TopolQuadQua(8);
        TPZManVector<REAL,3> coord(3,0.);
        TPZVec<REAL> xc(3,0.);
        REAL cphi = atan(sqrt(2.0));
        
        int nodeindex = 0;
        long id = 0;
        int matid = sim_data.omega_ids[0];
        
        TPZManVector< TPZManVector<REAL,3> , 8 > points(basenodes,0.);
        TPZManVector< TPZManVector<REAL,3> , 12 > points_q(12,0.);
        for (int il = 0; il < nl; il++) {
            
            if (il==0) {
                matid = sim_data.gamma_ids[0];
            }
            
            if (il==nl-1) {
                matid = sim_data.gamma_ids[1];
            }
            
            REAL radius = radius_o - REAL(il)*dr;
            
            /* 0 */
            points[0].Resize(3, 0.0);
            points[0][0]=radius;
            points[0][1]=cphi;
            points[0][2]=M_PI/4.0;

            /* 1 */
            points[1].Resize(3, 0.0);
            points[1][0]=radius;
            points[1][1]=cphi;
            points[1][2]=-M_PI/4.0;
            
            /* 2 */
            points[2].Resize(3, 0.0);
            points[2][0]=radius;
            points[2][1]=cphi;
            points[2][2]=-3.0*M_PI/4.0;
            
            /* 3 */
            points[3].Resize(3, 0.0);
            points[3][0]=radius;
            points[3][1]=cphi;
            points[3][2]=3.0*M_PI/4.0;
            
            /* 4 */
            points[4].Resize(3, 0.0);
            points[4][0]=radius;
            points[4][1]=M_PI-cphi;
            points[4][2]=M_PI/4.0;
            
            /* 5 */
            points[5].Resize(3, 0.0);
            points[5][0]=radius;
            points[5][1]=M_PI-cphi;
            points[5][2]=-M_PI/4.0;
            
            /* 6 */
            points[6].Resize(3, 0.0);
            points[6][0]=radius;
            points[6][1]=M_PI-cphi;
            points[6][2]=-3.0*M_PI/4.0;
            
            /* 7 */
            points[7].Resize(3, 0.0);
            points[7][0]=radius;
            points[7][1]=M_PI-cphi;
            points[7][2]=3.0*M_PI/4.0;
            

            for (int i = 0; i < basenodes; i++) {
                coord = ParametricSphere(points[i][0],points[i][1],points[i][2]);
                geomesh->NodeVec()[nodeindex].SetCoord(coord);
                geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
                nodeindex++;
            }
            
            if (il==0 || il==nl-1)
            {
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 1+il*basenodes;
                TopolQuad[2] = 2+il*basenodes;
                TopolQuad[3] = 3+il*basenodes;
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
                quad1->Geom().SetData(radius, xc);
                id++;
                
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 3+il*basenodes;
                TopolQuad[2] = 7+il*basenodes;
                TopolQuad[3] = 4+il*basenodes;
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
                quad2->Geom().SetData(radius, xc);
                id++;
                
                TopolQuad[0] = 0+il*basenodes;
                TopolQuad[1] = 1+il*basenodes;
                TopolQuad[2] = 5+il*basenodes;
                TopolQuad[3] = 4+il*basenodes;
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
                quad3->Geom().SetData(radius, xc);
                id++;
                
    //            TopolQuad[0] = 3+il*basenodes;
    //            TopolQuad[1] = 2+il*basenodes;
    //            TopolQuad[2] = 6+il*basenodes;
    //            TopolQuad[3] = 7+il*basenodes;
    //            TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
    //            quad4->Geom().SetData(radius, xc);
    //            id++;
                
                TopolQuad[0] = 1+il*basenodes;
                TopolQuad[1] = 2+il*basenodes;
                TopolQuad[2] = 6+il*basenodes;
                TopolQuad[3] = 5+il*basenodes;
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
                quad5->Geom().SetData(radius, xc);
                id++;
                
                TopolQuad[0] = 4+il*basenodes;
                TopolQuad[1] = 5+il*basenodes;
                TopolQuad[2] = 6+il*basenodes;
                TopolQuad[3] = 7+il*basenodes;
                TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*geomesh);
                quad6->Geom().SetData(radius, xc);
                id++;
            }
            
        }
        

        points_q[0].Resize(3, 0.0);
        points_q[0][0]= radius_o - 0.0*dr;
        points_q[0][1]=M_PI/4.0;
        points_q[0][2]=4.0*M_PI/4.0;
        
        points_q[1].Resize(3, 0.0);
        points_q[1][0]= radius_o - 0.5*dr;
        points_q[1][1]=cphi;
        points_q[1][2]=-3.0*M_PI/4.0;
        
        points_q[2].Resize(3, 0.0);
        points_q[2][0]= radius_o - 1.0*dr;
        points_q[2][1]=M_PI/4.0;
        points_q[2][2]=4.0*M_PI/4.0;
        
        points_q[3].Resize(3, 0.0);
        points_q[3][0]= radius_o - 0.5*dr;
        points_q[3][1]=cphi;
        points_q[3][2]=3.0*M_PI/4.0;
        
        /**/
        
        points_q[4].Resize(3, 0.0);
        points_q[4][0]= radius_o - 0.0*dr;
        points_q[4][1]=0.5*M_PI;
        points_q[4][2]=3.0*M_PI/4.0;
        
        
        points_q[5].Resize(3, 0.0);
        points_q[5][0]= radius_o - 0.5*dr;
        points_q[5][1]=M_PI-cphi;
        points_q[5][2]=3.0*M_PI/4.0;
        
        
        points_q[6].Resize(3, 0.0);
        points_q[6][0]= radius_o - 1.0*dr;
        points_q[6][1]=0.5*M_PI;
        points_q[6][2]=3.0*M_PI/4.0;
        
        /**/
        
        points_q[7].Resize(3, 0.0);
        points_q[7][0]= radius_o - 0.0*dr;
        points_q[7][1]=3.0*M_PI/4.0;
        points_q[7][2]=4.0*M_PI/4.0;
        
        points_q[8].Resize(3, 0.0);
        points_q[8][0]= radius_o - 0.5*dr;
        points_q[8][1]=M_PI-cphi;
        points_q[8][2]=-3.0*M_PI/4.0;
        
        points_q[9].Resize(3, 0.0);
        points_q[9][0]= radius_o - 1.0*dr;
        points_q[9][1]=3.0*M_PI/4.0;
        points_q[9][2]=4.0*M_PI/4.0;
        
        /**/
        
        points_q[10].Resize(3, 0.0);
        points_q[10][0]= radius_o - 0.0*dr;
        points_q[10][1]=0.5*M_PI;
        points_q[10][2]=-3.0*M_PI/4.0;
        
        points_q[11].Resize(3, 0.0);
        points_q[11][0]= radius_o - 1.0*dr;
        points_q[11][1]=0.5*M_PI;
        points_q[11][2]=-3.0*M_PI/4.0;
        
        for (int i = 0; i < points_q.size(); i++) {
            coord = ParametricSphere(points_q[i][0],points_q[i][1],points_q[i][2]);
            geomesh->NodeVec()[nodeindex].SetCoord(coord);
            geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
        
        matid = sim_data.omega_ids[0];
        
        for (int il = 0; il < nl - 1 ; il++) {
            //      Inserting blend elements
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 1+il*basenodes;
            TopolCube[2] = 2+il*basenodes;
            TopolCube[3] = 3+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 1+(il+1)*basenodes;
            TopolCube[6] = 2+(il+1)*basenodes;
            TopolCube[7] = 3+(il+1)*basenodes;
            geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
            id++;
            
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 3+il*basenodes;
            TopolCube[2] = 7+il*basenodes;
            TopolCube[3] = 4+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 3+(il+1)*basenodes;
            TopolCube[6] = 7+(il+1)*basenodes;
            TopolCube[7] = 4+(il+1)*basenodes;
            geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
            id++;
            
            TopolCube[0] = 0+il*basenodes;
            TopolCube[1] = 1+il*basenodes;
            TopolCube[2] = 5+il*basenodes;
            TopolCube[3] = 4+il*basenodes;
            TopolCube[4] = 0+(il+1)*basenodes;
            TopolCube[5] = 1+(il+1)*basenodes;
            TopolCube[6] = 5+(il+1)*basenodes;
            TopolCube[7] = 4+(il+1)*basenodes;
            geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
            id++;
            
            
            // Adding hole boundary
            int auxiliar_id = -10;
            int external_id = -1;
            
            //////////////////////////////////////////////////////////////////////////////
            // outer arcs
            
            TopolArc[0] = 3+il*basenodes;
            TopolArc[1] = 2+il*basenodes;
            TopolArc[2] = 0 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            
            TopolArc[0] = 2+il*basenodes;
            TopolArc[1] = 6+il*basenodes;
            TopolArc[2] = 10 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            TopolArc[0] = 6+il*basenodes;
            TopolArc[1] = 7+il*basenodes;
            TopolArc[2] = 7 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            TopolArc[0] = 3+il*basenodes;
            TopolArc[1] = 7+il*basenodes;
            TopolArc[2] = 4 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            
            //////////////////////////////////////////////////////////////////////////////
            // inner arcs
            
            TopolArc[0] = 3+(il+1)*basenodes;
            TopolArc[1] = 2+(il+1)*basenodes;
            TopolArc[2] = 2 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            TopolArc[0] = 2+(il+1)*basenodes;
            TopolArc[1] = 6+(il+1)*basenodes;
            TopolArc[2] = 11 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            TopolArc[0] = 6+(il+1)*basenodes;
            TopolArc[1] = 7+(il+1)*basenodes;
            TopolArc[2] = 9 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            
            TopolArc[0] = 3+(il+1)*basenodes;
            TopolArc[1] = 7+(il+1)*basenodes;
            TopolArc[2] = 6 + 2*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc, auxiliar_id,*geomesh);
            id++;
            
            //////////////////////////////////////////////////////////////////////////////
            // bc of the hole
            
            TopolQuad[0] = 3+il*basenodes;
            TopolQuad[1] = 2+il*basenodes;
            TopolQuad[2] = 2+(il+1)*basenodes;
            TopolQuad[3] = 3+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (id,TopolQuad, external_id,*geomesh);
            id++;
            
            TopolQuad[0] = 2+il*basenodes;
            TopolQuad[1] = 6+il*basenodes;
            TopolQuad[2] = 6+(il+1)*basenodes;
            TopolQuad[3] = 2+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (id,TopolQuad, external_id,*geomesh);
            id++;
            
            TopolQuad[0] = 6+il*basenodes;
            TopolQuad[1] = 7+il*basenodes;
            TopolQuad[2] = 7+(il+1)*basenodes;
            TopolQuad[3] = 6+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (id,TopolQuad, external_id,*geomesh);
            id++;
            
            TopolQuad[0] = 3+il*basenodes;
            TopolQuad[1] = 7+il*basenodes;
            TopolQuad[2] = 7+(il+1)*basenodes;
            TopolQuad[3] = 3+(il+1)*basenodes;
            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (id,TopolQuad, external_id,*geomesh);
            id++;
            
            
    //        TopolQuadQua[0] = 3+il*basenodes;
    //        TopolQuadQua[1] = 2+il*basenodes;
    //        TopolQuadQua[2] = 2+(il+1)*basenodes;
    //        TopolQuadQua[3] = 3+(il+1)*basenodes;
    //        
    //        TopolQuadQua[4] = 0 + 2*basenodes;
    //        TopolQuadQua[5] = 1 + 2*basenodes;
    //        TopolQuadQua[6] = 2 + 2*basenodes;
    //        TopolQuadQua[7] = 3 + 2*basenodes;
    //        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad  > (id,TopolQuadQua,external_id,*geomesh);
    //        id++;
    //        
    //        
    //        TopolQuadQua[0] = 2+il*basenodes;
    //        TopolQuadQua[1] = 6+il*basenodes;
    //        TopolQuadQua[2] = 6+(il+1)*basenodes;
    //        TopolQuadQua[3] = 2+(il+1)*basenodes;
    //        
    //        TopolQuadQua[4] = 10 + 2*basenodes;
    //        TopolQuadQua[5] = 8 + 2*basenodes;
    //        TopolQuadQua[6] = 11 + 2*basenodes;
    //        TopolQuadQua[7] = 1 + 2*basenodes;
    //        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad  > (id,TopolQuadQua,external_id,*geomesh);
    //        id++;
    //        
    //
    //        TopolQuadQua[0] = 6+il*basenodes;
    //        TopolQuadQua[1] = 7+il*basenodes;
    //        TopolQuadQua[2] = 7+(il+1)*basenodes;
    //        TopolQuadQua[3] = 6+(il+1)*basenodes;
    //        
    //        TopolQuadQua[4] = 7 + 2*basenodes;
    //        TopolQuadQua[5] = 5 + 2*basenodes;
    //        TopolQuadQua[6] = 9 + 2*basenodes;
    //        TopolQuadQua[7] = 8 + 2*basenodes;
    //        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad  > (id,TopolQuadQua,external_id,*geomesh);
    //        id++;
    //
    //        TopolQuadQua[0] = 3+il*basenodes;
    //        TopolQuadQua[1] = 7+il*basenodes;
    //        TopolQuadQua[2] = 7+(il+1)*basenodes;
    //        TopolQuadQua[3] = 3+(il+1)*basenodes;
    //        
    //        TopolQuadQua[4] = 4 + 2*basenodes;
    //        TopolQuadQua[5] = 5 + 2*basenodes;
    //        TopolQuadQua[6] = 6 + 2*basenodes;
    //        TopolQuadQua[7] = 3 + 2*basenodes;
    //        new TPZGeoElRefPattern< pzgeom::TPZQuadraticQuad  > (id,TopolQuadQua,external_id,*geomesh);
    //        id++;
            
            
            TopolCube[0] = 1+il*basenodes;
            TopolCube[1] = 2+il*basenodes;
            TopolCube[2] = 6+il*basenodes;
            TopolCube[3] = 5+il*basenodes;
            TopolCube[4] = 1+(il+1)*basenodes;
            TopolCube[5] = 2+(il+1)*basenodes;
            TopolCube[6] = 6+(il+1)*basenodes;
            TopolCube[7] = 5+(il+1)*basenodes;
            geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
            id++;
            
            TopolCube[0] = 4+il*basenodes;
            TopolCube[1] = 5+il*basenodes;
            TopolCube[2] = 6+il*basenodes;
            TopolCube[3] = 7+il*basenodes;
            TopolCube[4] = 4+(il+1)*basenodes;
            TopolCube[5] = 5+(il+1)*basenodes;
            TopolCube[6] = 6+(il+1)*basenodes;
            TopolCube[7] = 7+(il+1)*basenodes;
            geomesh->CreateGeoBlendElement(ECube, TopolCube, matid, id);
            id++;
            
        }
    
    }
    
    geomesh->BuildConnectivity();
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = ndiv;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    int axis = 1;
    REAL angle = 0.0;//-45.0;
    RotateGeomesh(geomesh, angle, axis);
    
    return geomesh;
}

// Cylinder
TPZGeoMesh * MakeCylinderFromLinearFaces(int ndiv, SimulationCase  & sim_data){
    
    /* GID mesh */
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cylinder:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> inner surface ),  (-2 -> outer surface) } " << std::endl;
        DebugStop();
    }
#endif
    
    bool ExtrudeMeshQ = true;
    
    // changin id internally
    sim_data.gamma_ids[0] = 2;
    sim_data.gamma_ids[1] = 3;
    
    std::string dirname = PZSOURCEDIR;
    std::string grid_name;
//    grid_name = dirname + "/Projects/Hdiv3DCurved/gid_meshes/CircularMeshVerticalWell.dump";
//    grid_name = dirname + "/Projects/Hdiv3DCurved/gid_meshes/CircularMeshVerticalWell3D.dump";
//    grid_name = dirname + "/Projects/Hdiv3DCurved/gid_meshes/CircularMeshVerticalWellQ.dump";
    grid_name = dirname + "/Projects/Hdiv3DCurved/gid_meshes/CircularMeshVerticalWellQII.dump";
    
    TPZReadGIDGrid GeometryInfo;
    REAL s = 1.0;
    GeometryInfo.SetfDimensionlessL(s);
    TPZGeoMesh * geomesh = GeometryInfo.GeometricGIDMesh(grid_name);
    
    long last_node = geomesh->NNodes() - 1;
    long last_element = geomesh->NElements() - 1;
    long node_id = geomesh->NodeVec()[last_node].Id();
    long element_id = geomesh->Element(last_element)->Id();
    const std::string name("GID geometry");
    geomesh->SetName(name);
    geomesh->SetMaxNodeId(node_id);
    geomesh->SetMaxElementId(element_id);
    
    TPZGeoMesh * geomesh_3d;
    if (ExtrudeMeshQ) {
        geomesh->SetDimension(2);
        
        int nel_z = 1;
        TPZManVector<REAL,2> dz(2,nel_z);
        dz[0] = 20.0/REAL(nel_z)/s;
        geomesh_3d = ExtrudedGIDMesh(geomesh, sim_data, dz);
        
        TPZVec<TPZGeoEl *> sons;
        const int nref = ndiv;
        for (int iref = 0; iref < nref; iref++) {
            int nel = geomesh_3d->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *gel = geomesh_3d->ElementVec()[iel];
                if(!gel->HasSubElement())
                {
                    gel->Divide(sons);
                }
            }
        }
        
        geomesh_3d->BuildConnectivity();
        
        int axis = 1;
        REAL angle = 0.0;//-45.0;
        RotateGeomesh(geomesh_3d, angle, axis);
        return geomesh_3d;
        
    }
    else{
        
        TPZVec<TPZGeoEl *> sons;
        const int nref = ndiv;
        for (int iref = 0; iref < nref; iref++) {
            int nel = geomesh->NElements();
            for (int iel = 0; iel < nel; iel++) {
                TPZGeoEl *gel = geomesh->ElementVec()[iel];
                if(!gel->HasSubElement())
                {
                    gel->Divide(sons);
                }
            }
        }
        
        geomesh->BuildConnectivity();
        
        int axis = 1;
        REAL angle = 0.0;//-45.0;
        RotateGeomesh(geomesh, angle, axis);
        return geomesh;
        
    }
    

    
}

TPZGeoMesh * ExtrudedGIDMesh(TPZGeoMesh * gmesh, SimulationCase sim_data, TPZManVector<REAL,2> dz){
    
    
    
    REAL t=0.0;
    REAL dt;
    int n;
    bool IsTetrahedronMeshQ = false;
    
    int bc_B =  sim_data.gamma_ids[1];
    int bc_T =  sim_data.gamma_ids[1];
    
    
    TPZHierarquicalGrid CreateGridFrom2D(gmesh);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(bc_B,bc_T);
    if(IsTetrahedronMeshQ){
        CreateGridFrom2D.SetTriangleExtrusion();
        CreateGridFrom2D.SetTetrahedonExtrusion();
    }
//    else{
////        CreateGridFrom2D.SetTriangleExtrusion();
//        CreateGridFrom2D.SetPrismExtrusion();
//    }

    
    dt = dz[0];
    n = int(dz[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * gmesh_3d = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
    long last_node = gmesh_3d->NNodes() - 1;
    long last_element = gmesh_3d->NElements() - 1;
    long node_id = gmesh_3d->NodeVec()[last_node].Id();
    long element_id = gmesh_3d->Element(last_element)->Id();
    const std::string name("Reservoir with vertical extrusion");
    gmesh_3d->SetName(name);
    gmesh_3d->SetMaxNodeId(node_id);
    gmesh_3d->SetMaxElementId(element_id);
    gmesh_3d->SetDimension(3);
    return gmesh_3d;
    
}

void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}


void TransformToQuadratic(TPZGeoMesh *gmesh)
{
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement()) {
            continue;
        }
        
        int whichsubel =-1;
        TPZGeoEl *father = gel->Father();
        if (gel->HasSubElement()) {
            whichsubel = gel->WhichSubel();
        }

        gel = TPZChangeEl::ChangeToQuadratic(gmesh, el);
        if (whichsubel != -1) {
            father->SetSubElement(whichsubel, gel);
        }
    }
}


TPZManVector<STATE,3> ParametricSphere(REAL radius, REAL theta, REAL phi)
{
    TPZManVector<STATE,3> xcoor(3,0.0);
    xcoor[0] = radius * sin(theta) * cos(phi) ;
    xcoor[1] = radius * sin(theta) * sin(phi) ;
    xcoor[2] = radius * cos(theta) ;
    return xcoor;
}

TPZManVector<STATE,3> ParametricCylinder(REAL radius, REAL theta, REAL z)
{
    TPZManVector<STATE,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = z;
    return xcoor;
}

void PrintGeometry(TPZGeoMesh * gmesh, SimulationCase & sim_data){
    std::stringstream text_name;
    std::stringstream vtk_name;
    text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "mhm" << "_l" << level_mhm << ".txt";
    vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "mhm" << "_l" << level_mhm << ".vtk";
    ofstream textfile(text_name.str());
    gmesh->Print(textfile);
    
    std::ofstream vtkfile(vtk_name.str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

/** @brief Apply uniform refinement on the Geometric mesh */
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gmesh->ElementVec() [i];
            if (gel->Dimension() != 0) gel->Divide (filhos);
        }//for i
    }//ref
    gmesh->BuildConnectivity();
}

void UniformRefineTetrahedrons(TPZGeoMesh *gmesh, int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D;
    
    {// needed!
        char buf[] =
        "10 9 "
        "-50 Tet0000111111111 "
        "0 0 0 "
        "1 0 0 "
        "0 1 0 "
        "0 0 1 "
        "0.5 0 0 "
        "0 0.5 0 "
        "0 0 0.5 "
        "0.5 0.5 0 "
        "0 0.5 0.5 "
        "0.5 0 0.5 "
        "4 4 0  1  2  3 "
        "4 4 0  4  5  6 "
        "4 4 4  1  7  9 "
        "4 4 7  2  5  8 "
        "4 4 6  9  8  3 "
        "4 4 4  9  6  5 "
        "4 4 5  8  6  9 "
        "4 4 7  8  9  5 "
        "4 4 4  7  5  9 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    //    TPZAutoPointer<TPZRefPattern> refp3D = gRefDBase.FindRefPattern("UnifTet");
    
    if(!refp3D)
    {
        DebugStop();
    }
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < n_ref; r++)
    {
        int nels = gmesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==3)
            {
                gel->SetRefPattern(refp3D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            if(gel->Dimension()==2)
            {
                //                gel->SetRefPattern(refp3D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            
        }
    }
}

void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the axis -> i.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void ErrorH1(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_h1)
{
    
    long nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerror(3,0.   );
    for (long iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = cmesh->ElementVec()[iel];

        if (!cel) {
            continue;
        }
        
        if(cel->Reference()->Dimension()!=dim) {
            continue;
        }
        
        TPZManVector<STATE,10> elerror(3,0.);
        elerror.Fill(0.);
        cel->EvaluateError(Analytic, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerror[i] += elerror[i]*elerror[i];
            
        }
    }
    
    error_primal    = sqrt(globalerror[0]);
    error_dual      = sqrt(globalerror[1]);
    error_h1        = sqrt(globalerror[2]);
    
    
}

void ErrorHdiv(TPZCompMesh *cmesh, REAL &error_primal , REAL & error_dual, REAL & error_hdiv){
    
    long nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerror(3,0.0);
    
    for (long iel = 0; iel < nel; iel++) {
        
        TPZCompEl *cel = cmesh->ElementVec()[iel];
        if (!cel) {
            continue;
        }
        
        if(cel->Reference()->Dimension()!=dim) {
            continue;
        }
        
        TPZManVector<STATE,10> elerror(3,0.);
        elerror.Fill(0.);
        cel->EvaluateError(Analytic, elerror, NULL);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globalerror[i] += elerror[i]*elerror[i];
            
        }
        
    }
    
    error_primal    = sqrt(globalerror[0]);
    error_dual      = sqrt(globalerror[1]);
    error_hdiv      = sqrt(0.0*globalerror[1] + globalerror[2]);
    
}

/// uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (long el=0; el<nel; el++) {
            
            TPZCompEl *cel = cmesh->Element(el);
            TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
            if (condense) {
                condense->Unwrap();
                change = true;
            }
            cel = cmesh->Element(el);
            TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
            if (elgr) {
                elgr->Unwrap();
                change = true;
            }
        }
    }
}


/** @brief Sparated connects by hdiv connect neighborhood */
void SeparateConnectsByNeighborhood(TPZCompMesh * mixed_cmesh){
    
#ifdef PZDEBUG
    if(!mixed_cmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh *gmesh = mixed_cmesh->Reference();
    gmesh->ResetReference();
    mixed_cmesh->LoadReferences();
    mixed_cmesh->ComputeNodElCon();
    long nel = mixed_cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = mixed_cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (!gel || (gel->Dimension() != gmesh->Dimension()) ) {
            continue;
        }
        int nc = cel->NConnects();
        for (int ic =0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() && c.NElConnected() == 2) // @omar:: Hdiv connects have this invariant characteristic
            {
                // duplicate the connect
                long cindex = mixed_cmesh->AllocateNewConnect(c);
                TPZConnect &newc = mixed_cmesh->ConnectVec()[cindex];
                newc = c;
                c.DecrementElConnected();
                newc.DecrementElConnected();
                cel->SetConnectIndex(ic, cindex);
            }
        }
    }
    mixed_cmesh->ExpandSolution();
}

/** @brief Build MHM form the current hdvi mesh */
void InsertSkeletonInterfaces(TPZGeoMesh * gmesh){
    
#ifdef PZDEBUG
    if(!gmesh){
        DebugStop();
    }
#endif
    
    int Skeleton_material_Id = 100;
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Level() != level_mhm || gel->Dimension() != gmesh->Dimension()) {
            continue;
        }
        
//        if (!gel || gel->HasSubElement() || gel->Dimension() != gmesh->Dimension()) {
//            continue;
//        }
        
        int nsides = gel->NSides();
        for (int is = gel->NCornerNodes(); is<nsides; is++) {
            if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, Skeleton_material_Id);
            }
        }
    }
    gmesh->BuildConnectivity();
}

/** @brief Construc computational macro elements */
void BuildMacroElements(TPZCompMesh * mixed_cmesh)
{
    
    std::cout << "ndof Elliptic before MHM substructuring = " << mixed_cmesh->Solution().Rows() << std::endl;
    
#ifdef PZDEBUG
    if(!mixed_cmesh){
        DebugStop();
    }
#endif
    
    bool KeepOneLagrangian = true;
    typedef std::set<long> TCompIndexes;
    std::map<long, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = mixed_cmesh->Reference();
    gmesh->ResetReference();
    mixed_cmesh->LoadReferences();
    long nelg = gmesh->NElements();
    for (long el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Father() != NULL) {
            continue;
        }
        if (gel->Dimension() == gmesh->Dimension() - 1) {
            continue;
        }
        long mapindex = gel->Index();
        if (gel->Dimension() == gmesh->Dimension() - 1) {
            TPZGeoElSide neighbour = gel->Neighbour(gel->NSides()-1);
            if (neighbour.Element()->Dimension() != gmesh->Dimension()) {
                DebugStop();
            }
            mapindex= neighbour.Element()->Index();
        }
        TPZStack<TPZCompElSide> highlevel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.HigherLevelCompElementList3(highlevel, 0, 0);
        long nelst = highlevel.size();
        for (long elst=0; elst<nelst; elst++) {
            ElementGroups[mapindex].insert(highlevel[elst].Element()->Index());
        }
        if (gel->Reference()) {
            if (nelst) {
                DebugStop();
            }
            ElementGroups[mapindex].insert(gel->Reference()->Index());
        }
    }
    
#ifdef PZDEBUG
    std::map<long,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " with size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<long>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
#endif
    
    std::set<long> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(mixed_cmesh, ElementGroups, submeshindices, KeepOneLagrangian);
    
    std::cout << "Inserting " << ElementGroups.size()  <<  " macro elements into MHM substructures" << std::endl;
    mixed_cmesh->ComputeNodElCon();
    mixed_cmesh->CleanUpUnconnectedNodes();
    for (std::set<long>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = mixed_cmesh->Element(*it);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        subcmesh->SetAnalysisSkyline(0, 0, 0);
    }
    mixed_cmesh->ComputeNodElCon();
    mixed_cmesh->CleanUpUnconnectedNodes();
    mixed_cmesh->ExpandSolution();
    
}
