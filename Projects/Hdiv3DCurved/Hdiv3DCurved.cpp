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
#include "pzintel.h"


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


struct SimulationCase {
    bool            IsHdivQ;
    bool            IsMHMQ;
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    bool            UseGmshMeshQ;
    bool            NonAffineQ;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             n_acc_terms;
    int             int_order;
    int             n_threads;
    int             perturbation_type;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   gamma_ids;
    
    SimulationCase() : IsHdivQ(false), IsMHMQ(false), UsePardisoQ(true), UseFrontalQ(false), UseGmshMeshQ(false), NonAffineQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),gamma_ids()
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) : IsHdivQ(copy.IsHdivQ), IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
        UseGmshMeshQ(copy.UseGmshMeshQ), NonAffineQ(copy.NonAffineQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
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
        NonAffineQ = copy.NonAffineQ;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        n_acc_terms = copy.n_acc_terms;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        perturbation_type = copy.perturbation_type;
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
#define Solution6
//#define Thiem

// Solutions for non-affine meshes
//#define Solution7

// Arctan solution
//#define Solution8

// MHM rates subtructuring level
static int level_mhm = 0;

static void Analytic(const TPZVec<REAL> &x, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu);
static void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &f);   ///Jorhge 2017. It is not used: , TPZFMatrix<STATE> &gradf);
static void f(const TPZVec<REAL> &p, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf);

TPZGeoMesh * GeomtricMesh(int ndiv, SimulationCase  & sim_data);
void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, REAL & rho_min);
void PrintGeometry(TPZGeoMesh * gmesh, SimulationCase & sim_data);
void PrintGeometryVols(TPZGeoMesh * gmesh, std::stringstream & file_name);
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref);
void UniformRefineTetrahedrons(TPZGeoMesh * gmesh, int n_ref);
void RefineHexahedronsToTetrahedrons(TPZGeoMesh * gmesh, int n_ref);
void RefineTetrahedronsToHexahedrons(TPZGeoMesh * gmesh, int n_ref);
void RefineHexahedronsToNonAffineHexahedrons(TPZGeoMesh * gmesh, int n_ref);
void RefineHexahedronsToPrisms(TPZGeoMesh * gmesh, int n_ref);
void RefineHexahedronsToNonAffinePrisms(TPZGeoMesh * gmesh, int n_ref);

TPZGeoMesh * MakeCubeFromLinearQuadrilateralFaces(int ndiv, SimulationCase  & sim_data);
TPZGeoMesh * MakeCubeFromLinearTriangularFaces(int ndiv, SimulationCase  & sim_data);
void Parametricfunction_x(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction_y(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Parametricfunction_z(const TPZVec<REAL> &par, TPZVec<REAL> &X);

// Meshes for the publication
// H(div) finite elements based on non-affine meshes for three dimensional mixed formulations
// of flow problems with arbitrary high order flux divergence accuracy
TPZGeoMesh * MakeCube(SimulationCase  & sim_data);
TPZGeoMesh * MakeCubeFromTetrahedrons(int ndiv, SimulationCase  & sim_data);
TPZGeoMesh * MakeCubeFromHexahedrons(int ndiv, SimulationCase  & sim_data);
TPZGeoMesh * MakeCubeFromPrisms(int ndiv, SimulationCase  & sim_data);

TPZGeoMesh * MakeSphereFromLinearQuadrilateralFaces(int ndiv, SimulationCase  & sim_data);
TPZGeoMesh * MakeSphereFromQuadrilateralFaces(int ndiv, SimulationCase  & sim_data);

TPZGeoMesh * MakeSphereFromQuadrilateralFacesR(int ndiv, SimulationCase  & sim_data);

// Cylinder
TPZGeoMesh * MakeCylinderFromLinearFaces(int ndiv, SimulationCase & sim_data);
TPZGeoMesh * ExtrudedGIDMesh(TPZGeoMesh * gmesh, SimulationCase sim_data, TPZManVector<REAL,2> dz);
void ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X);

TPZManVector<REAL,3> ParametricSphere(REAL radius,REAL theta,REAL phi);

TPZManVector<REAL,3> ParametricCylinder(REAL radius,REAL theta,REAL z);


void TransformToQuadratic(TPZGeoMesh *gmesh);
void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
void PerturbateNodes_I(TPZGeoMesh *gmesh, int ndiv);
void PertubationMatrix_I(TPZManVector<REAL> CoordX, REAL pertub_param, REAL ElSize, TPZManVector<REAL> &CoordsPertubated);

TPZCompMesh *DualMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, int64_t &ndof, TPZVec<TPZCompMesh *> &meshvec);

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, int64_t &ndof); //  Primal approximation
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
void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &h_size, TPZVec<REAL> &convergence);

STATE IntegrateVolume(TPZGeoMesh * geometry, SimulationCase sim_data);
STATE IntegrateSolution(TPZCompMesh * cmesh,  SimulationCase sim_data);

// MHM utilities
/** @brief Sparated connects by hdiv connect neighborhood */
void SeparateConnectsByNeighborhood(TPZCompMesh * mixed_cmesh);

/** @brief Insert low dimentional elements defining the skeleton */
void InsertSkeletonInterfaces(TPZGeoMesh * gmesh);

/** @brief Build mhm macro elements following the mixed sense (space constrains) */
void BuildMacroElements(TPZCompMesh * mixed_cmesh);

void ErrorH1(TPZAnalysis * analysis, REAL &error_primal , REAL & error_dual, REAL & error_h1);
void ErrorHdiv(TPZAnalysis * analysis, REAL &error_primal , REAL & error_dual, REAL & error_hdiv);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.main"));
#endif

/**
 *  Configuration for the publication:
 *  Three dimensional hierarchical mixed finite element approximations with
 *  enhanced primal variable accuracy
 */
void Configuration_Affine();

/**
 *  Configuration for the publication:
 *  Enhanced mixed finite element approximations for 3D elliptic problems 
 *  based on non-affine hexahedral and prismatic meshes
 */
void Configuration_Non_Affine();

int main()
{

  HDivPiola = 1;
    
  gRefDBase.InitializeAllUniformRefPatterns();
    //	gRefDBase.InitializeRefPatterns();
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    
    
    bool IsAffineSettingQ = false;
    
    
    if (IsAffineSettingQ) {
        // Runing the affine meshes
        Configuration_Affine();
    }
    else{
        // Runing the non affine meshes
        Configuration_Non_Affine();
    }
    
    
    
    
    return 0;
    
}

void Configuration_Affine(){
    
    TPZStack<SimulationCase> simulations;
    
    // Formulations over the cube
    struct SimulationCase common;
    
    common.UsePardisoQ = true;
    common.UseFrontalQ = false;
    common.UseGmshMeshQ = true;
    common.n_h_levels = 3;
    common.n_p_levels = 2;
    common.int_order  = 10;
    common.n_threads  = 10;
    common.NonAffineQ = false;
    common.domain_type = "cube";
    common.conv_summary = "convergence_summary";
    common.omega_ids.Push(1);     // Domain
    common.gamma_ids.Push(-1);    // Gamma_D outer surface
    common.gamma_ids.Push(-2);    // Gamma_D inner surface
    
    //     // Primal Formulation over the solid cube
    struct SimulationCase H1Case_1 = common;
    H1Case_1.IsHdivQ = false;
    H1Case_1.mesh_type = "linear";
    H1Case_1.elemen_type = 0;
    H1Case_1.dump_folder = "H1_T_affine_cube";
//    simulations.Push(H1Case_1);
    H1Case_1.elemen_type = 1;
    H1Case_1.dump_folder = "H1_H_affine_cube";
//    simulations.Push(H1Case_1);
    H1Case_1.elemen_type = 2;
    H1Case_1.dump_folder = "H1_P_affine_cube";
    simulations.Push(H1Case_1);
    H1Case_1.elemen_type = 3;
    H1Case_1.dump_folder = "H1_Hybrid_affine_cube";
//    simulations.Push(H1Case_1);
    
//    //    // Dual Formulation n = 1
//    struct SimulationCase HdivCase_1 = common;
//    HdivCase_1.IsHdivQ = true;
//    HdivCase_1.mesh_type = "linear";
//    HdivCase_1.n_acc_terms = 1;
//    HdivCase_1.elemen_type = 2;
//    HdivCase_1.dump_folder = "Hdiv_n_1_P_affine_cube";
//    simulations.Push(HdivCase_1);
    
    ComputeCases(simulations);
    
    
}

void Configuration_Non_Affine(){
    
    TPZStack<SimulationCase> simulations;
    bool IsNonAffineQ = true;
    // Formulations over the cube
    struct SimulationCase common;

    common.UsePardisoQ = true;
    common.UseFrontalQ = false;
    common.UseGmshMeshQ = true;
    common.n_h_levels = 1;
    common.n_p_levels = 1;
    common.int_order  = 4;
    common.n_threads  = 8;
    common.NonAffineQ = true;
    common.domain_type = "cube";
    common.conv_summary = "convergence_summary";
    common.omega_ids.Push(1);     // Domain
    common.gamma_ids.Push(-1);    // Gamma_D outer surface
    common.gamma_ids.Push(-2);    // Gamma_D inner surface
    
    if (IsNonAffineQ) {
        
//        //     // Primal Formulation over the solid cube
//        struct SimulationCase H1Case_1 = common;
//        H1Case_1.IsHdivQ = false;
//        H1Case_1.mesh_type = "linear";
//        H1Case_1.elemen_type = 1;
//        H1Case_1.dump_folder = "H1_H_non_affine_cube";
//        simulations.Push(H1Case_1);

        //    // Dual Formulation n = 0
        struct SimulationCase HdivCase_1 = common;
        HdivCase_1.IsHdivQ = true;
        HdivCase_1.mesh_type = "linear";
        HdivCase_1.n_acc_terms = 0;
        HdivCase_1.elemen_type = 1;
        HdivCase_1.dump_folder = "Hdiv_n_0_H_non_affine_cube";
        simulations.Push(HdivCase_1);
//
//        //    // Dual Formulation n = 1
//        struct SimulationCase HdivCase_2 = common;
//        HdivCase_2.IsHdivQ = true;
//        HdivCase_2.mesh_type = "linear";
//        HdivCase_2.n_acc_terms = 1;
//        HdivCase_2.elemen_type = 1;
//        HdivCase_2.dump_folder = "Hdiv_n_1_H_non_affine_cube";
//        simulations.Push(HdivCase_2);
//
//        //    // Dual Formulation n = 2
//        struct SimulationCase HdivCase_3 = common;
//        HdivCase_3.IsHdivQ = true;
//        HdivCase_3.mesh_type = "linear";
//        HdivCase_3.n_acc_terms = 2;
//        HdivCase_3.elemen_type = 1;
//        HdivCase_3.dump_folder = "Hdiv_n_2_H_non_affine_cube";
//        simulations.Push(HdivCase_3);
        
    }
    else{
        
        common.UsePardisoQ = true;
        common.UseFrontalQ = false;
        common.UseGmshMeshQ = true;
        common.n_h_levels = 2;
        common.n_p_levels = 2;
        common.int_order  = 10;
        common.n_threads  = 10;
        common.NonAffineQ = IsNonAffineQ;
        common.domain_type = "cube";
        common.conv_summary = "convergence_summary";
        common.omega_ids.Push(1);     // Domain
        common.gamma_ids.Push(-1);    // Gamma_D outer surface
        common.gamma_ids.Push(-2);    // Gamma_D inner surface
        
//        //     // Primal Formulation over the solid cube
//        struct SimulationCase H1Case_1 = common;
//        H1Case_1.IsHdivQ = false;
//        H1Case_1.mesh_type = "linear";
//        H1Case_1.elemen_type = 0;
//        H1Case_1.dump_folder = "H1_T_affine_cube";
//        simulations.Push(H1Case_1);
//        H1Case_1.elemen_type = 1;
//        H1Case_1.dump_folder = "H1_H_affine_cube";
//        simulations.Push(H1Case_1);
//        H1Case_1.elemen_type = 2;
//        H1Case_1.dump_folder = "H1_P_affine_cube";
//        simulations.Push(H1Case_1);
//
//        //    // Dual Formulation n = 0
//        struct SimulationCase HdivCase_1 = common;
//        HdivCase_1.IsHdivQ = true;
//        HdivCase_1.mesh_type = "linear";
//        HdivCase_1.n_acc_terms = 0;
//        HdivCase_1.elemen_type = 0;
//        HdivCase_1.dump_folder = "Hdiv_n_0_T_affine_cube";
//        simulations.Push(HdivCase_1);
//        HdivCase_1.elemen_type = 1;
//        HdivCase_1.dump_folder = "Hdiv_n_0_H_affine_cube";
//        simulations.Push(HdivCase_1);
//        HdivCase_1.elemen_type = 2;
//        HdivCase_1.dump_folder = "Hdiv_n_0_P_affine_cube";
//        simulations.Push(HdivCase_1);
//
//        //    // Dual Formulation n = 1
//        struct SimulationCase HdivCase_2 = common;
//        HdivCase_2.IsHdivQ = true;
//        HdivCase_2.mesh_type = "linear";
//        HdivCase_2.n_acc_terms = 1;
//        HdivCase_2.elemen_type = 0;
//        HdivCase_2.dump_folder = "Hdiv_n_1_T_affine_cube";
//        simulations.Push(HdivCase_2);
//        HdivCase_2.elemen_type = 1;
//        HdivCase_2.dump_folder = "Hdiv_n_1_H_affine_cube";
//        simulations.Push(HdivCase_2);
//        HdivCase_2.elemen_type = 2;
//        HdivCase_2.dump_folder = "Hdiv_n_1_P_affine_cube";
//        simulations.Push(HdivCase_2);
//
//        //    // Dual Formulation n = 2
//        struct SimulationCase HdivCase_3 = common;
//        HdivCase_3.IsHdivQ = true;
//        HdivCase_3.mesh_type = "linear";
//        HdivCase_3.n_acc_terms = 2;
//        HdivCase_3.elemen_type = 0;
//        HdivCase_3.dump_folder = "Hdiv_n_2_T_affine_cube";
//        simulations.Push(HdivCase_3);
//        HdivCase_3.elemen_type = 1;
//        HdivCase_3.dump_folder = "Hdiv_n_2_H_affine_cube";
//        simulations.Push(HdivCase_3);
//        HdivCase_3.elemen_type = 2;
//        HdivCase_3.dump_folder = "Hdiv_n_2_P_affine_cube";
//        simulations.Push(HdivCase_3);
        
        //    // Dual Formulation n = 3
        struct SimulationCase HdivCase_4 = common;
        HdivCase_4.IsHdivQ = true;
        HdivCase_4.mesh_type = "linear";
        HdivCase_4.n_acc_terms = 3;
        HdivCase_4.elemen_type = 0;
        HdivCase_4.dump_folder = "Hdiv_n_3_T_affine_cube";
//        simulations.Push(HdivCase_4);
        HdivCase_4.elemen_type = 1;
        HdivCase_4.dump_folder = "Hdiv_n_3_H_affine_cube";
        simulations.Push(HdivCase_4);
        HdivCase_4.elemen_type = 2;
        HdivCase_4.dump_folder = "Hdiv_n_3_P_affine_cube";
//        simulations.Push(HdivCase_4);
        
    }
    
    ComputeCases(simulations);
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
    
#ifdef USING_BOOST
    std::cout << "End:: Overal time = " << int_t2-int_t1 << std::endl;
#endif

    
}

void ComputeApproximation(SimulationCase & sim_data){
    
#ifdef USING_BOOST
    boost::posix_time::ptime int_case_t1 = boost::posix_time::microsec_clock::local_time();
#endif
    
    // Creating the directory
    std::string command = "mkdir " + sim_data.dump_folder;
    int result = system(command.c_str());
    if(result != 0)
    {
        std::cout << "Executing " << command << " returned result = " << result << std::endl;
    }
    std::stringstream summary;
    summary   << sim_data.dump_folder << "/" "conv" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << ".txt";

    std::ofstream convergence(summary.str().c_str(),std::ios::app);
    
    TPZManVector<REAL,10> p_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> d_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> h_error(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> h_size(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> rho_size(sim_data.n_h_levels+1,1.0);
    TPZManVector<REAL,10> sigma_size(sim_data.n_h_levels+1,1.0);
    
    TPZManVector<REAL,10> p_conv(sim_data.n_h_levels,0.0);
    TPZManVector<REAL,10> d_conv(sim_data.n_h_levels,0.0);
    TPZManVector<REAL,10> h_conv(sim_data.n_h_levels,0.0);
    
    int n_h_levels = sim_data.n_h_levels;
    int n_p_levels = sim_data.n_p_levels;

    using namespace std;
    
    for (int p = 1; p <= n_p_levels; p++) {
        
        convergence << std::endl;        
        convergence << " Polynomial order  =  " << p << std::endl;
        convergence << std::setw(5)  << " ndiv" << std::setw(20)  << " h" << std::setw(15) << " rho" << std::setw(15) << " sigma" << std::setw(15) << " ndof" << std::setw(15) << " ndof_cond" << std::setw(25) << " assemble_time (msec)" << std::setw(25) << " solving_time (msec)" << std::setw(25) << " error_time (msec)" << std::setw(25) << " Primal l2 error" << std::setw(25) << " Dual l2 error"  << std::setw(25) << " H error (H1 or Hdiv)" << std::endl;
        
        int h_base = 0;
        if (sim_data.IsMHMQ) {
            h_base = n_h_levels;
        }
        for (int h = 0; h <= n_h_levels; h++) {
            
            // Compute the geometry
            TPZGeoMesh * gmesh;
            gmesh = GeomtricMesh(h, sim_data);
    
#ifdef PZDEBUG2
            TPZCheckGeom check(gmesh);
            int checkQ = check.PerformCheck();
            if (checkQ) {
                DebugStop();
            }
#endif
            
            if (sim_data.IsMHMQ) {
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
            std::stringstream vols_name;
            std::stringstream vtk_name;
        
#ifdef PZDEBUG
            text_name   << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".txt";
            vols_name   << sim_data.dump_folder << "/" "geo_vols" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".txt";
            vtk_name    << sim_data.dump_folder << "/" "geo" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";
            ofstream textfile(text_name.str().c_str());
            gmesh->Print(textfile);
            
            std::ofstream vtkfile(vtk_name.str().c_str());
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
            PrintGeometryVols(gmesh, vols_name);
#endif
            
            // Compute the geometry
            int64_t ndof, ndof_cond;
            TPZManVector<TPZCompMesh *,5> meshvec;
            
            std::cout << "Allocating computational mesh\n";
            TPZCompMesh * cmesh = ComputationalMesh(gmesh, p, sim_data, ndof, meshvec);

            std::cout << "Create analysis\n";
            // Create Analysis
            TPZAnalysis * analysis = CreateAnalysis(cmesh,sim_data);
            
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "BEFORE\n";
                cmesh->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef USING_BOOST
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
            std::cout << "Assembly neq = " << cmesh->NEquations() << "\n";
            analysis->Assemble();
            ndof_cond = cmesh->NEquations();
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "AFTER\n";
                cmesh->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            

#ifdef USING_BOOST
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
            std::cout << "Solution of the system\n";
            analysis->Solve();
            
#ifdef USING_BOOST
            boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
            REAL assemble_time = boost::numeric_cast<double>((t2-t1).total_milliseconds());
            REAL solving_time  = boost::numeric_cast<double>((t3-t2).total_milliseconds());
#endif
            
#ifdef USING_BOOST
                boost::posix_time::ptime int_unwrap_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            std::cout << "Post processing\n";
                UnwrapMesh(cmesh);
                analysis->LoadSolution();
                cmesh->Solution() *= -1.0; /* consequence of newton correction */
                analysis->LoadSolution(cmesh->Solution());
                TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmesh);
                
#ifdef USING_BOOST
                boost::posix_time::ptime int_unwrap_t2 = boost::posix_time::microsec_clock::local_time();
#endif
                

#ifdef USING_BOOST
            std::cout << "StaticCondensation::Time for uncondense equations = " << int_unwrap_t2-int_unwrap_t1 <<std::endl;
#endif

#ifdef PZDEBUG
            std::stringstream file_name;
            file_name   << sim_data.dump_folder << "/" << "Primal_cmesh" << ".txt";
            std::ofstream sout(file_name.str().c_str());
            cmesh->Print(sout);
#endif
            
            // PostProccessing
            std::stringstream sol_vtk_name;
#ifdef PZDEBUG
            sol_vtk_name    << sim_data.dump_folder << "/" "sol" << "_" << sim_data.mesh_type << "_" << sim_data.domain_type << "_" << "p" << p << "h" <<  h << ".vtk";

            std::string file(sol_vtk_name.str());
            PosProcess(analysis, file, sim_data);
#endif

            
            
            // compute the error
#ifdef USING_BOOST
            boost::posix_time::ptime error_t1 = boost::posix_time::microsec_clock::local_time();
#endif
            
            if (sim_data.IsHdivQ) {
                ErrorHdiv(analysis, p_error[h], d_error[h], h_error[h]);
            }
            else{
                ErrorH1(analysis, p_error[h], d_error[h], h_error[h]);
            }
            
            ComputeCharacteristicHElSize(gmesh,h_size[h],rho_size[h]);
            sigma_size[h] = h_size[h]/rho_size[h];
            
#ifdef USING_BOOST
            boost::posix_time::ptime error_t2 = boost::posix_time::microsec_clock::local_time();
#endif

            
#ifdef USING_BOOST
            REAL error_time = boost::numeric_cast<double>((error_t2 - error_t1).total_milliseconds());
            
            // current summary
            convergence << setw(5) << h << setw(20) << h_size[h] << setw(15) << rho_size[h] << setw(15) << sigma_size[h] << setw(15) << ndof << setw(15) << ndof_cond << setw(25) << assemble_time << setw(25) << solving_time << setw(25) << error_time << setw(25) << p_error[h] << setw(25) << d_error[h]  << setw(25) << h_error[h] << endl;
#endif
            
            
            analysis->CleanUp();
            delete cmesh;
            for (int i = 0; i < meshvec.size(); i++) {
                meshvec[i]->CleanUp();
                delete meshvec[i];
            }
            delete gmesh;
            
        }
        
        // compute rates
        ComputeConvergenceRates(p_error,h_size,p_conv);
        ComputeConvergenceRates(d_error,h_size,d_conv);
        ComputeConvergenceRates(h_error,h_size,h_conv);
        
        
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
    REAL case_solving_time = boost::numeric_cast<double>((int_case_t2-int_case_t1).total_milliseconds());
#endif
    
#ifdef USING_BOOST
    convergence <<  "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << "Case closed in :" << setw(10) <<  case_solving_time/1000.0 << setw(5)   << " seconds." << std::endl;
    std::cout << std::endl;
#endif
    
}

void ComputeConvergenceRates(TPZVec<REAL> &error, TPZVec<REAL> &h_size, TPZVec<REAL> &convergence){
    
    int ndata = error.size();
    STATE log0p5 = log(0.5);
    for (int i = 1; i < ndata; i++) {
        STATE logerror = log(error[i-1]);
        STATE logerrori = log(error[i]);
        convergence[i-1] = (logerrori - logerror)/(log(h_size[i])-log(h_size[i-1]));
    }
}

void ComputeCharacteristicHElSize(TPZGeoMesh * geometry, REAL & h_min, REAL & rho_min){
    
    h_min   = 1.0;
    rho_min = 1.0;
    
    REAL h;
    REAL rho;
    int nel = geometry->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = geometry->Element(iel);
        
#ifdef PZDEBUG
        if(!gel){
            DebugStop();
        }
#endif
        if (gel->Dimension() != geometry->Dimension() || gel->HasSubElement() == 1) {
            continue;
        }
        
        h = gel->CharacteristicSize();
        rho = 2.0*gel->ElementRadius();
        
        if (h < h_min) {
            h_min = h;
        }
        
        if (rho < rho_min) {
            rho_min = rho;
        }
    }
    
}

void Analytic(const TPZVec<REAL> &p, TPZVec<STATE> &u,TPZFMatrix<STATE> &gradu){
    
    gradu.Resize(4,1);
    
    STATE x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
/*    STATE r = sqrt(x*x+y*y+z*z);
    STATE theta = acos(z/r);
    STATE phi = atan2(y,x);
    
    STATE costheta = cos(theta);
    STATE cosphi = cos(phi);
    STATE sintheta = sin(theta);
    STATE sinphi = sin(phi);
    
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
  */
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
    

#ifdef Solution7
    
    REAL cos_pi_x = cos(M_PI*x);
    REAL cos_pi_y = cos(M_PI*y);
    REAL cos_pi_z = cos(M_PI*z);
    
    REAL sin_pi_x = sin(M_PI*x);
    REAL sin_pi_y = sin(M_PI*y);
    REAL sin_pi_z = sin(M_PI*z);
    
    u[0] = sin_pi_x * sin_pi_y * sin_pi_z;
    
    gradu(0,0) = -1.0*(M_PI * cos_pi_x * sin_pi_y * sin_pi_z);
    gradu(1,0) = -1.0*(M_PI * sin_pi_x * cos_pi_y * sin_pi_z);
    gradu(2,0) = -1.0*(M_PI * sin_pi_x * sin_pi_y * cos_pi_z);
    
    gradu(3,0) = 3.0 * M_PI * M_PI * sin_pi_x * sin_pi_y * sin_pi_z;
    
#endif
    

#ifdef Solution8 // arctan solution
    
    gradu(0,0)=0., gradu(1,0)=0., gradu(2,0)=0., gradu(3,0)=0.;
    
    REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
    REAL r0 = M_PI/3.;
    REAL alpha = 5.;
    
    REAL temp1;
    REAL temp2;
    REAL r;
    REAL grad;
    temp1 = 0.,temp2=0., r=0., grad=0.;
    
    //solucao u
    temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    r = sqrt(temp1);
    temp2 = (r - r0)*alpha;
    u[0] = M_PI/2. - atan(temp2);
    
    //fluxo em x
    temp2 =  r*(1./alpha + (r - r0)*(r - r0)*alpha);
    temp1 = (x0-x);
    grad=temp1/temp2;
    grad *=-1.;
    gradu(0,0) = grad;
    
    //fluxo em y
    temp1 = (y0-y);
    grad=temp1/temp2;
    grad *=-1.;
    gradu(1,0) = grad;
    
    
    //fluxo em z
    temp1 = (z0-z);
    grad=temp1/temp2;
    grad *=-1.;
    gradu(2,0) = grad;
    
    //Solucao do divergente de u
    REAL temp3=0., temp4 = 0., div=0.;
    temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    r = sqrt(temp1);
    temp2 = (2./alpha) + 2.*r0*(r0-r)*alpha;
    temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
    temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
    div = temp2/(temp3*temp4);
    gradu(3,0) = div; //valor do divergente
    
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

void Solution(const TPZVec<REAL> &p, TPZVec<STATE> &f) {   //Jorge 2017    It is not used:, TPZFMatrix<STATE> &gradf){

    REAL x,y,z;
    x = p[0];
    y = p[1];
    z = p[2];
    
//    REAL r = sqrt(x*x+y*y+z*z);
    
#ifdef Solution1
    
    f[0] = r*r;
    
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
    
    
//    REAL denomfactor1 = -9.0 + d*d*(-9.0*rad+M_PI*(-M_PI+6.0*sqrt_rad));
//    REAL numfactro1   = 18.0*d*(-9.0+d*d*M_PI*(-M_PI+3.0*sqrt_rad));
    
    f[0] = piover2 - artan_arg;
    
#endif
    
#ifdef Solution7
    
    REAL sin_pi_x = sin(M_PI*x);
    REAL sin_pi_y = sin(M_PI*y);
    REAL sin_pi_z = sin(M_PI*z);
    
    f[0] = sin_pi_x * sin_pi_y * sin_pi_z;
    
#endif
    
    
#ifdef Solution8 //arctan solution
    
    REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
    REAL r0 = M_PI/3.;
    REAL alpha = 5.;
    
    REAL temp1;
    REAL temp2;
    REAL r;
    REAL grad;
    temp1 = 0.,temp2=0., r=0., grad=0.;
    
    //solucao u
    temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    r = sqrt(temp1);
    temp2 = (r - r0)*alpha;
    f[0] = M_PI/2. - atan(temp2);
    
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
    
//    REAL r = sqrt(x*x+y*y+z*z);
    
#ifdef Solution1
    
    f[0] = -6.0;
    
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
    

#ifdef Solution7
    
    REAL sin_pi_x = sin(M_PI*x);
    REAL sin_pi_y = sin(M_PI*y);
    REAL sin_pi_z = sin(M_PI*z);
    
    f[0] = 3.0 * M_PI * M_PI * sin_pi_x * sin_pi_y * sin_pi_z;
    
#endif
    
    
#ifdef Solution8 // arctan solution
    
    f[0] = 0.;
    REAL x0 = 1.25, y0 = -0.25, z0 = -0.25;
    REAL r0 = M_PI/3.;
    REAL alpha = 5.;
    
    REAL temp1=0., temp2=0.,temp3=0., temp4=0., r=0., sol=0.;
    
    temp1 = (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0);
    r = sqrt(temp1);
    
    temp2 = (2./alpha) + 2.*r0*(r0-r)*alpha;
    temp3 = r*((r - r0)*(r - r0) + 1./(alpha*alpha));
    temp4 = 1. + (r - r0)*(r - r0)*alpha*alpha;
    
    sol = temp2/(temp3*temp4);
    f[0] = sol;
    
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
        
        TPZManVector<REAL,3> triplet(3,0.0);
        REAL w;
        
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
    
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<STATE,10> globalerror(3,0.   );
    for (int64_t iel = 0; iel < nel; iel++) {
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
        
        TPZManVector<REAL,3> triplet(3,0.0);
        REAL w;
        
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
    
    if (sim_data.domain_type == "cube") {
        
        if (sim_data.mesh_type == "linear") {
            
            switch (sim_data.elemen_type) {
                case 0:
                    { // T
                        geometry = MakeCubeFromTetrahedrons(ndiv, sim_data);
                    }
                    break;
                case 1:
                    { // H
                        geometry = MakeCubeFromHexahedrons(ndiv, sim_data);
                    }
                    break;
                case 2:
                    { // P
                        geometry = MakeCubeFromPrisms(ndiv, sim_data);
                    }
                    break;
                default:
                {
                    std::cout << "Element type not implemented. " << std::endl;
                    DebugStop();
                }
                    break;
            }
            return geometry;
        }

        std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
        DebugStop();
        
    }
    
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
            
            if (sim_data.UseGmshMeshQ) {
                
                std::string dirname = PZSOURCEDIR;
                std::string grid;
                
                switch (sim_data.elemen_type) {
                    case 0:
                    { // T
                        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/vertical_wellbore_Te_l_" + std::to_string(ndiv) +".msh";
                    }
                        break;
                    case 1:
                    { // H
                        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/vertical_wellbore_He_l_" + std::to_string(ndiv) +".msh";
                    }
                        break;
                    case 2:
                    { // P
                        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/vertical_wellbore_Pe_l_" + std::to_string(ndiv) +".msh";
                    }
                        break;
                    case 3:
                    { // P
                        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/vertical_wellbore_hybrid_l_" + std::to_string(ndiv) +".msh";
                    }
                        break;
                    default:
                    {
                        std::cout << "Element type not implemented. " << std::endl;
                        DebugStop();
                    }
                        break;
                }
                
                TPZGmshReader Geometry;
                REAL s = 1.0;
                Geometry.SetfDimensionlessL(s);
                geometry = Geometry.GeometricGmshMesh(grid);
                const std::string name("Geometry and mesh from gmsh script");
                geometry->SetName(name);
                
                // changin id internally
                sim_data.gamma_ids[0] = 2;
                sim_data.gamma_ids[1] = 3;
                
            }
            else{
                geometry = MakeCylinderFromLinearFaces(ndiv, sim_data);
            }
            
            return geometry;
        }
        
        std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
        DebugStop();
        
    }
    
    std::cout << "Error:: unable to compute geometry for the given case = " << &sim_data << std::endl;
    DebugStop();
    return geometry;
}

TPZCompMesh * ComputationalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, int64_t &ndof, TPZVec<TPZCompMesh *> &meshvec){
    
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
        step.SetDirect(ECholesky);
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

    int div = 0;
    an->DefineGraphMesh(dim,scalnames,vecnames,file);
    an->PostProcess(div,dim);
    
}

TPZCompMesh * PrimalMesh(TPZGeoMesh * geometry, int p, SimulationCase sim_data, int64_t &ndof){
    
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
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic, 5);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution, 5);
            analytic_bc->SetPolynomialOrder(sim_data.int_order);
            TPZAutoPointer< TPZFunction<STATE> > solution = analytic_bc;
            
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            face->SetForcingFunction(solution);
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
    
    TPZCompMeshTools::CreatedCondensedElements(cmesh, false);
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Primal_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
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
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic, 5);
        analytic->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > analytic_full = analytic;
        volume->SetForcingFunctionExact(analytic_full);
        
        cmesh->InsertMaterialObject(volume);
        
        for (int ib = 0; ib < nboundaries; ib++) {
            
            TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution, 5);
            analytic_bc->SetPolynomialOrder(sim_data.int_order);
            TPZAutoPointer< TPZFunction<STATE> > solution = analytic_bc;
            
            TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ib],dirichlet,val1,val2);
            face->SetForcingFunction(solution);
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
    
    std::cout << "Created multi physics mesh\n";
    if (sim_data.IsMHMQ) {
//        BuildMacroElements(cmesh);
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    else{
        TPZCompMeshTools::GroupElements(cmesh);
        std::cout << "Created grouped elements\n";
        bool keepmatrix = false;
        bool keeponelagrangian = true;
        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
        std::cout << "Created condensed elements\n";
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }

#ifdef PZDEBUG2
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    meshvec = meshvector;
    
    return cmesh;
    

}

/// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus)
{
//    int dim = fluxmesh->Dimension();
//    /// loop over all the elements
//    int64_t nel = fluxmesh->NElements();
//    for (int64_t el=0; el<nel; el++) {
//        TPZCompEl *cel = fluxmesh->Element(el);
//        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
//        if (!intel) {
//            continue;
//        }
//        TPZGeoEl *gel = intel->Reference();
//        if (gel->Dimension() != dim) {
//            continue;
//        }
//        // compute the maxorder
//        int maxorder = -1;
//        int ncon = intel->NConnects();
//        for (int i=0; i<ncon-1; i++) {
//            int conorder = intel->Connect(i).Order();
//            maxorder = maxorder < conorder ? conorder : maxorder;
//        }
//        int nsides = gel->NSides();
//        int nconside = intel->NSideConnects(nsides-1);
//        // tive que tirar para rodar H1
//        //        if (nconside != 1 || maxorder == -1) {
//        //            DebugStop();
//        //        }
//        int64_t cindex = intel->SideConnectIndex(nconside-1, nsides-1);
//        TPZConnect &c = fluxmesh->ConnectVec()[cindex];
//        if (c.NElConnected() != 1) {
//            DebugStop();
//        }
//        if (c.Order()+hdivplusplus != maxorder) {
//            //            std::cout << "Changing the order of the central connect " << cindex << " from " << c.Order() << " to " << maxorder+hdivplusplus << std::endl;
//            // change the internal connect order to be equal do maxorder
//            intel->SetSideOrder(nsides-1, maxorder+hdivplusplus);
//        }
//    }
//    fluxmesh->ExpandSolution();
    
    int nEl= fluxmesh-> NElements();
    int dim = fluxmesh->Dimension();

    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = fluxmesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;

        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();

            int neworder = corder + hdivplusplus;//Aqui = +1
            int64_t cindex = cel->ConnectIndex(ncon-1);
            conel.SetOrder(neworder,cindex);

            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape = intel->NConnectShapeF(ncon-1,neworder);

            conel.SetNShape(nshape);
            fluxmesh->Block().Set(conel.SequenceNumber(),nshape);
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
    TPZManVector<int64_t> pressorder(pressuremesh->NElements(),-1);
    int64_t nel = fluxmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
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
        int64_t cindex = intel->SideConnectIndex(0, nsides-1);
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
    for (int64_t el=0; el<nel; el++) {
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
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution, 5);
        analytic_bc->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > solution = analytic_bc;
        volume->SetForcingFunction(solution);
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic, 5);
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
    std::ofstream sout(file_name.str().c_str());
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
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(f, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
        
        TPZDummyFunction<STATE> * analytic_bc = new TPZDummyFunction<STATE>(Solution, 5);
        analytic_bc->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > solution = analytic_bc;
        volume->SetForcingFunction(solution);
        
        
        TPZDummyFunction<STATE> * analytic = new TPZDummyFunction<STATE>(Analytic, 5);
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
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}

TPZGeoMesh * MakeCubeFromLinearQuadrilateralFaces(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    REAL t=0.0;
    int n = pow(2,ndiv+1);
    REAL dt = 1.0/REAL(n);

    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh_point = new TPZGeoMesh;
    GeoMesh_point->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh_point->NodeVec()[0]=Node;
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    int front = -1;
    int back  = -2;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh_point);
    GeoMesh_point->BuildConnectivity();
    GeoMesh_point->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh_point);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction_x, 5);
    CreateGridFrom.SetParametricFunction(ParFunc);
    CreateGridFrom.SetFrontBackMatId(front, front);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh_line = CreateGridFrom.ComputeExtrusion(t, dt, n);

    TPZHierarquicalGrid CreateGridFrom2(GeoMesh_line);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction_y, 5);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(front, front);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_surface = CreateGridFrom2.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom3(GeoMesh_surface);
    GeoMesh_surface->SetDimension(2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction_z, 5);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    CreateGridFrom3.SetFrontBackMatId(back, back);
    
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_cube = CreateGridFrom3.ComputeExtrusion(t, dt, n);

//    switch (sim_data.perturbation_type) {
//        case 1:
//        {
//            PerturbateNodes_I(GeoMesh_cube, ndiv);
//        }
//            break;
//        default:
//            std::cout << "The mesh is not perturbed" << std::endl;
//            break;
//    }
    
    return GeoMesh_cube;
}

TPZGeoMesh * MakeCubeFromLinearTriangularFaces(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    

    
    TPZGeoMesh * GeoMesh_cube = new TPZGeoMesh;
    
    if (sim_data.UseGmshMeshQ) {
        
        std::string dirname = PZSOURCEDIR;
        std::string grid;
        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryH_adapted.msh";
        
        TPZGmshReader Geometry;
        REAL s = 1.0;
        Geometry.SetfDimensionlessL(s);
        GeoMesh_cube = Geometry.GeometricGmshMesh(grid);
        const std::string name("Geometry and mesh from gmsh script");
        GeoMesh_cube->SetName(name);
        
        // changin id internally
        sim_data.gamma_ids[0] = 2;
        sim_data.gamma_ids[1] = 3;
    
    }
    else
    {
    
        REAL t=0.0;
        int n = pow(2,ndiv);
        REAL dt = 1.0/REAL(n);
        
        
        // Creating a 0D element to be extruded
        TPZGeoMesh * GeoMesh_point = new TPZGeoMesh;
        GeoMesh_point->NodeVec().Resize(1);
        TPZGeoNode Node;
        TPZVec<REAL> coors(3,0.0);
        Node.SetCoord(coors);
        Node.SetNodeId(0);
        GeoMesh_point->NodeVec()[0]=Node;
        
        TPZVec<int64_t> Topology(1,0);
        int elid=0;
        int matid=1;
        int front = -1;
        int back  = -2;
        
        new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh_point);
        GeoMesh_point->BuildConnectivity();
        GeoMesh_point->SetDimension(0);
        
        TPZHierarquicalGrid CreateGridFrom(GeoMesh_point);
        TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction_x, 5);
        CreateGridFrom.SetParametricFunction(ParFunc);
        CreateGridFrom.SetFrontBackMatId(front, front);
        
        // Computing Mesh extruded along the parametric curve Parametricfunction
        TPZGeoMesh * GeoMesh_line = CreateGridFrom.ComputeExtrusion(t, dt, n);
        
        TPZHierarquicalGrid CreateGridFrom2(GeoMesh_line);
        TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction_y, 5);
        CreateGridFrom2.SetParametricFunction(ParFunc2);
        CreateGridFrom2.SetFrontBackMatId(front, front);
        CreateGridFrom2.SetTriangleExtrusion();
        
        // Computing Mesh extruded along the parametric curve Parametricfunction2
        TPZGeoMesh * GeoMesh_surface = CreateGridFrom2.ComputeExtrusion(t, dt, n);
        
        TPZHierarquicalGrid CreateGridFrom3(GeoMesh_surface);
        GeoMesh_surface->SetDimension(2);
        TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction_z, 5);
        CreateGridFrom3.SetParametricFunction(ParFunc3);
        CreateGridFrom3.SetFrontBackMatId(back, back);
        CreateGridFrom3.SetTriangleExtrusion();
        CreateGridFrom3.SetTetrahedonExtrusion();
        
        // Computing Mesh extruded along the parametric curve Parametricfunction2
        TPZGeoMesh * GeoMesh_cube = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    
    }

    UniformRefinement(GeoMesh_cube, ndiv);
    RefineHexahedronsToTetrahedrons(GeoMesh_cube, 1);
    
    switch (sim_data.perturbation_type) {
        case 1:
        {
            RefineTetrahedronsToHexahedrons(GeoMesh_cube, 1); // Divide each tetrahedron into 4 hexahedra.
//            UniformRefinement(GeoMesh_cube, ndiv);
        }
            break;
        default:
            std::cout << "The mesh is not perturbed" << std::endl;
            break;
    }
    
    
    return GeoMesh_cube;
}

TPZGeoMesh * MakeCube(SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif

    REAL t=0.0;
    int n = pow(2,0);
    REAL dt = 1.0/REAL(n);

    // Create the Base Hexahedron

    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh_point = new TPZGeoMesh;
    GeoMesh_point->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh_point->NodeVec()[0]=Node;

    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    int front = -1;
    int back  = -2;

    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh_point);
    GeoMesh_point->BuildConnectivity();
    GeoMesh_point->SetDimension(0);

    TPZHierarquicalGrid CreateGridFrom(GeoMesh_point);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc = new TPZDummyFunction<REAL>(Parametricfunction_x, 5);
    CreateGridFrom.SetParametricFunction(ParFunc);
    CreateGridFrom.SetFrontBackMatId(front, front);

    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh_line = CreateGridFrom.ComputeExtrusion(t, dt, n);

    TPZHierarquicalGrid CreateGridFrom2(GeoMesh_line);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc2 = new TPZDummyFunction<REAL>(Parametricfunction_y, 5);
    CreateGridFrom2.SetParametricFunction(ParFunc2);
    CreateGridFrom2.SetFrontBackMatId(front, front);

    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_surface = CreateGridFrom2.ComputeExtrusion(t, dt, n);

    TPZHierarquicalGrid CreateGridFrom3(GeoMesh_surface);
    GeoMesh_surface->SetDimension(2);
    TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction_z, 5);
    CreateGridFrom3.SetParametricFunction(ParFunc3);
    CreateGridFrom3.SetFrontBackMatId(back, back);

    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh_cube = CreateGridFrom3.ComputeExtrusion(t, dt, n);
    
    
//#ifdef PZDEBUG
//    //  Print Geometrical Base Mesh
//    std::ofstream argument("TheCube.txt");
//    GeoMesh_cube->Print(argument);
//    std::ofstream Dummyfile("TheCube.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh_cube,Dummyfile, true);
//
//#endif
    
    return GeoMesh_cube;
}

TPZGeoMesh * MakeCubeFromTetrahedrons(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * GeoMesh_cube = MakeCube(sim_data);
    UniformRefinement(GeoMesh_cube, ndiv+1);
    RefineHexahedronsToTetrahedrons(GeoMesh_cube, 1);
    return GeoMesh_cube;
}

TPZGeoMesh * MakeCubeFromHexahedrons(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * GeoMesh_cube = MakeCube(sim_data);
    
    if (sim_data.NonAffineQ) {
        UniformRefinement(GeoMesh_cube, ndiv);
        RefineHexahedronsToTetrahedrons(GeoMesh_cube, 1);
        RefineTetrahedronsToHexahedrons(GeoMesh_cube, 1);
    }else{
        UniformRefinement(GeoMesh_cube, ndiv+1);
    }

//#ifdef PZDEBUG
//    //  Print Geometrical Base Mesh
//    std::ofstream argument("TheCube.txt");
//    GeoMesh_cube->Print(argument);
//    std::ofstream Dummyfile("TheCube.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh_cube,Dummyfile, true);
//
//#endif
    
    return GeoMesh_cube;
}

TPZGeoMesh * MakeCubeFromPrisms(int ndiv, SimulationCase  & sim_data){
    
#ifdef PZDEBUG
    if (sim_data.omega_ids.size() != 1 || sim_data.gamma_ids.size() != 2) {
        std::cout << "Cube:: Please pass materials ids for volume (1) and boundary domains = { (-1 -> first surface ),  (-2 -> second surface) }" << std::endl;
        DebugStop();
    }
#endif
    
    TPZGeoMesh * GeoMesh = new TPZGeoMesh;
    
    if (sim_data.NonAffineQ) {
        
        GeoMesh = MakeCube(sim_data);
        UniformRefinement(GeoMesh, ndiv+1);
//        RefineHexahedronsToTetrahedrons(GeoMesh, 1);
//        RefineTetrahedronsToHexahedrons(GeoMesh, 1);
//        RefineHexahedronsToPrisms(GeoMesh, 1);
        
        RefineHexahedronsToNonAffinePrisms(GeoMesh, 1);
        
//        std::string dirname = PZSOURCEDIR;
//        std::string grid;
//        grid = dirname + "/Projects/Hdiv3DCurved/gmsh_meshes/msh/GeometryDisc.msh";
//
//        TPZGmshReader Geometry;
//        REAL s = 1.0;
//        Geometry.SetfDimensionlessL(s);
//        TPZGeoMesh * GeoMesh2D = Geometry.GeometricGmshMesh(grid);
//
//        REAL t=0.0;
//        int n = pow(2,2);
//        REAL dt = 1.0/REAL(n);
//
//        int front = -1;
//        int back  = -2;
//
//        // Create the Prism extrusion
//        TPZHierarquicalGrid CreateGridFromTriangles(GeoMesh2D);
//        GeoMesh2D->SetDimension(2);
//        TPZAutoPointer<TPZFunction<REAL> > ParFunc3 = new TPZDummyFunction<REAL>(Parametricfunction_z);
//        CreateGridFromTriangles.SetParametricFunction(ParFunc3);
//        CreateGridFromTriangles.SetFrontBackMatId(back, back);
//        CreateGridFromTriangles.SetPrismExtrusion();
//
//        // Computing Mesh extruded along the parametric curve Parametricfunction2
//        GeoMesh = CreateGridFromTriangles.ComputeExtrusion(t, dt, n);
//        UniformRefinement(GeoMesh, ndiv);
////        RefinePrismsToNonAffinePrisms(GeoMesh_cube, 1);
    }else{
        GeoMesh = MakeCube(sim_data);
        RefineHexahedronsToPrisms(GeoMesh, 1);
        UniformRefinement(GeoMesh, ndiv+1);
    }
    
//#ifdef PZDEBUG
//    //  Print Geometrical Base Mesh
//    std::ofstream argument("TheCylinder.txt");
//    GeoMesh->Print(argument);
//    std::ofstream Dummyfile("TheCylinder.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(GeoMesh,Dummyfile, true);
//
//#endif
    
    return GeoMesh;
}

void Parametricfunction_x(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0;
}

void Parametricfunction_y(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void Parametricfunction_z(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
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
        const std::string name("Geometry and mesh from gmsh script");
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
        
        TPZManVector<int64_t,4> TopolQuad(4);
        TPZManVector<int64_t,8> TopolCube(8);
        TPZManVector<REAL,3> coord(3,0.);
        TPZVec<REAL> xc(3,0.);
        REAL cphi = atan(sqrt(2.0));
        
        int nodeindex = 0;
        int64_t id = 0;
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
    
    TPZManVector<int64_t,4> TopolQuad(4);
    TPZManVector<int64_t,8> TopolCube(8);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    int64_t id = 0;
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
        
        TPZManVector<int64_t,1> TopolPoint(1);
        TPZManVector<int64_t,4> TopolQuad(4),TopolArc(3);
        TPZManVector<int64_t,8> TopolCube(8),TopolQuadQua(8);
        TPZManVector<REAL,3> coord(3,0.);
        TPZVec<REAL> xc(3,0.);
        REAL cphi = atan(sqrt(2.0));
        
        int nodeindex = 0;
        int64_t id = 0;
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
    
    int64_t last_node = geomesh->NNodes() - 1;
    int64_t last_element = geomesh->NElements() - 1;
    int64_t node_id = geomesh->NodeVec()[last_node].Id();
    int64_t element_id = geomesh->Element(last_element)->Id();
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
    TPZAutoPointer<TPZFunction<REAL> > ParFuncZ = new TPZDummyFunction<REAL>(ParametricfunctionZ, 5);
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
    
    int64_t last_node = gmesh_3d->NNodes() - 1;
    int64_t last_element = gmesh_3d->NElements() - 1;
    int64_t node_id = gmesh_3d->NodeVec()[last_node].Id();
    int64_t element_id = gmesh_3d->Element(last_element)->Id();
    const std::string name("Reservoir with vertical extrusion");
    gmesh_3d->SetName(name);
    gmesh_3d->SetMaxNodeId(node_id);
    gmesh_3d->SetMaxElementId(element_id);
    gmesh_3d->SetDimension(3);
    return gmesh_3d;
    
}

void ParametricfunctionZ(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}


void TransformToQuadratic(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++)
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


TPZManVector<REAL,3> ParametricSphere(REAL radius, REAL theta, REAL phi)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * sin(theta) * cos(phi) ;
    xcoor[1] = radius * sin(theta) * sin(phi) ;
    xcoor[2] = radius * cos(theta) ;
    return xcoor;
}

TPZManVector<REAL,3> ParametricCylinder(REAL radius, REAL theta, REAL z)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
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

    std::ofstream textfile(text_name.str().c_str());
    gmesh->Print(textfile);
    std::ofstream vtkfile(vtk_name.str().c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
}

void PrintGeometryVols(TPZGeoMesh * gmesh, std::stringstream & file_name){
    
    std::ofstream mesh_data(file_name.str().c_str(),std::ios::app);
    
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    TPZStack<int64_t> gel_indexes;
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZGeoEl * gel = gmesh->Element(iel);
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->Dimension() != dim || gel->HasSubElement()) {
            continue;
        }
        gel_indexes.Push(gel->Index());
    }
    
    // Writing element nodes
    int64_t n_vols = gel_indexes.size();
    mesh_data << std::setw(5) << n_vols << std::endl;
    TPZFMatrix<REAL> xcoor;
    
    for (int64_t ivol = 0; ivol < n_vols; ivol++) {
        TPZGeoEl * gel = gmesh->Element(gel_indexes[ivol]);
        gel->NodesCoordinates(xcoor);
        int nr = xcoor.Rows();
        int nc = xcoor.Cols();
        for (int c = 0; c < nc; c++) {
            for (int r = 0; r < nr; r++) {
                mesh_data << std::setw(10) << xcoor(r,c) ;
            }
        }
        mesh_data << std::endl;
    }
    mesh_data.flush();
    return;
}

/** @brief Apply uniform refinement on the Geometric mesh */
void UniformRefinement(TPZGeoMesh * gmesh, int n_ref){
    for ( int ref = 0; ref < n_ref; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gmesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
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

void RefineTetrahedronsToHexahedrons(TPZGeoMesh *gmesh, int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D, refp2D;
    
    {// needed!
        char buf[] =
        "15 5 "
        "-60 TetToHex "
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
        "0.333333333333333 0.333333333333333 0. "
        "0. 0.333333333333333 0.333333333333333 "
        "0.333333333333333 0. 0.333333333333333 "
        "0.333333333333333 0.333333333333333 0.333333333333333 "
        "0.25 0.25 0.25 "
        "4 4 0  1  2  3 "
        "7 8 0 4 10 5 6 12 14 11 "
        "7 8 4 1 7 10 12 9 13 14 "
        "7 8 5 10 7 2 11 14 13 8 "
        "7 8 6 12 14 11 3 9 13 8 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    {// needed!
        char buf[] =
        "7 4 "
        "-60 TriToQuad "
        "0 0 0 "
        "1 0 0 "
        "0 1 0 "
        "0.5 0 0 "
        "0 0.5 0 "
        "0.5 0.5 0 "
        "0.333333333333333 0.333333333333333 0. "
        "2 3 0  1  2 "
        "3 4 0 3 6 4 "
        "3 4 3 1 5 6 "
        "3 4 4 6 5 2 ";
        std::istringstream str(buf);
        refp2D = new TPZRefPattern(str);
        refp2D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp2D);
        if(!refp2D)
        {
            DebugStop();
        }
    }
    
    if(!refp3D || !refp2D)
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
                gel->SetRefPattern(refp2D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            
        }
    }
}

void RefineHexahedronsToTetrahedrons(TPZGeoMesh *gmesh, int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D, refp2D;
    
    {// needed!
        char buf[] =
        "8 7 "
        "-70 HexToTet "
        "-1 -1 -1 "
        "1 -1 -1 "
        "1 1 -1 "
        "-1 1 -1 "
        "-1 -1 1 "
        "1 -1 1 "
        "1 1 1 "
        "-1 1 1 "
        "7 8 0 1 2 3 4 5 6 7 "
        "4 4 0 1 3 4 "
        "4 4 4 1 3 5 "
        "4 4 1 2 3 5 "
        "4 4 5 7 2 6 "
        "4 4 5 2 3 7 "
        "4 4 4 5 3 7 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    {// needed!
        char buf[] =
        "4 3 "
        "-60 QuadTo "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0 "
        "-1 1 0 "
        "3 4 0 1 2 3 "
        "2 3 0 1 3 "
        "2 3 1 2 3 ";
        std::istringstream str(buf);
        refp2D = new TPZRefPattern(str);
        refp2D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp2D);
        if(!refp2D)
        {
            DebugStop();
        }
    }
    
    if(!refp3D || !refp2D)
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
                gel->SetRefPattern(refp2D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            
        }
    }
}

void RefineHexahedronsToNonAffineHexahedrons(TPZGeoMesh *gmesh, int n_ref){
    
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
    DebugStop();
}

void RefineHexahedronsToPrisms(TPZGeoMesh *gmesh, int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D, refp2D_T;
    
    {// needed!
        char buf[] =
        "8 3 "
        "-70 HexToTet "
        "-1 -1 -1 "
        "1 -1 -1 "
        "1 1 -1 "
        "-1 1 -1 "
        "-1 -1 1 "
        "1 -1 1 "
        "1 1 1 "
        "-1 1 1 "
        "7 8 0 1 2 3 4 5 6 7 "
        "6 6 0 1 3 4 5 7 "
        "6 6 1 2 3 5 6 7 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    {// needed!
        char buf[] =
        "4 3 "
        "-60 QuadTo "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0 "
        "-1 1 0 "
        "3 4 0 1 2 3 "
        "2 3 0 1 3 "
        "2 3 1 2 3 ";
        std::istringstream str(buf);
        refp2D_T = new TPZRefPattern(str);
        refp2D_T->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp2D_T);
        if(!refp2D_T)
        {
            DebugStop();
        }
    }
    
    if(!refp3D || !refp2D_T)
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
        }
    }
    
    for(int r = 0; r < n_ref; r++)
    {
        int nels = gmesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==3)
            {
                continue;
            }
            
            TPZStack<TPZGeoElSide> elsides;
            if(gel->Dimension()==2)
            {

                if (gel->MaterialId() == -2) {
                    gel->SetRefPattern(refp2D_T);
                    TPZVec<TPZGeoEl*> sons;
                    gel->Divide(sons);
                }
            }
            
        }
    }
    
}

void RefineHexahedronsToNonAffinePrisms(TPZGeoMesh *gmesh, int n_ref){
    
    TPZAutoPointer<TPZRefPattern> refp3D, refp2D_T, refp2D_Q;
    
    {// needed!
        char buf[] =
        "27 17 "
        "-70 HexToPrism "
        "-1 -1 -1 "
        "1 -1 -1 "
        "1 1 -1 "
        "-1 1 -1 "
        "-1 -1 1 "
        "1 -1 1 "
        "1 1 1 "
        "-1 1 1 "
        "0 -1 -1 "
        "1 0 -1 "
        "0 1 -1 "
        "-1 0 -1 "
        "0 0 -1 "
        "0 -1 0.5 "
        "1 0 0.5 "
        "0 1 0.5 "
        "-1 0 0.5 "
        "0 0 0 "
        "0 -1 1 "
        "1 0 1 "
        "0 1 1 "
        "-1 0 1 "
        "0 0 1 "
        "-1 -1 -0.5 "
        "1 -1 -0.5 "
        "1 1 -0.5 "
        "-1 1 -0.5 "
        "7 8 0 1 2 3 4 5 6 7 "
        "6 6 1 9 12 24 14 17 "
        "6 6 1 12 8 24 17 13 "
        "6 6 2 9 12 25 14 17 "
        "6 6 2 12 10 25 17 15 "
        "6 6 3 11 12 26 16 17 "
        "6 6 3 12 10 26 17 15 "
        "6 6 0 11 12 23 16 17 "
        "6 6 0 12 8 23 17 13 "
        "6 6 24 14 17 5 19 22 "
        "6 6 24 17 13 5 22 18 "
        "6 6 25 14 17 6 19 22 "
        "6 6 25 17 15 6 22 20 "
        "6 6 26 16 17 7 21 22 "
        "6 6 26 17 15 7 22 20 "
        "6 6 23 16 17 4 21 22 "
        "6 6 23 17 13 4 22 18 ";
        std::istringstream str(buf);
        refp3D = new TPZRefPattern(str);
        refp3D->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp3D);
        if(!refp3D)
        {
            DebugStop();
        }
    }
    
    {// needed!
        char buf[] =
        "9 9 "
        "-60 QuadTo8Tri "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0 "
        "-1 1 0 "
        "0 -1 0 "
        "1 0 0 "
        "0 1 0 "
        "-1 0 0 "
        "0 0 0 "
        "3 4 0 1 2 3 "
        "2 3 0 8 7 "
        "2 3 0 4 8 "
        "2 3 4 1 8 "
        "2 3 1 5 8 "
        "2 3 8 5 2 "
        "2 3 8 2 6 "
        "2 3 6 3 8 "
        "2 3 3 7 8 ";
        std::istringstream str(buf);
        refp2D_T = new TPZRefPattern(str);
        refp2D_T->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp2D_T);
        if(!refp2D_T)
        {
            DebugStop();
        }
    }
    
    {// needed!
        char buf[] =
        "9 5 "
        "-60 QuadTo4Quad "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0 "
        "-1 1 0 "
        "0 -1 0 "
        "1 -0.5 0 "
        "0 1 0 "
        "-1 -0.5 0 "
        "0 0.5 0 "
        "3 4 0 1 2 3 "
        "3 4 0 4 8 7 "
        "3 4 4 1 5 8 "
        "3 4 5 2 6 8 "
        "3 4 6 3 7 8 ";
        std::istringstream str(buf);
        refp2D_Q = new TPZRefPattern(str);
        refp2D_Q->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp2D_Q);
        if(!refp2D_Q)
        {
            DebugStop();
        }
    }
    
    if(!refp3D || !refp2D_T || !refp2D_Q)
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
        }
    }
    
    for(int r = 0; r < n_ref; r++)
    {
        int nels = gmesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==3)
            {
                continue;
            }
            
            TPZStack<TPZGeoElSide> elsides;
            if(gel->Dimension()==2)
            {
                
                if (gel->MaterialId() == -2) {
                    gel->SetRefPattern(refp2D_T);
                }else{
                    gel->SetRefPattern(refp2D_Q);
                }
                
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
    
    TPZVec<REAL> iCoords(3,0.0);
    TPZVec<REAL> iCoordsRotated(3,0.0);
    
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

void PerturbateNodes_I(TPZGeoMesh *gmesh, int ndiv)
{
    TPZManVector<REAL> iCoords(3,0.);
    TPZManVector<REAL> iCoordsRotated(3,0.);
    TPZStack<int64_t> inter_nodes;
    TPZStack<int64_t> exter_nodes;
    int NumberofGeoNodes = gmesh->NNodes();
    REAL pertub_param = (ndiv+2)*10;
    
    //maior comprimento da malha
    REAL Lx = 1.;
    REAL hsize = Lx;//Lx/ndiv;
    
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
        
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        REAL iCoords0 = iCoords[0];
        REAL iCoords1 = iCoords[1];
        REAL iCoords2 = iCoords[2];
        
        REAL maxcoord = Max(iCoords0, iCoords1);
        maxcoord = Max(maxcoord, iCoords2);
        REAL mincoord = Min(iCoords0, iCoords1);
        mincoord = Min(mincoord, iCoords2);
        
        if(mincoord==0. || maxcoord==Lx)
        {
            exter_nodes.Push(inode);
        }else
        {
            inter_nodes.Push(inode);
            hsize = Lx*pow(-1.,inode);
            PertubationMatrix_I(iCoords, pertub_param, hsize, iCoordsRotated);
            
            GeoNode.SetCoord(iCoordsRotated);
            gmesh->NodeVec()[inode] = GeoNode;
        }
    }
    int nintnodes = inter_nodes.size();
    int nextnodes = exter_nodes.size();
    if(nintnodes+nextnodes!= NumberofGeoNodes){
        DebugStop();
    }
    inter_nodes.Print();
    exter_nodes.Print();
}

void PertubationMatrix_I(TPZManVector<REAL> CoordX, REAL pertub_param, REAL ElSize, TPZManVector<REAL> &CoordsPertubated)
{
    //Rotation matrix: Rotation_x*Rotation_y*Rotation_z
    TPZFMatrix<REAL> matM(3,3,0.);
    
    //Dilation or contraction matrix
    TPZFMatrix<REAL> matM_alpha(3,3,0.);
    
    //Rotation angle
    REAL teta1 = CoordX[0];
    REAL teta2 = CoordX[1];
    REAL teta3 = CoordX[2];
    REAL sin1 = sin(teta1), sin2 = sin(teta2), sin3 = sin(teta3);
    REAL cos1 = cos(teta1), cos2 = cos(teta2), cos3 = cos(teta3);
    
    //rotation matrix
    matM(0,0) = cos2*cos3;
    matM(0,1) = -cos2*sin3;
    matM(0,2) = sin2;
    matM(1,0) = cos3*sin1*sin2 + cos1*sin3;
    matM(1,1) = cos1*cos3 - sin1*sin2*sin3;
    matM(1,2) = -cos2*sin1;
    matM(2,0) = -cos1*cos3*sin2 + sin1*sin3;
    matM(2,1) = cos3*sin1 + cos1*sin2*sin3;
    matM(2,2) = cos1*cos2;
    
    //Dilation or contraction matrix
    REAL temp1 = CoordX[0]*CoordX[0] + CoordX[1]*CoordX[1] + CoordX[2]*CoordX[2];
    REAL normacoord = sqrt(temp1);
    REAL alpha = (ElSize/pertub_param)*sin(normacoord);
    matM_alpha(0,0) = alpha, matM_alpha(1,1) = alpha, matM_alpha(2,2) = alpha;
    
    //coordinate pertubated
    REAL temp2 = sin1*sin1 + sin2*sin2 + sin3*sin3;
    //REAL temp2 = cos1*cos1 + cos2*cos2 + cos3*cos3;
    TPZFMatrix<REAL> uniVec(3,1,0.);
    uniVec(0,0) = sin1/temp2, uniVec(1,0) = sin2/temp2, uniVec(2,0) = sin3/temp2;
    //uniVec(0,0) = cos1/temp2, uniVec(1,0) = cos2/temp2, uniVec(2,0) = cos3/temp2;
    TPZFMatrix<REAL> matRes1(3,1,0.);
    
   // matRes1.Print("matres0=");
    //matM*uniVec
    matM.MultAdd(uniVec, matM, matRes1);
   // matRes1.Print("matres1=");
    
    //matM_alpha*(matM*uniVec)
    TPZFMatrix<REAL> matRes2(3,1,0.);
    matM_alpha.MultAdd(matRes1, matM_alpha, matRes2);
    //matM_alpha.Print("matalpha=");
    //matRes2.Print("matres2=");
    
    //matM_alpha*(matM*uniVec) + coordX
    CoordsPertubated[0] = matRes2(0,0) + CoordX[0];
    CoordsPertubated[1] = matRes2(1,0) + CoordX[1];
    CoordsPertubated[2] = matRes2(2,0) + CoordX[2];
 }

void ErrorH1(TPZAnalysis * analysis, REAL &error_primal , REAL & error_dual, REAL & error_h1)
{
    bool Serial_ErrorQ = false;
    int nthreads = 12;
    
    TPZCompMesh * cmesh = analysis->Mesh();
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<REAL,10> globalerror(3,0.);
    
    if (Serial_ErrorQ) {
        for (int64_t iel = 0; iel < nel; iel++) {
            TPZCompEl *cel = cmesh->ElementVec()[iel];
            
            if (!cel) {
                continue;
            }
            
            if(cel->Reference()->Dimension()!=dim) {
                continue;
            }
            
            TPZManVector<REAL,10> elerror(3,0.);
            elerror.Fill(0.);
            cel->EvaluateError(Analytic, elerror, 0);
            int nerr = elerror.size();
            for (int i=0; i < nerr; i++) {
                globalerror[i] += elerror[i]*elerror[i];
                
            }
        }
        error_primal    = sqrt(globalerror[0]);
        error_dual      = sqrt(globalerror[1]);
        error_h1        = error_primal + error_dual;
    }
    else{
        analysis->SetThreadsForError(nthreads);
        analysis->SetExact(Analytic);
        analysis->PostProcessError(globalerror);
        
        error_primal    = globalerror[0];
        error_dual      = globalerror[1];
        error_h1        = error_primal + error_dual;
    }
    

    
    
}

void ErrorHdiv(TPZAnalysis * analysis, REAL &error_primal , REAL & error_dual, REAL & error_hdiv){
    
    bool Serial_ErrorQ = false;
    int nthreads = 12;
    
    TPZCompMesh * cmesh = analysis->Mesh();
    int64_t nel = cmesh->NElements();
    int dim = cmesh->Dimension();
    TPZManVector<REAL,10> globalerror(3,0.0);
    
    if (Serial_ErrorQ) {
        for (int64_t iel = 0; iel < nel; iel++) {
            
            TPZCompEl *cel = cmesh->ElementVec()[iel];
            if (!cel) {
                continue;
            }
            
            if(cel->Reference()->Dimension()!=dim) {
                continue;
            }
            
            TPZManVector<REAL,10> elerror(3,0.);
            elerror.Fill(0.);
            cel->EvaluateError(Analytic, elerror, 0);
            int nerr = elerror.size();
            for (int i=0; i<nerr; i++) {
                globalerror[i] += elerror[i]*elerror[i];
                
            }
            
        }
        
        error_primal    = sqrt(globalerror[0]);
        error_dual      = sqrt(globalerror[1]);
        error_hdiv      = sqrt(globalerror[2]);
    }
    else{
        analysis->SetThreadsForError(nthreads);
        analysis->SetExact(Analytic);
        analysis->PostProcessError(globalerror);
        
        error_primal    = globalerror[0];
        error_dual      = globalerror[1];
        error_hdiv      = globalerror[2];
    }
    
}

/// uncondense the elements unwrap the elements
void UnwrapMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    bool change = true;
    while(change)
    {
        change = false;
        for (int64_t el=0; el<nel; el++) {
            
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
    int64_t nel = mixed_cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
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
                int64_t cindex = mixed_cmesh->AllocateNewConnect(c);
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
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
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
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = mixed_cmesh->Reference();
    gmesh->ResetReference();
    mixed_cmesh->LoadReferences();
    int64_t nelg = gmesh->NElements();
    for (int64_t el=0; el<nelg; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Father() != NULL) {
            continue;
        }
        if (gel->Dimension() == gmesh->Dimension() - 1) {
            continue;
        }
        int64_t mapindex = gel->Index();
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
        int64_t nelst = highlevel.size();
        for (int64_t elst=0; elst<nelst; elst++) {
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
    std::map<int64_t,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " with size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<int64_t>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
#endif
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(mixed_cmesh, ElementGroups, submeshindices, KeepOneLagrangian);
    
    std::cout << "Inserting " << ElementGroups.size()  <<  " macro elements into MHM substructures" << std::endl;
    mixed_cmesh->ComputeNodElCon();
    mixed_cmesh->CleanUpUnconnectedNodes();

    for (std::map<int64_t,int64_t>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = mixed_cmesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        TPZGuiInterface *gui = 0;
        subcmesh->SetAnalysisSkyline(0, 0, gui);
    }
    mixed_cmesh->ComputeNodElCon();
    mixed_cmesh->CleanUpUnconnectedNodes();
    mixed_cmesh->ExpandSolution();
    
}
