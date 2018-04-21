
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
#include "pzcompel.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

// Materials
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TPZPoroPermCoupling.h"
#include "TPZPoroPermCoupling3D.h"

// Elasticity
#include "TPZElasticCriterion.h"

// Plasticity
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"

// Analysis
#include "pzanalysis.h"
#include "TPZPoroPermAnalysis.h"

// Matrix
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"

// Simulation data structure
#include "TPZSimulationData.h"

// Methods declarations
//#define USING_Pardiso

// Apply the mesh refinement
void UniformRefinement(TPZGeoMesh *gmesh, int nh);
void UniformRefinement(TPZGeoMesh * gmesh, int nh, int mat_id);

// Create some functions
static void Sigma(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);
static void u_y(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);
static void u_xy(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP);

// Create a computational mesh for Deformation
TPZCompMesh * CMesh_Deformation(TPZSimulationData * sim_data);

// Create a computational mesh for Pore Pressure
TPZCompMesh * CMesh_PorePressure(TPZSimulationData * sim_data);

// Create a computational mesh for PorePerm Coupling
TPZCompMesh * CMesh_PorePermCoupling(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPZSimulationData * sim_data);

#ifdef LOG4CXX
static LoggerPtr log_data(Logger::getLogger("pz.PorePermCoupling"));
#endif


// Shear-enhanced compaction and strain localization:
// Inelastic deformation and constitutive modeling of four porous sandstones

/**
 * This function generate the data associated to the Figure 2a. Darley Dale Sandstone Data for Cap Model
 *
 */
void LEDSPorosityReductionPlot();



// Restructuring implementation

int main(int argc, char *argv[])
{
    
    LEDSPorosityReductionPlot();
    
    return 0;
    
    
#ifdef LOG4CXX
    if(log_data->isInfoEnabled())
    {
        std::stringstream sout;
        sout << " Defining the case MS -> Review... " << std::endl;
        LOGPZ_DEBUG(log_data,sout.str())
    }
#endif
    
    //    Reading arguments
    char *simulation_file = NULL;
    {
        using namespace std;
        if (argc != 2)
        {
            cout << "Size: " << argc << " Number of Arguments " << endl;
            cout << "Usage: " << argv[0] << " Myinputfile.xml " << endl;
            cout <<    "Program stop: not xml file found \n" << endl;
            DebugStop();
        }
        
        if (argc == 2)
        {
            cout << "Control File used : " << argv[1] << "\n" << endl;
            simulation_file        = argv[1];
        }
    }
    
    // Simulation data to be configurated
    // First a linear poroelasticity kernel.
    TPZSimulationData * sim_data = new TPZSimulationData;
    sim_data->ReadSimulationFile(simulation_file);
    
#ifdef PZDEBUG
    sim_data->PrintGeometry();
#endif
    
    // Create multiphysisc mesh
    TPZManVector<TPZCompMesh * , 2 > mesh_vector(2);
    TPZCompMesh * cmesh_poro_perm_coupling = CMesh_PorePermCoupling(mesh_vector,sim_data);
    
    // The initial condition is set up to zero for Deformation and Pore Pressure
    // Create and run the Transient analysis
    
    bool mustOptimizeBandwidth = true;
    int number_threads = 4;
    TPZPoroPermAnalysis * time_analysis = new TPZPoroPermAnalysis;
    time_analysis->SetCompMesh(cmesh_poro_perm_coupling,mustOptimizeBandwidth);
    time_analysis->SetSimulationData(sim_data);
    time_analysis->SetMeshvec(mesh_vector);
    time_analysis->AdjustVectors();
    
#ifdef USING_Pardiso
    TPZSymetricSpStructMatrix struct_mat(cmesh_poro_perm_coupling); // Pardiso MKL flag
#else
    
    //TPZSkylineNSymStructMatrix struct_mat(cmesh_poro_perm_coupling);
    //    TPZSkylineStructMatrix struct_mat(cmesh_poro_perm_coupling);
    
    TPZParFrontStructMatrix<TPZFrontSym<STATE> > struct_mat(cmesh_poro_perm_coupling);
    struct_mat.SetDecomposeType(ELDLt);
    
#endif

    
    TPZStepSolver<STATE> step;
    struct_mat.SetNumThreads(number_threads);
    step.SetDirect(ELDLt);
    time_analysis->SetSolver(step);
    time_analysis->SetStructuralMatrix(struct_mat);
    
    TPZVec<REAL> x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    std::string file_ss_name("plot.nb");
    std::string file_sp_name("porosity.nb");
    std::string file_sk_name("permeability.nb");
    std::string file_spex_name("porepressure.nb");
    
    // Run Transient analysis
    time_analysis->Run_Evolution(x);
//    time_analysis->PlotStrainStress(file_ss_name);
//    time_analysis->PlotStrainPorosity(file_sp_name);
//    time_analysis->PlotStrainPermeability(file_sk_name);
//    time_analysis->PlotStrainPressure(file_spex_name);
    std::cout << " Execution finished" << std::endl;

    
    
	return EXIT_SUCCESS;
}


// Apply the mesh refinement
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


// Create some functions
void Sigma(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP)
{
    
    REAL MPa = 1.0e6;
    REAL scale = 2.0e4;
//    REAL s_n = 18.0*(scale*time)*MPa;
    REAL s_n = 25.0*(scale*time)*MPa;
    
    f[0] = 0.0;
    f[1] = -s_n;
    f[2] = 0.0;
    return;
}

void u_y(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& GradP)
{
    REAL scale = 1.0e5;//2.0e4;
    REAL uy = (2.0*(0.00028)*time*scale) + 0.0002;
    
    
    f[0] = 0.0;
    f[1] = -uy;
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


// Create a computational mesh for Deformation
TPZCompMesh * CMesh_Deformation(TPZSimulationData * sim_data){
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = dim;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, TPZManVector<int,20>>,20>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int dirichlet = 0;
        TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
        int n_bc = material_ids[iregion].second.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second [ibc];
            TPZMaterial * bc = material->CreateBC(material, bc_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->ElasticityOrder());
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshDeformation.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}


// Create a computational mesh for Pore Pressure
TPZCompMesh * CMesh_PorePressure(TPZSimulationData * sim_data){
    
    
    // Getting mesh dimension
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    int nstate = 1;
    TPZVec<STATE> sol;
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, TPZManVector<int,20>>,20>  material_ids = sim_data->MaterialIds();
    for (int iregion = 0; iregion < n_regions; iregion++)
    {
        int matid = material_ids[iregion].first;
        TPZL2Projection * material = new TPZL2Projection(matid,dim,nstate,sol);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int dirichlet = 0;
        TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);
        int n_bc = material_ids[iregion].second.size();
        for (int ibc = 0; ibc < n_bc; ibc++)
        {
            int bc_id = material_ids[iregion].second [ibc];
            TPZMaterial * bc = material->CreateBC(material, bc_id, dirichlet, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    // Setting H1 approximation space
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(sim_data->DiffusionOrder());
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef PZDEBUG
    std::ofstream out("CmeshPorePressure.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}


// Create a computational mesh for PorePerm Coupling
TPZCompMesh * CMesh_PorePermCoupling(TPZManVector<TPZCompMesh * , 2 > & mesh_vector, TPZSimulationData * sim_data){
    
    mesh_vector[0] = CMesh_Deformation(sim_data);
    mesh_vector[1] = CMesh_PorePressure(sim_data);
    
    // Plane strain assumption
    int planestress = 0;
    int dim = sim_data->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(sim_data->Geometry());
    
    std::map<int, std::string>::iterator it_bc_id_to_type;
    std::map< std::string,std::pair<int,std::vector<std::string> > >::iterator  it_condition_type_to_index_value_names;
    std::map<int , std::vector<REAL> >::iterator it_bc_id_to_values;
    
    REAL to_MPa = 1.0e6;
    REAL to_rad = M_PI/180.0;
    
    int n_regions = sim_data->NumberOfRegions();
    TPZManVector<std::pair<int, TPZManVector<int,20>>,20>  material_ids = sim_data->MaterialIds();
    TPZManVector<TPZManVector<REAL,20>,20> mat_props = sim_data->MaterialProps();
    for (int iregion = 0; iregion < n_regions; iregion++) {
        int matid = material_ids[iregion].first;
        TPZPoroPermCoupling * material = new TPZPoroPermCoupling(matid,dim);
        
        int eyoung = 0, nu = 1, phi = 2, kappa = 3, alpha = 4, m = 5, rho = 6, mu = 7;
        int n_parameters = mat_props[iregion].size();
        
#ifdef PZDEBUG
        if (n_parameters != 8) { // 8 for linear poroelasticity
            DebugStop();
        }
#endif
        
        int kmodel = 0;
        REAL Ey_r = mat_props[iregion][eyoung];
        REAL nu_r = mat_props[iregion][nu];
        REAL alpha_r = mat_props[iregion][alpha];
        REAL Se = 1.0/mat_props[iregion][m];
        REAL k = mat_props[iregion][kappa];
        REAL porosity = mat_props[iregion][phi];
        REAL eta = mat_props[iregion][mu];
        
        REAL c = 1010.0*to_MPa;
        REAL phi_f = 45.0*to_rad;
        
        material->SetSimulationData(sim_data);
        material->SetPlaneProblem(planestress);
        material->SetPorolasticParametersEngineer(Ey_r, nu_r);
        material->SetBiotParameters(alpha_r, Se);
        material->SetParameters(k, porosity, eta);
        material->SetKModel(kmodel);
        material->SetDruckerPragerParameters(phi_f, c);
        cmesh->InsertMaterialObject(material);
        
        // Inserting boundary conditions
        int n_bc = material_ids[iregion].second.size();
        for (int ibc = 0; ibc < n_bc; ibc++) {
            int bc_id = material_ids[iregion].second [ibc];
            
            it_bc_id_to_type = sim_data->BCIdToConditionType().find(bc_id);
            it_bc_id_to_values = sim_data->BCIdToBCValues().find(bc_id);
            it_condition_type_to_index_value_names = sim_data->ConditionTypeToBCIndex().find(it_bc_id_to_type->second);
            
            int bc_index = it_condition_type_to_index_value_names->second.first;
            int n_bc_values = it_bc_id_to_values->second.size();
            TPZFMatrix<STATE> val1(0,0,0.), val2(3,1,0.);
            for (int i = 0; i < n_bc_values; i++) {
                REAL value = it_bc_id_to_values->second[i];
                val2(i,0) = value;
            }
            TPZMaterial * bc = material->CreateBC(material, bc_id, bc_index, val1, val2);
            cmesh->InsertMaterialObject(bc);
        }
    }
    
    // Setting up multiphysics functions
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    // Transfer to multiphysic mesh
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
    std::ofstream out("PorePermCoupling.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}


// The function to generate the Cap Model
void LEDSPorosityReductionPlot(){
    
    
    // DS Dimaggio Sandler PV
    TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse> LEDS;
    
    // LE Linear elastic response
    TPZElasticResponse ER;
    
    
    /**
     * Input data for MCormick Ranch sand data:
     *
     */
    
//    // MCormick Ranch sand data:
//    REAL K = 66.67; // ksi
//    REAL G = 40.00; // ksi
//    
//    REAL E       = (9.0*K*G)/(3.0*K+G);
//    REAL nu      = (3.0*K - 2.0*G)/(2.0*(3.0*K+G));
//    REAL CA      = 0.250;
//    REAL CB      = 0.670;
//    REAL CC      = 0.180;
//    REAL CD      = 0.670;
//    REAL CR      = 1.3;
//    REAL CW      = 0.066;
//    REAL phi = 0, psi = 1., N = 0.;
    
//        REAL Pc = -0.03; // MPa
//    
//    //    REAL alpha = 1.0;
//    //    REAL Pp = 10.0; // MPa
//    
//    REAL epsilon_rate = 1.0e-4;
//    REAL epsilon_lateral = 0.01;
    
//        int64_t n_steps = 100;
    
    
    
//    /**
//     * Input data of Darley Dale Sanstone for Cap Model (First data):
//     *
//     */
//    
//    // Darley Dale Sandstone Data for Cap Model:
//    REAL E       = 19.45e3; // MPa
//    REAL nu      = 0.31;
//    
//    REAL K = E/3*(1-2*nu); // MPa
//    REAL G = E/2*(1+nu); // MPa
//    
//    REAL CA      = 182.44;
//    REAL CB      = -0.00236;
//    REAL CC      = 0.002441;
//    REAL CD      = 0.452;
//    REAL CR      = 3.48927;
//    REAL CW      = 5.5651e-55;
//    REAL phi = 0, psi = 1., N = 0.;
//    
//    REAL Pc = -0.03; // MPa
//    
//    //    REAL alpha = 1.0;
//    //    REAL Pp = 10.0; // MPa
//    
//    REAL epsilon_rate = 3.0e-4;
//    REAL epsilon_lateral = 0.01;
//    
//        int64_t n_steps = 30;
    
    
    
    /**
     * Input data of Darley Dale Sanstone for Cap Model (Second data):
     *
     */
    
    // Data for Cap Model: https://onlinelibrary.wiley.com/doi/pdf/10.1002/nag.1610110208
    REAL E       = 22466.8; // MPa
    REAL nu      = 0.253167;
    
    REAL K = E/(3*(1-2*nu)); // MPa
    REAL G = E/(2*(1+nu)); // MPa
    
    REAL CA      = 690;
    REAL CB      = 3.92e-4;
    REAL CC      = 674.0;
    REAL CD      = 1.45e-3;
    REAL CR      = 4.0;
    REAL CW      = 0.08;
    REAL phi = 0, psi = 1., N = 0.;
    
    REAL Pc = -467.0/3.0; // MPa
    
    //    REAL alpha = 1.0;
    //    REAL Pp = 10.0; // MPa
    
    REAL epsilon_rate = 1.0e-4;
    REAL epsilon_lateral = 0.01;
    
    int64_t n_steps = 100;
    
    
    ER.SetUp(E, nu);
    LEDS.SetElasticResponse(ER);
    LEDS.fYC.SetUp(CA, CB, CC, CD, K, G, CW, CR, phi, N, psi);
    
    TPZTensor<REAL> epsilon_t,sigma;
    sigma.Zero();
    epsilon_t.Zero();
    
    
    sigma.XX() = Pc;
    sigma.YY() = Pc;
    sigma.ZZ() = Pc;
    
    // Initial damage data
    REAL k_0;
    LEDS.InitialDamage(sigma, k_0);
    LEDS.fN.fAlpha = k_0;
    

    

    TPZFNMatrix<80,STATE> LEDS_epsilon_stress(n_steps,2);
    
    TPZPlasticState<STATE> plastic_state;
    for (int64_t t = 0; t < n_steps; t++) {
        
        epsilon_t.XX() = - t * epsilon_rate;
        epsilon_t.YY() = - 1.1 * t * epsilon_rate;
        epsilon_t.ZZ() = - 1.1 * t * epsilon_rate;
        
        LEDS.ApplyStrainComputeSigma(epsilon_t, sigma);
        
        LEDS_epsilon_stress(t,0) = epsilon_t.I1();
        LEDS_epsilon_stress(t,1) = sigma.I1();
//        LEDS_epsilon_stress(t,2) = sqrt(sigma.J2());
        
    }

    LEDS_epsilon_stress.Print("data = ", std::cout,EMathematicaInput);
    
}
