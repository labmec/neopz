
// --------------------- std includes ---------------------
#include <iostream>
#include <chrono>

// --------------------- PZ includes ---------------------
#include <pzgmesh.h>
#include <TPZGenGrid2D.h>
#include <TPZVTKGeoMesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzlog.h>

// --------------------- Global variables ---------------------
constexpr int nx = 220;
constexpr int ny = 60;
constexpr int n_cells = nx * ny;
TPZManVector<REAL, n_cells> perm_vec(n_cells, 1);
enum EMatid {ENone,EDomain,EPressureLeft,EPressureRight,ENoFlux};
//const HDivFamily ghdivfam = HDivFamily::EHDivStandard;
const HDivFamily ghdivfam = HDivFamily::EHDivConstant;
//const HDivFamily ghdivfam = HDivFamily::EHDivKernel;

// --------------------- Namespaces ---------------------
using namespace std;

// --------------------- Functions ---------------------
void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);
TPZGeoMesh *CreateSPE10CoarseGeoMesh();
STATE PermeabilityFunction(const TPZVec<REAL> &x);
TPZCompMesh* CreateCompMeshFlux(TPZGeoMesh* gmesh);
TPZCompMesh* CreateCompMeshPressure(TPZGeoMesh* gmesh);
TPZMultiphysicsCompMesh* CreateCompMeshMultiphysics(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> meshvector);
void PrintResultsVTK(const int dim, TPZLinearAnalysis &an, const std::string &plotfile);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main(){
    cout << "\n--------------------- Starting SPE10 simulations ---------------------\n" << endl;
    auto start_time = std::chrono::steady_clock::now();
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    cout << "\n--------------------- Reading Permeability Data ---------------------\n" << endl;
    constexpr int layer = 36; // Same as in repo ErrorEstimation/Projects/SPE10. I suppose is a layer with interesting permeability
    constexpr int dim{2};
    ReadSPE10CellPermeabilities(&perm_vec, layer);
    REAL min = std::numeric_limits<int>::max(), max = -std::numeric_limits<int>::max();
    for (REAL perm : perm_vec) {
        if (perm > max) {
            max = perm;
        }
        if (perm < min) {
            min = perm;
        }
    }
    cout << "\nMinimum permeability = " << min << endl;
    cout << "Maximum permeability = " << max << endl;
    
    cout << "\n--------------------- Creating GeoMesh ---------------------\n" << endl;
    TPZGeoMesh *gmesh = CreateSPE10CoarseGeoMesh();
    ofstream outvtk("spe10gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outvtk);
    
    cout << "\n--------------------- Creating CompMeshes ---------------------\n" << endl;
    TPZCompMesh* cmeshflux = CreateCompMeshFlux(gmesh);
    TPZCompMesh* cmeshpressure = CreateCompMeshPressure(gmesh);
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshflux;
    meshvector[1] = cmeshpressure;
    TPZMultiphysicsCompMesh* mpmesh = CreateCompMeshMultiphysics(gmesh,meshvector);
    
    cout << "\n--------------------- Creating Analysis ---------------------\n" << endl;
    auto start_time_anal = std::chrono::steady_clock::now();
    TPZLinearAnalysis an(mpmesh, true);
    auto total_time_anal = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_anal).count()/1000.;
    cout << "\nTotal time opt band = " << total_time_anal << " seconds" << endl;
        
    cout << "\n--------------------- Solving system ---------------------\n" << endl;
    SolveProblemDirect(an, mpmesh);
    
    cout << "\n--------------------- Post processing ---------------------\n" << endl;
    auto start_time_pp = std::chrono::steady_clock::now();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, mpmesh);
    string outres = "HDivResults.vtk";
    PrintResultsVTK(dim, an, outres);
    auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count()/1000.;
    cout << "Total time post process = " << total_time_pp << " seconds" << endl;
        
    cout << "\n--------------------- End of execution ---------------------\n" << endl;
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "Total time = " << total_time << " seconds" << endl << endl;
    
    return 0;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {
    // Fuction copied from ErrorEstimation/Projects/SPE10
    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("../InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const auto n_cells = perm_vec->size();
    const auto start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

TPZGeoMesh *CreateSPE10CoarseGeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {220., 60., 0.}; // size of domain
    
//    const TPZManVector<int, 3> ndiv = {13, 3, 0};
    const TPZManVector<int, 3> ndiv = {nx, ny, 0}; // 1 unit cell per permeability
//    const TPZManVector<int, 3> ndiv = {1, 1, 0};

    TPZGenGrid2D gen(ndiv, x0, x1);

    auto gmesh = new TPZGeoMesh;
    gen.Read(gmesh);

    gen.SetBC(gmesh, 4, ENoFlux); // bot
    gen.SetBC(gmesh, 5, EPressureRight); // right
    gen.SetBC(gmesh, 6, ENoFlux); // top
    gen.SetBC(gmesh, 7, EPressureLeft); // left
    
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

STATE PermeabilityFunction(const TPZVec<REAL> &x) {
    auto rounded_x = static_cast<int>(x[0]);
    auto rounded_y = static_cast<int>(x[1]);
    if (rounded_x == 220) rounded_x = 219;
    if (rounded_y == 60) rounded_y = 59;
    return perm_vec[rounded_x * 60 + rounded_y];
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

TPZCompMesh* CreateCompMeshFlux(TPZGeoMesh* gmesh){
    
    // ===> Creating compmesh
    constexpr int dim{2};
    const int pOrder = ghdivfam == HDivFamily::EHDivConstant ? 2 : 1;
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // Create domain materials
    {
        const int nstate = 1;
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain,dim,nstate);
        cmesh->InsertMaterialObject(mat);
    }
    // Create bc materials
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(ENoFlux);
        mat->SetDimension(dim);
        cmesh->InsertMaterialObject(mat);
    }
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(EPressureLeft);
        mat->SetDimension(dim);
        cmesh->InsertMaterialObject(mat);
    }
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(EPressureRight);
        mat->SetDimension(dim);
        cmesh->InsertMaterialObject(mat);
    }

    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().HDivFam() = ghdivfam;
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    
//    ofstream out("cmesh.txt");
//    cmesh->Print(out);
    
    return cmesh;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

TPZCompMesh* CreateCompMeshPressure(TPZGeoMesh* gmesh){
    
    // ===> Creating compmesh
    constexpr int dim{2};
    const int pOrder = ghdivfam == HDivFamily::EHDivStandard ? 1 : 0;
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // ===> Domain material
    const int nstate = 1;
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain,dim,nstate);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->SetDefaultOrder(pOrder);
    
    // ===> P order distinction
    if(pOrder == 0)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    cmesh->AutoBuild();
    
    // ===> Set lagrange multiplier for order of assembly
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
//    ofstream out("pmesh.txt");
//    cmesh->Print(out);
    
    return cmesh;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

TPZMultiphysicsCompMesh* CreateCompMeshMultiphysics(TPZGeoMesh* gmesh, TPZVec<TPZCompMesh *> meshvector){
    
    // ===> Creating mp cmesh
    constexpr int dim{2};
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    
    // ===> Domain mat
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
//    mat->SetConstantPermeability(1.);
    std::function<STATE(const TPZVec<REAL> &coord)> func = PermeabilityFunction;
    mat->SetPermeabilityFunction(func);
    cmesh->InsertMaterialObject(mat);
    
    // ===> Boundary Conditions
    constexpr int neumann{1};
    constexpr int dirichlet{0};
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2_0(1,0.), val2_1(1,1.);
    // domain bcs
    auto * BCondNoFlux = mat->CreateBC(mat, ENoFlux, neumann, val1, val2_0);
    cmesh->InsertMaterialObject(BCondNoFlux);
    
    auto * BCondLeft = mat->CreateBC(mat,EPressureLeft, dirichlet, val1, val2_0);
    cmesh->InsertMaterialObject(BCondLeft);
    
    auto * BCondRight = mat->CreateBC(mat,EPressureRight, dirichlet, val1, val2_1);
    cmesh->InsertMaterialObject(BCondRight);
        
    // ===> Build space
    TPZManVector<int> active(2,1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    
    return cmesh;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void PrintResultsVTK(const int dim, TPZLinearAnalysis &an, const std::string &plotfile){
    TPZManVector<std::string,2> scalnames(2), vecnames(1);
    
    scalnames[0] = "Permeability";
    scalnames[1] = "Pressure";
    vecnames[0]= "Flux";
    
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    constexpr int nThreads{8};
    //    TPZFStructMatrix<STATE> matskl(cmesh); // slowest - good for debugging
    //    TPZSkylineStructMatrix<STATE> matskl(cmesh); // medium speed - pz only
    TPZSSpStructMatrix<STATE> matskl(cmesh); // fast - works great with mkl

    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    auto start_time_ass = std::chrono::steady_clock::now();
    cout << "Doing assemble..." << endl;
    an.Assemble();
    auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count()/1000.;
    cout << "Total time assemble = " << total_time_ass << " seconds" << endl;

//  {
//    std::ofstream outmat("mat.nb");
//    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
//    matsol->Matrix()->Print("singmat=",outmat,EMathematicaInput);
//    std::ofstream outrhs("rhs.nb");
//      TPZFMatrix<STATE> rhs = an.Rhs();
//    rhs.Print("rhs=",outrhs,EMathematicaInput);
//  }
    
    ///solves the system
    auto start_time_solve = std::chrono::steady_clock::now();
    cout << "\nDoing solve..." << endl;
    an.Solve();
    auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
    cout << "Total time solve = " << total_time_solve << " seconds" << endl;
    
//    {
//        std::ofstream outsol("sol.nb");
//        TPZFMatrix<STATE> sol = an.Solution();
//        sol.Print("sol = ",outsol,EMathematicaInput);
//    }
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
