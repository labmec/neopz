
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

#include "TPZMatElastoPlastic2D.h"
#include "TPZMatElastoPlastic.h"
#include "TPZElastoPlasticMem.h"
#include "TPZElasticCriterion.h"
#include "TPZPlasticStepPV.h"
#include "TPZSandlerExtended.h"
#include "TPZYCMohrCoulombPV.h"
#include "pzstepsolver.h"

#include "Headers/TPZSolveMatrix.h"

#ifdef USING_TBB
#include "tbb/parallel_for_each.h"
#endif

TPZGeoMesh *Geometry2D(int nelem_x, int nelem_y, REAL len, int ndivide);

TPZCompMesh *CmeshElasticity(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CmeshElasticityNoBoundary(TPZGeoMesh *gmesh, int pOrder);

TPZCompMesh *CmeshElastoplasticity(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *CmeshElastoplasticityNoBoundary(TPZGeoMesh *gmesh, int pOrder);

void SolMatrix(TPZFMatrix<REAL> residual, TPZCompMesh *cmesh);
void SolVector(TPZFMatrix<REAL> residual, TPZCompMesh *cmesh);

TPZFMatrix<REAL>  Residual(TPZCompMesh *cmesh, TPZCompMesh *cmesh_noboundary);

int main(int argc, char *argv[]) {
//// ------------------------ DATA INPUT ------------------------------
    int nelem_x = atoi(argv[1]); // Number of elements in x direction
    int nelem_y = atoi(argv[1]); // Number of elements in y direction
    REAL len = 1; // Domain length
    int pOrder = 1; // Computational mesh order
    int ndivide = 0; // Subdivision of elements

// Generates the geometry
    TPZGeoMesh *gmesh = Geometry2D(nelem_x, nelem_y, len, ndivide);
// Creates the computational mesh
    TPZCompMesh *cmesh = CmeshElasticity(gmesh, pOrder);
    TPZCompMesh *cmesh_noboundary = CmeshElasticityNoBoundary(gmesh, pOrder);

// Defines the analysis
    bool optimizeBandwidth = true;
    int n_threads = 4;
    TPZAnalysis an(cmesh, optimizeBandwidth);
    TPZSymetricSpStructMatrix strskyl(cmesh);
    strskyl.SetNumThreads(n_threads);
    an.SetStructuralMatrix(strskyl);

// Solve
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
//    an.Solver().Matrix()->Print("K = ", std::cout, EMathematicaInput);
    an.Solve();
//    an.Solution().Print("U = ", std::cout, EMathematicaInput);


// Post process
//    TPZManVector<std::string> scalarnames(0), vecnames(1);
//    scalarnames[0] = "SigmaX";
//    scalarnames[1] = "SigmaY";
//    scalarnames[0] = "StressI1";
//    vecnames[0] = "Displacement";
//    std::string namefile = "Elasticity_teste";
//    an.DefineGraphMesh(2, scalarnames, vecnames, namefile + "ElasticitySolutions.vtk");
//    an.PostProcess(0,2);

// Calculates residual without boundary conditions
    TPZFMatrix<REAL> residual = Residual(cmesh, cmesh_noboundary);

// Calculates residual using matrix operations and check if the result is ok
    SolMatrix(residual, cmesh);
//    SolVector(residual, cmesh);
    return 0;
}

TPZGeoMesh *Geometry2D(int nelem_x, int nelem_y, REAL len, int ndivide) {
// Creates the geometric mesh
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    int dim = 2;
    gmesh->SetDimension(dim);

// Geometry definitions
    int nnodes_x = nelem_x + 1; //Number of nodes in x direction
    int nnodes_y = nelem_y + 1; //Number of nodes in x direction
    int64_t nelem = nelem_x * nelem_y; //Total number of elements

// Nodes initialization
// Enumeration: vertical order - from the below to the top, and from the left to the right
    TPZManVector<REAL> coord(3, 0.);
    int64_t id, index;
    for (int i = 0; i < nnodes_x; i++) {
        for (int j = 0; j < nnodes_y; j++) {
            id = i * nnodes_y + j;
            coord[0] = (i) * len / (nnodes_x - 1);
            coord[1] = (j) * len / (nnodes_y - 1);
            index = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[index] = TPZGeoNode(id, coord, *gmesh);
        }
    }

// Element connectivities
// Enumeration: vertical order - from the below to the top, and from the left to the right
    TPZManVector<int64_t> connect(4, 0);
    for (int i = 0; i < (nnodes_x - 1); i++) {
        for (int j = 0; j < (nnodes_y - 1); j++) {
            index = (i) * (nnodes_y - 1) + (j);
            connect[0] = (i) * nnodes_y + (j);
            connect[1] = connect[0] + (nnodes_y);
            connect[2] = connect[1] + 1;
            connect[3] = connect[0] + 1;
            gmesh->CreateGeoElement(EQuadrilateral, connect, 1, id); //Allocates and define the geometric element
        }
    }

// Generates neighborhood information
    gmesh->BuildConnectivity();

// Creates the boundary conditions
// Dirichlet
    for (int64_t i = 0; i < nelem_y; i++) {
        TPZGeoEl *gelem = gmesh->Element(i);
        TPZGeoElBC el_boundary(gelem, 7, -1); //Left side of the plane
    }
    for (int64_t i = 0; i < nelem_x; i++) {
        int64_t n = nelem_y * (i + 1) - 1;
        TPZGeoEl *gelem = gmesh->Element(n);
        TPZGeoElBC el_boundary(gelem, 6, -2); //Top side of the plane
    }

// Neumann
    for (int64_t i = 0; i < nelem_x; i++) {
        int64_t n = nelem_y * (i + 1) - (nelem_y);
        TPZGeoEl *gelem = gmesh->Element(n);
        TPZGeoElBC el_boundary(gelem, 4, -3); //Bottom side of the plane - tension
    }
    for (int64_t i = nelem - nelem_y; i < nelem; i++) {
        TPZGeoEl *gelem = gmesh->Element(i);
        TPZGeoElBC el_boundary(gelem, 5, -4); //Right side of the plane - tension
    }

// HP adaptativity
    if (ndivide != 0) {
        // Finding the elements which will be subdivided
        TPZGeoEl *gel; // Defining the element
        TPZVec<REAL> x(3, 0.); // Defining the coordinate at the end of the node
        x[0] = 0;
        x[1] = len;
        TPZVec<REAL> qsi(3, 0.); // Defining the parametric coordinate
        int64_t InitialElIndex = 0;
        int targetDim = 2;
        gel = gmesh->FindElement(x, qsi, InitialElIndex,
                                 targetDim); // Finding the element which is related to the coordinate
        int64_t elid = gel->Index(); // Atention: this procedure catchs the first element which is related to the coordinate

        TPZVec<TPZGeoEl *> subelindex;

        gel = gmesh->Element(elid);
        gel->Divide(subelindex);
        for (int i = 0; i < ndivide - 1; i++) {
            subelindex[3]->Divide(subelindex);
        }

    }
    return gmesh;
}

TPZCompMesh *CmeshElasticity(TPZGeoMesh *gmesh, int pOrder) {

    // Creating the computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);

// Creating elasticity material
    TPZMatElasticity2D *mat = new TPZMatElasticity2D(1);
    mat->SetElasticParameters(200000000., 0.3);
    mat->SetPlaneStrain();

// Setting the boundary conditions
    TPZMaterial *bcBottom, *bcRight, *bcTop, *bcLeft;
    TPZFMatrix<REAL> val1(2, 1, 0.);
    TPZFMatrix<REAL> val2(2, 1, 0.);

    bcLeft = mat->CreateBC(mat, -1, 7, val1, val2); // X displacement = 0
    bcTop = mat->CreateBC(mat, -2, 8, val1, val2); // Y displacement = 0

    val2(1, 0) = -1000000.;
    bcBottom = mat->CreateBC(mat, -3, 1, val1, val2); // Tension in y

    val2(0, 0) = 1000000.;
    val2(1, 0) = 0.0;
    bcRight = mat->CreateBC(mat, -4, 1, val1, val2); // Tension in x

    cmesh->InsertMaterialObject(mat);

    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcRight);
    cmesh->InsertMaterialObject(bcTop);
    cmesh->InsertMaterialObject(bcLeft);

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;


//// Creates the computational mesh
//    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//    cmesh->SetDefaultOrder(pOrder);
//
//// Creates elastic material
//    TPZMatElasticity2D *material = new TPZMatElasticity2D(1);
//    material->SetElasticParameters(200000., 0.3);
//    material->SetPlaneStrain();
//
//// Set the boundary conditions
//    TPZMaterial *bcBottom, *bcRight, *bcTop, *bcLeft;
//    TPZFMatrix<REAL> val1(2, 2), val2(2, 2);
//
//    val2(0, 0) = 0;
//    val2(1, 0) = 0;
//    bcLeft = material->CreateBC(material, -1, 0, val1, val2); // X displacement = 0
//
//    val2(0,0) = 0;
//    val2(1,0) = 0;
//    bcTop = material->CreateBC(material, -2, 0, val1, val2); // Y displacement = 0
//
//    val2(0, 0) = 0.0;
//    val2(1, 0) = -1000.;
//    bcBottom = material->CreateBC(material, -3, 1, val1, val2); // Tension in y
//
//    val2(0, 0) = 0.0;
//    val2(1, 0) = 0.0;
//    bcRight = material->CreateBC(material, -4, 0, val1, val2); // Tension in x
//
//    cmesh->InsertMaterialObject(material);
//    cmesh->InsertMaterialObject(bcBottom);
//    cmesh->InsertMaterialObject(bcRight);
//    cmesh->InsertMaterialObject(bcTop);
//    cmesh->InsertMaterialObject(bcLeft);
//
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->AutoBuild();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//
//    return cmesh;
}

TPZCompMesh *CmeshElasticityNoBoundary(TPZGeoMesh *gmesh, int pOrder) {

    // Creating the computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);

    // Creating elasticity material
    TPZMatElasticity2D *mat = new TPZMatElasticity2D(1);
    mat->SetElasticParameters(200000000., 0.3);
    mat->SetPlaneStrain();
    cmesh->InsertMaterialObject(mat);

    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    return cmesh;
}

TPZCompMesh *CmeshElastoplasticity(TPZGeoMesh * gmesh, int p_order) {

// Creates the computational mesh
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(p_order);

// Mohr Coulomb data
    REAL mc_cohesion    = 100.0;
    REAL mc_phi         = (50.0*M_PI/180);
    REAL mc_psi         = mc_phi;

// ElastoPlastic Material using Mohr Coulomb
// Elastic predictor
    TPZElasticResponse ER;
    REAL nu = 0.3;
    REAL E = 200000;

    TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse> LEMC;
    ER.SetUp(E, nu);
    LEMC.SetElasticResponse(ER);
    LEMC.fYC.SetUp(mc_phi, mc_psi, mc_cohesion, ER);
    int PlaneStrain = 1;
    int matid = 1;

// Creates elastoplatic material
    TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem > * material = new TPZMatElastoPlastic2D < TPZPlasticStepPV<TPZYCMohrCoulombPV, TPZElasticResponse>, TPZElastoPlasticMem >(matid,PlaneStrain);
    material->SetPlasticityModel(LEMC);

// Set the boundary conditions
    TPZMaterial *bcBottom, *bcRight, *bcTop, *bcLeft;
    TPZFMatrix<REAL> val1(2,2), val2(2,2);

    val2(0,0) = 0;
    val2(1,0) = 0;
    bcLeft = material->CreateBC(material, -1, 0, val1, val2);

    val2(0,0) = 0;
    val2(1,0) = 0;
    bcTop = material->CreateBC(material, -2, 0, val1, val2);

    val2(0,0) = 0;
    val2(1,0) = -1000;
    bcBottom = material->CreateBC(material, -3, 1, val1, val2);

    val2(0,0) = 0;
    val2(1,0) = 0;
    bcRight = material->CreateBC(material, -4, 0, val1, val2);

    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcTop);
    cmesh->InsertMaterialObject(bcLeft);
    cmesh->InsertMaterialObject(bcRight);

    cmesh->SetAllCreateFunctionsContinuousWithMem();
    cmesh->AutoBuild();

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
    ER.SetUp(E, nu);
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

void SolVector(TPZFMatrix<REAL> residual, TPZCompMesh *cmesh) {

    int dim_mesh = (cmesh->Reference())->Dimension(); // Mesh dimension
    int64_t nelem_c = cmesh->NElements(); // Number of computational elements
    std::vector<int64_t> cel_indexes;

// Number of domain geometric elements
    for (int64_t i = 0; i < nelem_c; i++) {
        TPZCompEl *cel = cmesh->Element(i);
        if (!cel) continue;
        TPZGeoEl *gel = cmesh->Element(i)->Reference();
        if (!gel || gel->Dimension() != dim_mesh) continue;
        cel_indexes.push_back(cel->Index());
    }

    if (cel_indexes.size() == 0) {
        DebugStop();
    }

// RowSizes and ColSizes vectors
    int64_t nelem = cel_indexes.size();
    TPZVec<int64_t> rowsizes(nelem);
    TPZVec<int64_t> colsizes(nelem);

    int64_t npts_tot = 0;
    int64_t nf_tot = 0;

    for (auto iel : cel_indexes) {
        //Verification
        TPZCompEl *cel = cmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        rowsizes[iel] = dim * npts;
        colsizes[iel] = nf;

        npts_tot += npts;
        nf_tot += nf;
    }

    TPZSolveVector *SolVec = new TPZSolveVector(dim_mesh * npts_tot, nf_tot, rowsizes, colsizes);

// Dphi matrix, weight and indexes vectors
    TPZFMatrix<REAL> elmatrix;
    TPZVec<REAL> weight(npts_tot);
    TPZManVector<MKL_INT> indexes(dim_mesh * nf_tot);
    int cont = 0;
    for (auto iel : cel_indexes) {
        int64_t cont1 = 0;
        int64_t cont2 = 0;
        //Verification
        TPZCompEl *cel = cmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        elmatrix.Resize(dim * npts, nf);
        for (int64_t inpts = 0; inpts < npts; inpts++) {
            TPZManVector<REAL> qsi(dim, 1);
            REAL w;
            int_rule->Point(inpts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);
            weight[iel + nelem*inpts] = w * std::abs(data.detjac);

            TPZFMatrix<REAL> &dphix = data.dphix;
            for (int inf = 0; inf < nf; inf++) {
                for (int idim = 0; idim < dim; idim++)
                    elmatrix(inpts * dim + idim, inf) = dphix(idim, inf);
            }
        }
        SolVec->SetElementMatrix(iel, elmatrix);

        int64_t ncon = cel->NConnects();
        for (int64_t icon = 0; icon < ncon; icon++) {
            int64_t id = cel->ConnectIndex(icon);
            TPZConnect &df = cmesh->ConnectVec()[id];
            int64_t conid = df.SequenceNumber();
            if (df.NElConnected() == 0 || conid < 0 || cmesh->Block().Size(conid) == 0) continue;
            else {
                int64_t pos = cmesh->Block().Position(conid);
                int64_t nsize = cmesh->Block().Size(conid);
                for (int64_t isize = 0; isize < nsize; isize++) {
                    if (isize % 2 == 0) {
                        indexes[cont1*nelem + cont] = pos + isize;
                        cont1++;
                    } else {
                        indexes[cont2*nelem + nf_tot + cont] = pos + isize;
                        cont2++;
                    }
                }
            }
        }
        cont++;
    }
    SolVec->SetIndexes(indexes);
    SolVec->ColoringElements(cmesh);

    TPZFMatrix<REAL> coef_sol = cmesh->Solution();
    int neq = cmesh->NEquations();

    TPZFMatrix<REAL> nodal_forces_global1(neq, 1, 0.);
    TPZFMatrix<REAL> nodal_forces_global2(neq, 1, 0.);
    TPZFMatrix<REAL> result;
    TPZFMatrix<REAL> sigma;
    TPZFMatrix<REAL> nodal_forces_vec;

#ifdef __CUDACC__
    std::cout << "\n\nSOLVING WITH GPU" << std::endl;
    SolVec->AllocateMemory(cmesh);
    SolVec->MultiplyCUDA(coef_sol,result);
    SolVec->ComputeSigmaCUDA(weight, result, sigma);
    SolVec->MultiplyTransposeCUDA(sigma,nodal_forces_vec);
    SolVec->ColoredAssembleCUDA(nodal_forces_vec,nodal_forces_global1);
    SolVec->FreeMemory();

#endif

    std::cout << "\n\nSOLVING WITH CPU" << std::endl;
    SolVec->Multiply(coef_sol, result);
    SolVec->ComputeSigma(weight, result, sigma);
    SolVec->MultiplyTranspose(sigma,nodal_forces_vec);
    SolVec->ColoredAssemble(nodal_forces_vec,nodal_forces_global2);

    //Check the result
    int rescpu = Norm(nodal_forces_global2 - residual);
    if(rescpu == 0){
        std::cout << "\nAssemble done in the CPU is ok." << std::endl;
    } else {
        std::cout << "\nAssemble done in the CPU is not ok." << std::endl;
    }

#ifdef __CUDACC__
    int resgpu = Norm(nodal_forces_global1 - residual);
    if(resgpu == 0){
        std::cout << "\nAssemble done in the GPU is ok." << std::endl;
    } else {
        std::cout << "\nAssemble done in the GPU is not ok." << std::endl;
    }
#endif
}

void SolMatrix(TPZFMatrix<REAL> residual, TPZCompMesh *cmesh) {

    int dim_mesh = (cmesh->Reference())->Dimension(); // Mesh dimension
    int64_t nelem_c = cmesh->NElements(); // Number of computational elements
    std::vector<int64_t> cel_indexes;

// Number of domain geometric elements
    for (int64_t i = 0; i < nelem_c; i++) {
        TPZCompEl *cel = cmesh->Element(i);
        if (!cel) continue;
        TPZGeoEl *gel = cmesh->Element(i)->Reference();
        if (!gel || gel->Dimension() != dim_mesh) continue;
        cel_indexes.push_back(cel->Index());
    }

    if (cel_indexes.size() == 0) {
        DebugStop();
    }

// RowSizes and ColSizes vectors
    int64_t nelem = cel_indexes.size();
    TPZVec<int> rowsizes(nelem);
    TPZVec<int> colsizes(nelem);

    int64_t npts_tot = 0;
    int64_t nf_tot = 0;

    for (auto iel : cel_indexes) {
        //Verification
        TPZCompEl *cel = cmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        rowsizes[iel] = dim * npts;
        colsizes[iel] = nf;

        npts_tot += npts;
        nf_tot += nf;
    }

    TPZSolveMatrix *SolMat = new TPZSolveMatrix(dim_mesh * npts_tot, nf_tot, rowsizes, colsizes);

// Dphi matrix, weight and indexes vectors
    TPZFMatrix<REAL> elmatrix;
    TPZStack<REAL> weight;
    TPZManVector<MKL_INT> indexes(dim_mesh * nf_tot);

    int64_t cont1 = 0;
    int64_t cont2 = 0;

    for (auto iel : cel_indexes) {
        //Verification
        TPZCompEl *cel = cmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        elmatrix.Resize(dim * npts, nf);
        for (int64_t inpts = 0; inpts < npts; inpts++) {
            TPZManVector<REAL> qsi(dim, 1);
            REAL w;
            int_rule->Point(inpts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);
            weight.Push(w * std::abs(data.detjac)); //weight = w * detjac

            TPZFMatrix<REAL> &dphix = data.dphix;
            for (int inf = 0; inf < nf; inf++) {
                for (int idim = 0; idim < dim; idim++)
                    elmatrix(inpts * dim + idim, inf) = dphix(idim, inf);
            }
        }
        SolMat->SetElementMatrix(iel, elmatrix);

        //Indexes vector
        int64_t ncon = cel->NConnects();
        for (int64_t icon = 0; icon < ncon; icon++) {
            int64_t id = cel->ConnectIndex(icon);
            TPZConnect &df = cmesh->ConnectVec()[id];
            int64_t conid = df.SequenceNumber();
            if (df.NElConnected() == 0 || conid < 0 || cmesh->Block().Size(conid) == 0) continue;
            else {
                int64_t pos = cmesh->Block().Position(conid);
                int64_t nsize = cmesh->Block().Size(conid);
                for (int64_t isize = 0; isize < nsize; isize++) {
                    if (isize % 2 == 0) {
                        indexes[cont1] = pos + isize;
                        cont1++;
                    } else {
                        indexes[cont2 + nf_tot] = pos + isize;
                        cont2++;
                    }
                }
            }
        }
    }
    SolMat->SetIndexes(indexes);
    SolMat->ColoringElements(cmesh);

    TPZFMatrix<REAL> coef_sol = cmesh->Solution();
    int neq = cmesh->NEquations();

    TPZFMatrix<REAL> nodal_forces_global1(neq, 1, 0.);
    TPZFMatrix<REAL> nodal_forces_global2(neq, 1, 0.);
    TPZFMatrix<REAL> result;
    TPZFMatrix<REAL> sigma;
    TPZFMatrix<REAL> nodal_forces_vec;


    #ifdef __CUDACC__
    std::cout << "\n\nSOLVING WITH GPU" << std::endl;
    SolMat->AllocateMemory(cmesh);
    SolMat->MultiplyCUDA(coef_sol, result);
    SolMat->MultiplyInThreadsCUDA(coef_sol, result);
//    SolMat->ComputeSigmaCUDA(weight, result, sigma);
//    SolMat->MultiplyTransposeCUDA(sigma, nodal_forces_vec);
//    SolMat->ColoredAssembleCUDA(nodal_forces_vec, nodal_forces_global1);
//    SolMat->FreeMemory();
    #endif

    std::cout << "\n\nSOLVING WITH CPU" << std::endl;
    SolMat->MultiplyInThreads(coef_sol, result);
    SolMat->Multiply(coef_sol, result);
    SolMat->ComputeSigma(weight, result, sigma);
    SolMat->MultiplyTranspose(sigma, nodal_forces_vec);
    SolMat->ColoredAssemble(nodal_forces_vec, nodal_forces_global2);

//    nodal_forces_global1.Print(std::cout);
//    nodal_forces_global2.Print(std::cout);
//    //Check the result
//    int rescpu = Norm(nodal_forces_global2 - residual);
//    if(rescpu == 0){
//        std::cout << "\nAssemble done in the CPU is ok." << std::endl;
//    } else {
//        std::cout << "\nAssemble done in the CPU is not ok." << std::endl;
//    }
//
//    #ifdef __CUDACC__
//    int resgpu = Norm(nodal_forces_global1 - residual);
//    if(resgpu == 0){
//        std::cout << "\nAssemble done in the GPU is ok." << std::endl;
//    } else {
//        std::cout << "\nAssemble done in the GPU is not ok." << std::endl;
//    }
//    #endif
}

TPZFMatrix<REAL> Residual(TPZCompMesh *cmesh, TPZCompMesh *cmesh_noboundary) {
//    bool optimizeBandwidth = true;
//    int n_threads = 4;
//
//    TPZAnalysis an_d(cmesh_noboundary, optimizeBandwidth);
//    TPZSymetricSpStructMatrix strskyl(cmesh_noboundary);
//    strskyl.SetNumThreads(n_threads);
//    an_d.SetStructuralMatrix(strskyl);
//
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt);
//    an_d.SetSolver(step);
//    an_d.Assemble();
//    an_d.Solve();

    TPZFMatrix<STATE> res(cmesh->NEquations(),1,0.);
//    TPZFMatrix<STATE> res;
//    an_d.Solver().Matrix()->Multiply(cmesh->Solution(), res);
    return res;
}
