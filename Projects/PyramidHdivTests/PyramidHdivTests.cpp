/**
 * @file
 * @brief Tests for hdiv pyramid
 * @author Nathan Shauer
 * @since 2016
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>

#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeotetrahedra.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZTimer.h"
#include "pzshapepiramHdiv.h"
#include "mixedpoisson.h"
#include "pzelctemp.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzgengrid.h"
#include "pzpoisson3d.h"
#include "TPZCompMeshTools.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzelchdiv.h"
#include "pzshapetetra.h"
#include "pzelementgroup.h"

#include "TPZVecL2.h"


#include <sys/time.h>

#include "run_stats_table.h"
#ifdef USING_TBB
#include <tbb/tbb.h>
#endif



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pyramtests"));
#endif


using namespace std;

TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref);
TPZGeoMesh *MalhaQuadrada(int &nelx, int &nely);
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
void InsertElasticityCubo(TPZCompMesh *mesh);
void InsertBidimensionalPoisson(TPZCompMesh *cmesh, int &dim);
TPZGeoMesh * CreateGeoMesh1Pir();
TPZGeoMesh * CreateGeoMesh1Tet();
TPZGeoMesh * CreateGeoMeshHexaOfPir();
TPZGeoMesh * CreateGeoMeshHexaOfPirTetra();
TPZGeoMesh * CreateGeoMeshPrism();
TPZCompMesh * CreateCmeshPressure(TPZGeoMesh *gmesh, int &p);
TPZCompMesh * CreateCmeshFlux(TPZGeoMesh *gmesh, int &p);
TPZCompMesh * CreateCmeshMulti(TPZVec<TPZCompMesh *> &meshvec);
void LoadSolution(TPZCompMesh *cpressure);
void ProjectFlux(TPZCompMesh *cfluxmesh);
void GroupElements(TPZCompMesh *cmesh);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &f) {

    f[0] = pt[0];
    return;
}

void FluxFunc(const TPZVec<REAL> &pt, TPZVec<STATE> &flux)
{
    flux[0] = -1.;
}

int main1(int argc, char *argv[]);
int main2(int argc, char *argv[]);



int main(int argc, char *argv[])
{
    bool assembleTests = false;
    if (assembleTests) {
        main1(argc, argv);
    }
    else{
        main2(argc,argv);
    }
    
    return 0;
}

using namespace pzshape;

int main2(int argc, char *argv[])
{
    string projectpath = "/Projects/PyramidHdivTests/";

#ifdef LOG4CXX
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    FileName = dirname + projectpath;
    FileName += "pyramlogfile.cfg";
    InitializePZLOG(FileName);
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream str;
        str << "\nRodando testes de piramede Hdiv" << std::endl;
        LOGPZ_DEBUG(logger,str.str())
    }
#endif

    gRefDBase.InitializeAllUniformRefPatterns();
    HDivPiola = 1;
//    TPZGeoMesh *gmesh = CreateGeoMesh1Pir();
//    TPZGeoMesh *gmesh = CreateGeoMeshHexaOfPir();
    TPZGeoMesh *gmesh = CreateGeoMeshHexaOfPirTetra();
    UniformRefine(gmesh, 1);
    {
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, "../PyramidGMesh.vtk", true);
    }
//    TPZGeoMesh *gmesh = CreateGeoMesh1Tet();
//    TPZGeoMesh *gmesh = CreateGeoMeshPrism();
    int pPressure = 1;
    int pFlux = 1;
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZManVector<TPZCompMesh*,2> meshvec(2);
    meshvec[1] = CreateCmeshPressure(gmesh, pPressure);
    LoadSolution(meshvec[1]);
    meshvec[0] = CreateCmeshFlux(gmesh, pFlux);
    TPZCompMeshTools::AddHDivPyramidRestraints(meshvec[0]);
    ProjectFlux(meshvec[0]);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        meshvec[0]->Print(sout);
        meshvec[1]->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZCompMesh *cmeshMult = CreateCmeshMulti(meshvec);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmeshMult);
//    GroupElements(cmeshMult);
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmeshMult->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZAnalysis an(cmeshMult,false);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    TPZSkylineStructMatrix skyl(cmeshMult);
//    TPZFStructMatrix skyl(cmeshMult);
    an.SetStructuralMatrix(skyl);
    an.SetSolver(step);
    
    an.Assemble();
    TPZAutoPointer<TPZMatrix<STATE> > mat = an.Solver().Matrix();
    TPZFMatrix<STATE> solution = cmeshMult->Solution();
    {
        std::ofstream sol("../SolInterpolate.txt");
        sol << "Solution obtained by interpolation\n";
        an.PrintVectorByElement(sol, solution,1.e-8);
    }
    TPZFMatrix<STATE> rhs;
    mat->Multiply(solution, rhs);
    {
        std::ofstream theory("../RhsbyMult.txt");
        theory << "Rhs obtained by matrix vector multiply\n";
        an.PrintVectorByElement(theory, rhs,1.e-8);
    }
    {
        std::ofstream practice("../RhsComputed.txt");
        practice << "Rhs obtained by assembly process\n";
        an.PrintVectorByElement(practice, an.Rhs(),1.e-8);
    }

    an.Solve();
    {
        std::ofstream sol("../SolComputed.txt");
        sol << "Solution obtained by matrix inversion\n";
        an.PrintVectorByElement(sol, cmeshMult->Solution(),1.e-8);
    }

    int dim = 3;
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    std::string plotfile = "../Pyramid_Solution.vtk";
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    
    int postprocessresolution = 0;
    an.PostProcess(postprocessresolution);
    
    return 0;
}

TPZGeoMesh * CreateGeoMesh1Pir()
{
    long nodeids[]={0,1,2,3,4};
//    REAL coords[5][3]={{-1.,-1.,0},{1.,-1.,0.},{1.,1.,0.},{-1.,1.,0.},{0.,0.,1.}};
    REAL coords[5][3]={{-1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.}};
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);

    // Setando os nohs
    int nnodes = 5;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    long index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
                    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 4
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // Criando elemento
    TPZManVector<long,5> topolPyr(5);
    for (int i = 0; i < 5; i++) {
        topolPyr[i] = i;
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    gel->CreateBCGeoEl(13, bc0); // fundo
    gel->CreateBCGeoEl(14, bc0); // frente (-1,-1 ateh 1,-1)
    gel->CreateBCGeoEl(15, bc0); // direita
    gel->CreateBCGeoEl(16, bc0); // atras
    gel->CreateBCGeoEl(17, bc0); // esquerda
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("1PyrGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}

TPZGeoMesh * CreateGeoMeshPrism()
{
//    long nodeids[]={0,1,2,3,4};
    long nodeids[]={3,0,1,2,4};
    long nodeidstet[] = {2,3,4,5};
    long boundaryids[7][4] = {{0,4,3},{3,4,5},{0,4,1},{3,5,2},{4,1,2},{4,2,5},{0,1,2,3}};
    //    REAL coords[5][3]={{-1.,-1.,0},{1.,-1.,0.},{1.,1.,0.},{-1.,1.,0.},{0.,0.,1.}};
    REAL coords[6][3]={{-1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,-1.,1.}};
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 6;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    long index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 4
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 5
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;

    // Criando elemento
    TPZManVector<long,5> topolPyr(5);
    for (int i = 0; i < 5; i++) {
        topolPyr[i] = nodeids[i];
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    
    TPZManVector<long,4> topolTet(4);
    for (int i=0; i<4; i++) {
        topolTet[i] = nodeidstet[i];
    }
    gel = gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    for (int el=0; el<6; el++) {
        TPZManVector<long,3> nodes(3);
        for (int i=0; i<3; i++) {
            nodes[i] = boundaryids[el][i];
        }
        long index;
        gmesh->CreateGeoElement(ETriangle, nodes, bc0, index);
    }
    TPZManVector<long,4> nodes(4);
    for (int i=0; i<4; i++) {
        nodes[i] = boundaryids[6][i];
    }
    gmesh->CreateGeoElement(EQuadrilateral, nodes, bc0, index);
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("../PrismGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}


TPZGeoMesh * CreateGeoMesh1Tet()
{
    long nodeids[]={0,1,2,3};
//    REAL coords[4][3]={{0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
    REAL coords[4][3]={{-1.,-1.,1.},{1.,-1.,1.},{1.,1.,1.},{1.,-1.,-1.}};


    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 5;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    long index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 1
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 2
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    // noh 3
    for (int i=0; i<3; i++) {
        nodecoord[i] = coords[ino][i];
    }
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(nodeids[ino]);
    ino++;
    
    
    // Criando elemento
    TPZManVector<long,5> topolTet(4);
    for (int i = 0; i < 4; i++) {
        topolTet[i] = i;
    }
    
    TPZGeoEl *gel = gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index);
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    gel->CreateBCGeoEl(10, bc0); // fundo
    gel->CreateBCGeoEl(11, bc0); // frente (-1,-1 ateh 1,-1)
    gel->CreateBCGeoEl(12, bc0); // direita
    gel->CreateBCGeoEl(13, bc0); // atras
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("../1TetGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
    
}


TPZGeoMesh * CreateGeoMeshHexaOfPir()
{
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 9;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    long index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 1
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 2
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 3
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 4
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;

    // noh 5
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 6
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 7
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 8
    nodecoord[0] = 0.;
    nodecoord[1] = 0.;
    nodecoord[2] = 0.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    
    // Criando elemento
    TPZManVector<long,5> topolPyr(5);
    int myels[6][5] = {{0,1,5,4,8},{1,2,6,5,8},{2,3,7,6,8},{0,3,7,4,8},{0,1,2,3,8},{4,5,6,7,8}};
    //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
    for (int iel = 0; iel < 6; iel++) {
        for (int i = 0; i < 5; i++) {
            topolPyr[i] = myels[iel][i];
        }
        gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    }
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    
    const long nel = gmesh->NElements();
    for (long iel = 0; iel < nel; iel++) {
        gmesh->Element(iel)->CreateBCGeoEl(13, bc0);
    }

    gmesh->BuildConnectivity();
    
    std::ofstream out("HexaPyrGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
}

TPZGeoMesh * CreateGeoMeshHexaOfPirTetra()
{
    const int dim = 3;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(dim);
    
    // Setando os nohs
    int nnodes = 9;
    gmesh->NodeVec().Resize(nnodes);
    int ino = 0;
    const int matid = 1;
    long index = 0;
    
    // noh 0
    TPZManVector<REAL, 3> nodecoord(3,0.);
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 1
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 2
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 3
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = -1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 4
    nodecoord[0] = -1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 5
    nodecoord[0] = 1.;
    nodecoord[1] = -1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 6
    nodecoord[0] = 1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 7
    nodecoord[0] = -1.;
    nodecoord[1] = 1.;
    nodecoord[2] = 1.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
    ino++;
    
    // noh 8
    nodecoord[0] = 0.;
    nodecoord[1] = 0.;
    nodecoord[2] = 0.;
    gmesh->NodeVec()[ino].SetCoord(nodecoord);
    gmesh->NodeVec()[ino].SetNodeId(ino);
//    ino++;
    gmesh->SetNodeIdUsed(ino);
    
    // Criando elemento
    TPZManVector<long,5> topolPyr(5), topolTet(4), topolTri(3);
    int myelsp[2][5] = {{4,0,2,6,1},{4,0,2,6,7}};
    int myelst[2][4] = {{4,6,5,1},{0,2,3,7}};
    //                          front           right          top             back            left            bottom
    int triangles[12][3] = {{0,1,4},{1,5,4},{1,2,6},{1,6,5},{4,5,6},{4,6,7},{2,6,7},{2,7,3},{0,3,7},{0,7,4},{0,1,2},{0,2,3} };
    //int myels[6][5] = {{0,1,5,4,8},{6,5,1,2,8},{2,3,7,6,8},{7,4,0,3,8},{0,1,2,3,8},{4,5,6,7,8}}; //Sequencia trocada soh para funcionar o AddHDivPyramidRestraints
    for (int iel = 0; iel < 2; iel++) {
        for (int i = 0; i < 5; i++) {
            topolPyr[i] = myelsp[iel][i];
        }
        gmesh->CreateGeoElement(EPiramide, topolPyr, matid, index);
    }
    for (int iel = 0; iel < 2; iel++) {
        for (int i = 0; i < 4; i++) {
            topolTet[i] = myelst[iel][i];
        }
        gmesh->CreateGeoElement(ETetraedro, topolTet, matid, index);
    }
    
    const int bc0 = -1;//, bc1 = -2, bc2 = -3, bc3 = -4, bc4 = -5;
    
    for (long iel = 0; iel < 12; iel++) {
        for (int i = 0; i < 3; i++) {
            topolTri[i] = triangles[iel][i];
        }
        gmesh->CreateGeoElement(ETriangle, topolTri, bc0, index);
    }
    
    gmesh->BuildConnectivity();
    
    std::ofstream out("../HexaPyrTetGmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    return gmesh;
}


TPZCompMesh * CreateCmeshPressure(TPZGeoMesh *gmesh, int &p)
{
    const int matid = 1;
    const int dim = 3;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(p);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    TPZMixedPoisson *mymat = new TPZMixedPoisson(matid, dim);
    cmesh->InsertMaterialObject(mymat);

    const long nel = gmesh->NElements();
    long index;
    for (long iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        if (!gel || gel->Type() != EPiramide){
            long index;
            if (gel->Dimension() == gmesh->Dimension()) {
                cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, index);
            }
        }
        else
        {
            new TPZIntelGen<TPZShapePiramHdiv>(*cmesh,gel,index);
        }
        gel->ResetReference();
    }
    cmesh->ExpandSolution();
    return cmesh;
}

TPZCompMesh * CreateCmeshFlux(TPZGeoMesh *gmesh, int &p)
{
    const int matid = 1, bc0 = -1;
    const int dim = 3;
    const int dirichlet = 0;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(p);
    
    TPZVecL2 *mymat = new TPZVecL2(matid);
    mymat->SetForcingFunction(FluxFunc);
    cmesh->InsertMaterialObject(mymat);
    
    TPZFMatrix<> val1(3,3,0.);
    TPZFMatrix<> val2(3,1,0.);
    TPZBndCond *bnd = mymat->CreateBC(mymat, bc0, dirichlet, val1, val2);
    cmesh->InsertMaterialObject(bnd);
    
    cmesh->SetAllCreateFunctionsHDiv();
    
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh * CreateCmeshMulti(TPZVec<TPZCompMesh *> &meshvec)
{
    const int matid = 1, bc0 = -1;
    const int dirichlet = 0;
    //Creating computational mesh for multiphysic elements
    TPZGeoMesh *gmesh = meshvec[0]->Reference();
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    //criando material
    int dim = gmesh->Dimension();
    mphysics->SetDimModel(dim);
    
    TPZMixedPoisson *mat = new TPZMixedPoisson(matid,dim);
    
    //inserindo o material na malha computacional
    mphysics->InsertMaterialObject(mat);
    
    
    //Criando condicoes de contorno
    TPZFMatrix<> val1(3,3,0.), val2(3,1,1.);
    TPZBndCond * BCond0 = mat->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(Forcing);
    BCond0->SetForcingFunction(0,force);
    
    mphysics->InsertMaterialObject(BCond0);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    
    //Fazendo auto build
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    //Creating multiphysic elements containing skeletal elements.
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
    
    meshvec[0]->CleanUpUnconnectedNodes();
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->ExpandSolution();


#ifdef LOG4CXX
    std::stringstream sout;
    mphysics->Print(sout);
    LOGPZ_DEBUG(logger,sout.str())
#endif
    
    return mphysics;
}

int main1(int argc, char *argv[])
{
    
//    TTimer tref;
//    tref.start();
    gRefDBase.InitializeUniformRefPattern(ETetraedro);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EPiramide);
//    tref.stop();
//    cout << "\ntempo de refpattern = " << tref.seconds() << endl;
    
    
    string projectpath = "/Projects/PyramidHdivTests/";

#ifdef LOG4CXX
    if (logger->isDebugEnabled()){
        std::string dirname = PZSOURCEDIR;
        std::string FileName = dirname;
        FileName = dirname + projectpath;
        FileName += "pyramlogfile.cfg";
        InitializePZLOG(FileName);
    }
#endif
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream str;
        str << "\nRodando testes de piramede Hdiv" << std::endl;
        LOGPZ_DEBUG(logger,str.str())
    }
#endif
    
    // Parametros
    int dim = 2;
    int plevel = 4;
    int numthreads = 2;

    // Para 2D
    int nelx = 100;
    int nely = 100;

    // Para 3D
    int nref = 2;
    
    if (argc == 1) {
        cout << "\nATENCAO: voce nao passou argumentos, rodando c/ parametros hardcode!" << endl;
    }
    else if (argc == 4) {
        nref = atoi(argv[1]);
        plevel = atoi(argv[2]);
        numthreads = atoi(argv[3]);
        cout << "\nRodando com:" << endl;
        cout << "nref = " << nref << endl;
        cout << "plevel = " << plevel << endl;
        cout << "numthreads = " << numthreads << endl;
    }
    else{
        cout << "\nERRO - Num de argumento nao especificado" << endl;
        DebugStop();
    }
    
    if (dim == 2) {
        cout << "\n-> Setado para rodar problema 2D com poisson simples. Parametros:" << endl;
        cout << "nelX = " << nelx << endl;
        cout << "nely = " << nely << endl;
        cout << "plevel = " << plevel << endl;
        cout << "numthreads = " << numthreads << endl;
    }
    else if (dim == 3){
        cout << "\n-> Setado para rodar problema do cubo 3D com elasticidade lienar. Parametros:" << endl;
        cout << "nref = " << nref << endl;
        cout << "plevel = " << plevel << endl;
        cout << "numthreads = " << numthreads << endl;
    }
    
    // Malha Geometrica
    cout << "\nCriando a gmesh... "; cout.flush();
//    TTimer timer;
//    timer.start();
    TPZGeoMesh *gmesh = NULL;
    if (dim == 3){
        gmesh = MalhaCubo(projectpath,nref);
    }
    else if (dim == 2){
        gmesh = MalhaQuadrada(nelx,nely);
    }
    else{
        DebugStop();
    }

    if (0) {
        std::cout << "CUIDADO - Voce esta fazendo print da gmesh em vtk !!!!!!!!!!!!!!!" << std::endl;
        std::ofstream out("GeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }
//    timer.stop();
//    cout << timer.seconds() << " s" << endl;

    // Malha computacional
    cout << "\nCriando a cmesh... "; cout.flush();
//    timer.start();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(plevel);
    if (dim == 3) {
        InsertElasticityCubo(cmesh);
    }
    else if (dim == 2){
        InsertBidimensionalPoisson(cmesh,dim);
    }
    else{
        DebugStop();
    }

    cmesh->AutoBuild();
//    timer.stop();
//    cout << timer.seconds() << " s" << endl;
    std::cout << "Numero de equacoes = " << cmesh->NEquations() << std::endl;
    cout << "Num elements = " << cmesh->NElements() << endl;
    
    // Analysis
    cout << "\nCriando o analysis... "; cout.flush();
//    timer.start();
    TPZAnalysis an(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZSkylineStructMatrix skyl(cmesh);
    skyl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(skyl);
    an.SetSolver(step);
//    timer.stop();
//    cout << timer.seconds() << " s" << endl;
    
    // Resolvendo
    cout << "\nComecando o assemble... "; cout.flush();
//    timer.start();
    an.Assemble();
//    timer.stop();
    cout << "\nRodando com:" << endl;
    cout << "nref = " << nref << endl;
    cout << "plevel = " << plevel << endl;
    cout << "numthreads = " << numthreads << endl;
    std::cout << "Numero de equacoes = " << cmesh->NEquations() << std::endl;
//    cout << "Assemble executado em " <<  timer.seconds() << " s" << endl;
    
    
    std::string exePath(argv[0]);
    stringstream ss;
    ss << numthreads;
    string strnthreads = ss.str();
    string filename = exePath + ".txt";// + "_nthreads_" + strnthreads + ".txt";
    ofstream out;
    if (numthreads == 0) {
        out.open(filename.c_str());
    }else{
        out.open(filename.c_str(), ofstream::out | ofstream::app);
    }

    out << "\nRodando com:" << endl;
    out << "nref = " << nref << endl;
    out << "plevel = " << plevel << endl;
    out << "numthreads = " << numthreads << endl;
//    out << "T assemble = " << timer.seconds() << endl;
    
    return 0;
    
    an.Solve();
    // Pos Processamento
    TPZStack<string> scalnames, vecnames;
    if (dim == 3) {
        scalnames.Push("StressX");
        vecnames.Push("state");
    }
    else if (dim == 2){
        scalnames.Push("Solution");
    }

    string plotfile = "ResultMesh.vtk";
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an.PostProcess(0);
    
    std::cout << "FINISHED" << std::endl;
	return 0;
}



// ------------------------ Para testes do assemble -----------------------------

TPZGeoMesh *MalhaQuadrada(int &nelx, int &nely)
{
    const int bcid1 = -1, bcid2 = -2;
    TPZManVector<int,2> nelxy(2,0);
    nelxy[0] = nelx;
    nelxy[1] = nely;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.); // dois pontos nos 2 cantos do quadrado
    x1[0] = 1.;
    x1[1] = 1.;
    TPZGenGrid gengrid(nelxy,x0,x1);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh);
    gengrid.SetBC(gmesh,7,bcid1); // na direita
    gengrid.SetBC(gmesh,5,bcid2); // na esquerda

    gmesh->BuildConnectivity();
    return gmesh;
}

void InsertBidimensionalPoisson(TPZCompMesh *cmesh, int &dim)
{
    int matid = 1;
    const int bcidLeft = -1, bcidRight = -2;
    const int diri = 0;
    TPZMatPoisson3d *mat = new TPZMatPoisson3d(matid,dim);
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bcLeft = mat->CreateBC(mat, bcidLeft, diri, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    val2(0,0) = 2.;
    TPZBndCond *bcRight = mat->CreateBC(mat, bcidRight, diri, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
}

TPZGeoMesh *MalhaCubo(string &projectpath, const int &nref)
{
    long numnodes=-1;
    long numelements=-1;
    
    string FileName, dirname = PZSOURCEDIR;
    FileName = dirname + projectpath;
    FileName += "cube1.txt";
    
    {
        bool countnodes = false;
        bool countelements = false;
        
        ifstream read (FileName.c_str());
        
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            if(str == "Coordinates") countnodes = true;
            if(str == "end coordinates") countnodes = false;
            if(countnodes) numnodes++;
            
            if(str == "Elements") countelements = true;
            if(str == "end elements") countelements = false;
            if(countelements) numelements++;
        }
    }
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    TPZManVector <long> TopolTetra(4);
    
    const long Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    //setting nodes coords
    long nodeId = 0, elementId = 0, matElId = 1;
    
    ifstream read;
    read.open(FileName.c_str());
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    
    char buf[1024];
    read.getline(buf, 1024);
    read.getline(buf, 1024);
    std::string str(buf);
    long in;
    for(in=0; in<numnodes; in++)
    {
        read >> nodeId;
        read >> nodecoordX;
        read >> nodecoordY;
        read >> nodecoordZ;
        Node[nodeId-1].SetNodeId(nodeId);
        Node[nodeId-1].SetCoord(0,nodecoordX);
        Node[nodeId-1].SetCoord(1,nodecoordY);
        Node[nodeId-1].SetCoord(2,nodecoordZ);
        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
    }
    
    {
        read.close();
        read.open(FileName.c_str());
        
        long l , m = numnodes+5;
        for(l=0; l<m; l++)
        {
            read.getline(buf, 1024);
        }
        
        
        long el;
        int neumann1 = -4, neumann2 = -5;
        //std::set<int> ncoordz; //jeitoCaju
        for(el=0; el<numelements; el++)
        {
            read >> elementId;
            read >> TopolTetra[0]; //node 1
            read >> TopolTetra[1]; //node 2
            read >> TopolTetra[2]; //node 3
            read >> TopolTetra[3]; //node 4
            
            // O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
            TopolTetra[0]--;
            TopolTetra[1]--;
            TopolTetra[2]--;
            TopolTetra[3]--;
            
            long index = el;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
        }
        
        gMesh->BuildConnectivity();
        
        // Colocando as condicoes de contorno
        for(el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gMesh->ElementVec()[el];
            
            // na face x = 1
            TPZVec<long> ncoordzVec(0); long sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                long pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == 1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann1);
            }
            
            // Na face x = -1
            ncoordzVec.Resize(0);
            sizeOfVec = 0;
            for (int i = 0; i < 4; i++) 
            {
                long pos = tetra->NodeIndex(i);
                Nodefinder[i] = gMesh->NodeVec()[pos];
                
                Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[0] == -1.)
                {
                    sizeOfVec++;
                    ncoordzVec.Resize(sizeOfVec);
                    ncoordzVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordzVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordzVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,neumann2);	
            }
            
        }
        
        TPZVec <REAL> xyz(3,-1.), yz(3,-1.), z(3,1.);
        yz[0] = 1.;
        z[2] = -1;
        int bcidxyz = -1, bcidyz = -2, bcidz = -3;
        SetPointBC(gMesh, xyz, bcidxyz);
        SetPointBC(gMesh, yz, bcidyz);
        SetPointBC(gMesh, z, bcidz);
    }
    
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        const int nel = gMesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gMesh->Element(iel);
            if (gel || !gel->HasSubElement() || gel->Dimension() > 0) {
                gel->Divide(sons);
            }
        }
    }
    
    return gMesh;
}

/// Generate a boundary geometric element at the indicated node
void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc)
{
    // look for an element/corner node whose distance is close to start
    TPZGeoNode *gn1 = gr->FindNode(x);
    long iel;
    long nelem = gr->ElementVec().NElements();
    TPZGeoEl *gel;
    for (iel = 0; iel<nelem; iel++) {
        gel = gr->ElementVec()[iel];
        if(!gel) continue;
        int nc = gel->NCornerNodes();
        int c;
        for (c=0; c<nc; c++) {
            TPZGeoNode *gn = gel->NodePtr(c);
            if (gn == gn1) {
                break;
            }
        }
        if (c<nc) {
            TPZGeoElBC(gel, c, bc);
            return;
        }
    }
}

void InsertElasticityCubo(TPZCompMesh *mesh)
{
    mesh->SetDimModel(3);
    int nummat = 1, neumann = 1, mixed = 2;
    //	int dirichlet = 0;
    int dir1 = -1, dir2 = -2, dir3 = -3, neumann1 = -4., neumann2 = -5;   //, dirp2 = -6;
    TPZManVector<STATE> force(3,0.);
    //force[1] = 0.;
    
    STATE ElaE = 1000., poissonE = 0.2;   //, poissonV = 0.1, ElaV = 100.;
    
    STATE lambdaV = 0, muV = 0, alpha = 0, deltaT = 0;
    lambdaV = 11.3636;
    muV = 45.4545;
    alpha = 1.;
    deltaT = 0.01;
    
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat);
    //viscoelast->SetMaterialDataHooke(ElaE, poissonE, ElaV, poissonV, alpha, deltaT, force);
    //TPZViscoelastic *viscoelast = new TPZViscoelastic(nummat, ElaE, poissonE, lambdaV, muV, alphaT, force);
    TPZElasticity3D *viscoelast = new TPZElasticity3D(nummat, ElaE, poissonE, force);
    
    TPZFNMatrix<6> qsi(6,1,0.);
    //viscoelast->SetDefaultMem(qsi); //elast
    //int index = viscoelast->PushMemItem(); //elast
    TPZMaterial * viscoelastauto(viscoelast);
    mesh->InsertMaterialObject(viscoelastauto);
    
    // Neumann em x = 1;
    TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
    val2(0,0) = 1.;
    TPZBndCond *bc4 = viscoelast->CreateBC(viscoelastauto, neumann1, neumann, val1, val2);
    TPZMaterial * bcauto4(bc4);
    mesh->InsertMaterialObject(bcauto4);
    
    // Neumann em x = -1;
    val2(0,0) = -1.;
    TPZBndCond *bc5 = viscoelast->CreateBC(viscoelastauto, neumann2, neumann, val1, val2);
    TPZMaterial * bcauto5(bc5);
    mesh->InsertMaterialObject(bcauto5);
    
    val2.Zero();
    // Dirichlet em -1 -1 -1 xyz;
    val1(0,0) = 1e4;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc1 = viscoelast->CreateBC(viscoelastauto, dir1, mixed, val1, val2);
    TPZMaterial * bcauto1(bc1);
    mesh->InsertMaterialObject(bcauto1);
    
    // Dirichlet em 1 -1 -1 yz;
    val1(0,0) = 0.;
    val1(1,1) = 1e4;
    val1(2,2) = 1e4;
    TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
    TPZMaterial * bcauto2(bc2);
    mesh->InsertMaterialObject(bcauto2);
    
    // Dirichlet em 1 1 -1 z;
    val1(0,0) = 0.;
    val1(1,1) = 0.;
    val1(2,2) = 1e4;
    TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
    TPZMaterial * bcauto3(bc3);
    mesh->InsertMaterialObject(bcauto3);
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    // Re-constructing connectivities
//    gmesh->ResetConnectivities();
//    gmesh->BuildConnectivity();
}

void GroupElements(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Type() == EPiramide) {
            std::list<TPZOneShapeRestraint> shape;
            shape = cel->GetShapeRestraints();
            if (shape.size() != 1) {
                DebugStop();
            }
            TPZOneShapeRestraint local = *shape.begin();
            TPZMultiphysicsElement *mult = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if (!mult) {
                DebugStop();
            }
            TPZCompElHDiv<pzshape::TPZShapePiram> *hdiv = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapePiram> *>(mult->Element(0));
            if (!hdiv) {
                DebugStop();
            }
            int face = hdiv->RestrainedFace();
            TPZGeoElSide gelside(gel,face);
            TPZGeoElSide neighbour = gelside.Neighbour();
            TPZCompEl *celneigh = neighbour.Element()->Reference();
            TPZMultiphysicsElement *multneigh = dynamic_cast<TPZMultiphysicsElement *>(celneigh);
            if (!multneigh) {
                DebugStop();
            }
            TPZCompElHDiv<pzshape::TPZShapeTetra> *hdivneigh = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(multneigh->Element(0));
            if (!hdivneigh) {
                DebugStop();
            }
            long index;
            TPZElementGroup *grp = new TPZElementGroup(*cmesh,index);
            grp->AddElement(cel);
            grp->AddElement(celneigh);
        }
    }
}

void LoadSolution(TPZCompMesh *cpressure)
{
    long nel = cpressure->NElements();
    for (long iel=0; iel<nel; iel++) {
        TPZCompEl *cel = cpressure->Element(iel);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != 3) {
            continue;
        }
        if (gel->Type() == EPiramide) {
            TPZGeoNode *top = gel->NodePtr(4);
            TPZManVector<REAL,3> topco(3),valvec(1);
            top->GetCoordinates(topco);
            Forcing(topco, valvec);
            STATE topval = valvec[0];
            TPZConnect &c = cel->Connect(0);
            long seqnum = c.SequenceNumber();
            for (int i=0; i<4; i++) {
                cpressure->Block()(seqnum,0,i,0) = topval;
            }
            for (int i=0; i<4; i++) {
                TPZConnect &c = cel->Connect(i+1);
                TPZGeoNode *no = gel->NodePtr(i);
                no->GetCoordinates(topco);
                Forcing(topco, valvec);
                STATE nodeval = valvec[0];
                seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,0,0) = nodeval-topval;
            }
        }
        else if(gel->Type() == ETetraedro)
        {
            for (int i=0; i<4; i++) {
                TPZConnect &c = cel->Connect(i);
                TPZGeoNode *no = gel->NodePtr(i);
                TPZManVector<REAL,3> topco(3);
                no->GetCoordinates(topco);
                TPZManVector<STATE,3> valvec(1);
                Forcing(topco, valvec);
                STATE nodeval = valvec[0];
                long seqnum = c.SequenceNumber();
                cpressure->Block()(seqnum,0,0,0) = nodeval;
            }
            
        }
    }
}

void ProjectFlux(TPZCompMesh *cfluxmesh)
{
    TPZAnalysis an(cfluxmesh,false);
    TPZSkylineStructMatrix str(cfluxmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
}
