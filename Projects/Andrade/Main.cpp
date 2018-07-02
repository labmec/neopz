#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"

#include <pzgengrid.h>
#include "TPZMatElasticity2D.h"
#include "TPZInterfaceEl.h"
#include "pzdiscgal.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include <tpzarc3d.h>

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzl2projection.h"

// renumbering options
#include "TPZCutHillMcKee.h"
#include "TPZSloanRenumbering.h"
#include "pzsloan.h"
#include "TPZBoostGraph.h"
#include "pzmetis.h"

#include "pzvisualmatrix.h"
#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzelasmat.h"
#include "pzelast3d.h"
#include "pzplaca.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSpStructMatrix.h"
#include "pzlog.h"
#include <iostream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzfunction.h"
#include "TPZReadGIDGrid.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZExtendGridDimension.h"


#include <cmath>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif


enum MMaterials {ENoMat, EFront, ELateral, EBottom, EPlateFront, EPlateLateral, EBeamLateral, ELagrange, EMatConcrete, EMatSteel, EBeam};
/// Transverse grid
TPZAutoPointer<TPZGeoMesh> TransverseGrid(REAL halfplatewidth, REAL width, REAL depth, TPZVec<int> &nx, TPZVec<REAL> &dx);


/// Extend the 2d mesh to volumetric elements
TPZAutoPointer<TPZGeoMesh> ExtendGrid(TPZAutoPointer<TPZGeoMesh> mesh2d, REAL platelength, REAL totalextent, int nelements);

// Defintions of Implemented Methods
TPZCompMesh *ComputationalElasticityMesh3D(TPZGeoMesh *gmesh,int pOrder);

// Defintions of Implemented Methods
TPZCompMesh *ComputationalElasticityMesh2D(TPZAutoPointer<TPZGeoMesh>  gmesh,int pOrder);

void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter = 20);

//	This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

//	These are tools for spatial and polynomial refinement and Postprocess of solutions
void PostProcessElasticity(TPZAnalysis &an, std::string plotfile);
void UniformRefinement(TPZGeoMesh  *gMesh, int nh);
void UniformRefinement(TPZGeoMesh *gMesh, int nh, int MatId);
void RefinElemComp(TPZCompMesh  *cMesh, int indexEl);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void ReservoirPressure(const TPZVec<STATE> &x, TPZVec<STATE> &p,  TPZFMatrix<STATE> &gradp);

void OldTrial();

void CutByPlane(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &center, TPZVec<REAL> &dir);

void PutEnvelope(TPZGeoMesh *gmesh, int matid, REAL zco, REAL yco);

void CutByPlaneAbove(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &center, TPZVec<REAL> &dir, REAL zco);

void InsertBoundaryElements(TPZGeoMesh *gmesh);

void RefineGMesh(TPZGeoMesh *gmesh);

void SwitchBCTrilho2(TPZGeoMesh *gmesh);

void DisconnectElementLayer(TPZGeoMesh *gmesh);

void AdaptPOrders(TPZCompMesh *cmesh, REAL zsmall, int porder);

void AddZZeroFaces(TPZGeoMesh *gmesh);


int bcbottom = 8;
int bcsidex = 9;
int bcsidey = 10;
int bcloadtop = 11;
int matchapa = 1;
int matgraut = 2;
int matenchimento = 3;
int mattrilho1 = 4;
int mattrilho2 = 5;
int bctrilho2 = 51;
int bctrilho1 = 41;
int bcloadtopTESTE = -1;
int bctopsurface = 12;

REAL xplane1 = 1435;
REAL xplane2 = -4535;

int main(int argc, char *argv[])
{
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    gRefDBase.InitializeRefPatterns();
    
    TPZReadGIDGrid readGid;
    
    TPZGeoMesh *gmesh =  readGid.GeometricGIDMesh("../rail.dump");
    
    {
        std::ofstream out("../GmeshOrig.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    InsertBoundaryElements(gmesh);
    // CASO NAO TEM TRILHO 2
    SwitchBCTrilho2(gmesh);
    
    AddZZeroFaces(gmesh);

    RefineGMesh(gmesh);
    
    {
        std::ofstream out("../GmeshRef.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    
    const int porder = 2;
    TPZCompMesh * cmesh = ComputationalElasticityMesh3D(gmesh,porder);
    
    AdaptPOrders(cmesh, -400., 2);

    // toto
//    {
//        std::set<int> mats;
//        mats.insert(bcbottom);
//        TPZManVector<STATE> force(1,0.);
//        force = cmesh->Integrate("StressZ", mats);
//        std::cout << "force = " << force << std::endl;
//    }
    
    {
        std::ofstream out("../CmeshRef.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    }
    {
        TPZFMatrix<STATE> visualf;
        cmesh->ComputeFillIn(150, visualf);
        VisualMatrix(visualf,"../VisualMatrixBefore.vtk");
    }

    
    TPZAnalysis an;
    TPZAutoPointer<TPZRenumbering> renumber;
//    renumber = new TPZSloan;
//    renumber = new TPZCutHillMcKee;
    renumber = new TPZMetis;
    an.SetRenumber(renumber);
    an.SetCompMesh(cmesh, true);
    
#ifdef PZDEBUG
    {
        std::ofstream out("../gmesh.txt");
        gmesh->Print(out);
    }
    
    {
        std::ofstream out("../cmesh.txt");
        cmesh->Print(out);
    }
    
#endif
    {
        TPZFMatrix<STATE> visualf;
        cmesh->ComputeFillIn(150, visualf);
        VisualMatrix(visualf,"../MatrixMetis.vtk");
    }

    std::cout << "NEquations " << an.Solution().Rows() << std::endl;
    
    SolveSist(&an,cmesh);
    
    std::set<int> mats;
    mats.insert(bcbottom);
    TPZManVector<STATE> Force(1,0.);
    Force = cmesh->Integrate("StressZ", mats);
    std::cout << "Integrated sigma_Z " << Force << std::endl;
    
    std::cout << "Post processing" << std::endl;
    
    PostProcessElasticity(an, "../postProc.vtk");
    
    std::cout << "Finished\n";
//    delete cmesh;
//    delete gmesh;
    return EXIT_SUCCESS;
}

void AdjustBoundary(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Dimension() == 3 || gel->HasSubElement()) {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        bool should_refine = false;
        int nsub = -1;
        int numneigh = 0;
        while (gelside != neighbour) {
            nsub = neighbour.Element()->NSideSubElements(neighbour.Side());
            if (neighbour.Element()->HasSubElement() && nsub > 1) {
                should_refine = true;
            }
            numneigh++;
            neighbour = neighbour.Neighbour();
        }
        if (should_refine == true) {
            TPZAutoPointer<TPZRefPattern> match = TPZRefPatternTools::PerfectMatchRefPattern(gel);
            if (!match) {
                DebugStop();
            }
            gel->SetRefPattern(match);
            TPZStack<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
}

void InsertBoundaryElements(TPZGeoMesh *gmesh)
{
    REAL xplane = 1435;
    REAL yplane = 1161;
    REAL zplane = -1530;
    REAL ztop = 122.8;
    REAL minx = 0;
    REAL minz = 0;
    REAL maxz = 0;
    
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            int bccreated = 999;
            TPZGeoElSide gelside(gel,is);
            if (gel->SideDimension(is) == 1) {
                TPZManVector<REAL,3> xcenter(3);
                gelside.CenterX(xcenter);
                if (xcenter[2] > maxz) {
                    maxz = xcenter[2];
                }
                if (fabs(xcenter[0]) < 1 && fabs(xcenter[2]-ztop) < 1) {
                    bccreated = bcloadtop;
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    while (neighbour != gelside) {
                        if (neighbour.Element()->MaterialId() == bcloadtop) {
                            bccreated = 999;
                        }
                        neighbour = neighbour.Neighbour();
                        
                    }
                    if (bccreated == bcloadtop)
                    {
                        std::cout << "Added boundary bcloadtop xcenter = " << xcenter << std::endl;
                    }
                }
            }
            if (gel->SideDimension(is) == 2)
            {
                TPZManVector<REAL,3> xcenter(3);
                gelside.CenterX(xcenter);
                if (xcenter[0] < minx) {
                    minx = xcenter[0];
                }
                if (xcenter[2] < minz) {
                    minz = xcenter[2];
                }
                if (fabs(xcenter[2] - zplane) < 1.) {
                    bccreated = bcbottom;
                }
                if (fabs(fabs(xcenter[0])-xplane) < 1) {
                    bccreated = bcsidex;
                }
                if (fabs(fabs(xcenter[1])-yplane) < 1) {
                    bccreated = bcsidey;
                }
                if (fabs(xcenter[2]-ztop) <1) {
                    bccreated = bctopsurface;
                }
            }
            if(bccreated != 999)
            {
                if (gelside.Dimension() == 2 && gelside.Neighbour() != gelside) {
                    DebugStop();;
                }
                gel->CreateBCGeoEl(is, bccreated);
            }
        }
    }
    std::cout << "minx = " << minx << " minz = " << minz << " maxz " << maxz << std::endl;
}

void RefineGMesh(TPZGeoMesh *gmesh)
{
    TPZManVector<REAL,3> center(3,0.), dir(3,0.);
    center[1] = 170.;
    dir[1] = 1.;
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    center[1] = -170.;
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    center[1] = 0.;
    center[2] = -35.;
    dir[1] = 0;
    dir[2] = 1.;
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    CutByPlane(gmesh, 0, center, dir);
    
    dir[2] = 0;
    dir[0] = 1.;
    center[2] = 0.;
    for (int ndiv = 2; ndiv<16; ndiv *=2) {
        TPZManVector<REAL> xco(ndiv,0.);
        REAL delta = 740./ndiv;
        for (int i=0; i<ndiv; i++) {
            xco[i] = -370.+delta/2.+i*delta;
        }
        for (int i=0; i<ndiv; i++) {
            center[0] = xco[i];
            CutByPlane(gmesh, 0, center, dir);
        }
    }
    
    center[0] = 0.;
    dir[0] = 0.;
    center[2] = 0.;
    dir[2] = 0.;
    dir[1] = 1.;
    CutByPlaneAbove(gmesh, 0, center, dir, -400.);
    center[1] = 3.;
    CutByPlaneAbove(gmesh, 0, center, dir, -100.);
    center[1] = -3.;
    CutByPlaneAbove(gmesh, 0, center, dir, -100.);
    
    
    PutEnvelope(gmesh, 6, -30., 161.);
    
    for (int pass = 0; pass <2; pass++)
    {
        int64_t nel = gmesh->NElements();
        std::set<int> matids;
        matids.insert(6);
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (gel->HasSubElement()) {
                continue;
            }
            TPZRefPatternTools::RefineDirectional(gel, matids);
        }
    }
    
}


TPZCompMesh *ComputationalElasticityMesh2D(TPZAutoPointer<TPZGeoMesh>  gmesh,int pOrder)
{
    
    // remove some connectivities 3, 5
    TPZGeoEl *gel = gmesh->Element(0);
    TPZGeoElSide gelside(gel,3);
    gelside.RemoveConnectivity();
    gelside.SetSide(5);
    gelside.RemoveConnectivity();
    gelside.SetSide(4);
    TPZGeoElSide neighbour = gelside.NNeighbours();
    int matid = neighbour.Element()->MaterialId();
    gel->SetMaterialId(matid);
    neighbour.Element()->RemoveConnectivities();
    int64_t index = neighbour.Element()->Index();
    delete neighbour.Element();
    gmesh->ElementVec()[index] = 0;
    
    
    // Plane strain assumption
    int planestress = 0;
    
    // Getting mesh dimension
    int dim = 2;
    
    TPZMatElasticity2D *materialConcrete;
    materialConcrete = new TPZMatElasticity2D(EMatConcrete);
    
    TPZMatElasticity2D *materialSteel;
    materialSteel = new TPZMatElasticity2D(EMatSteel);
    
    
    // Setting up paremeters
    materialConcrete->SetfPlaneProblem(planestress);
    materialConcrete->SetElasticity(25.e6, 0.25);
    materialSteel->SetElasticity(205.e6, 0.25);
    
    //material->SetBiotAlpha(Alpha);cade o metodo?
    
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val1(1,1) = 1.e12;
    TPZMaterial * BCond2 = materialConcrete->CreateBC(materialConcrete,EBottom,3, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val1.Zero();
    val1(0,0) = 1.e12;
    TPZMaterial * BCond3 = materialConcrete->CreateBC(materialConcrete,ELateral,3, val1, val2);
    
    val2(0,0) = 0.0;
    val2(1,0) = -1000.0;
    val1.Zero();
    TPZMaterial * BCond4 = materialSteel->CreateBC(materialSteel,EBeam,1, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(materialConcrete);
    cmesh->InsertMaterialObject(materialSteel);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->AutoBuild();
    return cmesh;
    
}

//int bcbottom = 8;
//int bcsidex = 9;
//int bcsidey = 10;
//int bcloadtop = 11;
//int matchapa = 1;
//int matgraut = 2;
//int matenchimento = 3;
//int mattrilho1 = 4;
//int mattrilho2 = 5;

TPZCompMesh * ComputationalElasticityMesh3D(TPZGeoMesh *gmesh,int pOrder)
{
    // Getting mesh dimension
    const int dim = 3;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    {//material da chapa
        const REAL Ey = 205000.;
        const REAL poisson = 0.3;
        const int matid = matchapa;
        TPZManVector<STATE,3> fx(3,0.);
        TPZElasticity3D * mat = new TPZElasticity3D(matid,Ey,poisson,fx);
        mat->SetVonMises(300.);
        cmesh->InsertMaterialObject(mat);
    }
    
    {//material da trilho1
        const REAL Ey = 205000.;
        const REAL poisson = 0.3;
        const int matid = mattrilho1;
        TPZManVector<STATE,3> fx(3,0.);
        TPZElasticity3D * mat = new TPZElasticity3D(matid,Ey,poisson,fx);
        mat->SetVonMises(690.);
        cmesh->InsertMaterialObject(mat);
        //int bcsidex = 9;
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        val2(0,0) = 1.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bctrilho1, 3, val1, val2));
    }
    
    if(1)
    {//material da trilho2
        const REAL Ey = 205000.;
        const REAL poisson = 0.3;
        const int matid = mattrilho2;
        TPZManVector<STATE,3> fx(3,0.);
        TPZElasticity3D * mat = new TPZElasticity3D(matid,Ey,poisson,fx);
        mat->SetVonMises(690.);
        cmesh->InsertMaterialObject(mat);
        
        //int bcsidex = 9;
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
        val2(0,0) = 1.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bctrilho2, 3, val1, val2));
    }
    
    REAL percTracao = 0.1;
    {//material do concreto de 40 MPa
        const REAL Ey = 35417.;
        const REAL poisson = 0.2;
        const int matid = matgraut;
        TPZManVector<STATE,3> fx(3,0.);
        TPZElasticity3D * mat = new TPZElasticity3D(matid,Ey,poisson,fx);
        mat->SetMohrCoulomb(40.,percTracao*40.);
        cmesh->InsertMaterialObject(mat);
    }
    
    {//material do concreto de 30 MPa
        const REAL Ey = 27000;
        const REAL poisson = 0.2;
        const int matid = matenchimento;
        TPZManVector<STATE,3> fx(3,0.);
        TPZElasticity3D * mat = new TPZElasticity3D(matid,Ey,poisson,fx);
        mat->SetMohrCoulomb(30.,percTracao*30.);
        cmesh->InsertMaterialObject(mat);
        
        //c.c.
        //int bcbottom = 8;
        TPZFNMatrix<9,STATE> val1(3,3,0.), val2(3,1,0.);
//        val1(0,0) = 1.e-3;
//        val1(1,1) = 1.e-3;
//        val1(2,2) = 1.e12;
        val2(2) = 1.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcbottom, 3, val1, val2));
        val1.Zero();
        val2.Zero();
        
        //int bcsidex = 9;
        val2.Zero();
        val2(0,0) = 1.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcsidex, 3, val1, val2));
        
        //int bcsidey = 10;
        val1.Zero();
        val2.Zero();
        val2(1,0) = 1.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcsidey, 3, val1, val2));
        
        //int bcloadtop = 11;
        val2.Zero();
        val2(2,0) = -800000./120.;
        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcloadtop, 1, val1, val2));
        
//        somente para teste de tensao uniforme
//        cmesh->InsertMaterialObject(mat->CreateBC(mat, bcloadtopTESTE, 1, val1, val2));//toto
        
        
    }
    
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    return cmesh;
    
}

#include "TPZParFrontStructMatrix.h"
#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

#include "pzskylmat.h"

#define VTK
void SolveSist(TPZAnalysis *an, TPZCompMesh *Cmesh)
{
//    TPZParFrontStructMatrix<TPZFrontSym<STATE> > strmat(Cmesh);
    TPZSkylineStructMatrix strmat(Cmesh);
//    TPZSymetricSpStructMatrix strmat(Cmesh);
    strmat.SetNumThreads(8);
    an->SetStructuralMatrix(strmat);

    int64_t neq = Cmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an->SetSolver(step);

    an->Assemble();
    
//    std::ofstream andrade("../Andrade.mtx");
//    andrade.precision(16);
//    an->Solver().Matrix()->Print("Andrade",andrade,EMatrixMarket);
//    std::cout << "Leaving Assemble\n";
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
#endif
    
//#define NONO
#ifdef NONO
    step.SetMatrix(an->Solver().Matrix());

    TPZAutoPointer<TPZMatrix<STATE> > matrix = an->Solver().Matrix();
    TPZSkylMatrix<STATE> *skylmat = dynamic_cast<TPZSkylMatrix<STATE> *>(matrix.operator->());
    TPZSkylMatrix<float> * floatmat = new TPZSkylMatrix<float>;
    floatmat->CopyFrom(*skylmat);
    floatmat->Decompose_Cholesky();
    TPZSkylMatrix<STATE> *floatdec = new TPZSkylMatrix<STATE>;
    floatdec->CopyFrom(*floatmat);
    TPZStepSolver<STATE> stepfloat;
    stepfloat.SetMatrix(floatdec);
    stepfloat.SetDirect(ECholesky);

    step.SetCG(10, stepfloat, 1.e-6, 0);
    
    an->SetSolver(step);

#endif

    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    an->Solve(); 
    
#ifdef USING_BOOST
    boost::posix_time::ptime t3 = boost::posix_time::microsec_clock::local_time();
    std::cout << "Time for assembly " << t2-t1 << " Time for solving " << t3-t2 << std::endl;
#endif

    
}

void PostProcessElasticity(TPZAnalysis &an, std::string plotfile)
{
    TPZManVector<std::string,10> scalnames(0), vecnames(0);
    
    
    scalnames.Resize(6);
    vecnames.Resize(1);
    scalnames[0] = "StressX";
    scalnames[1] = "StressY";
    scalnames[2] = "StressZ";
    scalnames[3] = "PlasticFunction";
    scalnames[4] = "MaterialId";
    scalnames[5] = "POrder";
    vecnames[0]= "Displacement";
    const int dim = 3;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void UniformRefinement(TPZGeoMesh * gMesh, int nh, int MatId)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        int64_t n = gMesh->NElements();
        for ( int64_t i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1){
                if (gel->MaterialId()== MatId){
                    gel->Divide (filhos);
                }
            }
        }//for i
    }//ref
}

void RefinElemComp(TPZCompMesh  *cMesh, int indexEl)
{
    
    TPZVec<int64_t > subindex;
    int64_t nel = cMesh->ElementVec().NElements();
    for(int64_t el=0; el < nel; el++){
        TPZCompEl * compEl = cMesh->ElementVec()[el];
        if(!compEl) continue;
        int64_t ind = compEl->Index();
        if(ind==indexEl){
            compEl->Divide(indexEl, subindex, 1);
        }
    }
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
    
    TPZVec<int64_t > subindex;
    for (int64_t iref = 0; iref < ndiv; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
        int64_t nel = elvec.NElements();
        for(int64_t el=0; el < nel; el++){
            TPZCompEl * compEl = elvec[el];
            if(!compEl) continue;
            int64_t ind = compEl->Index();
            compEl->Divide(ind, subindex, 0);
        }
    }
}

void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter)
{
    int iter = 0;
    REAL error = 1.e10, NormResLambdaLast = 1.e10;;
    const REAL tol = 1.e-5;
    
    int numeq = an->Mesh()->NEquations();
    
    TPZFMatrix<STATE> Uatk0(an->Solution());
    TPZFMatrix<STATE> Uatk(Uatk0),DeltaU(Uatk0);
    if(Uatk0.Rows() != numeq) Uatk0.Redim(numeq,1);
    
    an->Assemble();
    an->Rhs() *= -1.0; //- [R(U0)];
    
    TPZAutoPointer< TPZMatrix<STATE> > matK; // getting X(Uatn)
    
    bool converged = false;
    while(!converged && iter < numiter) {
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            matK=an->Solver().Matrix();
            matK->Print("matK = ", sout,EMathematicaInput);
            an->Rhs().Print("Rhs = ", sout, EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        // Computing Uatk = Uatn + DeltaU;
        an->Solve();
        DeltaU= an->Solution();
        Uatk = Uatk0 + DeltaU;
        
        //Computing ||DeltaU||
        REAL NormOfDeltaU = Norm(DeltaU);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaU.Print("DeltaU = ", sout,EMathematicaInput);
            Uatk.Print("Uatk = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        an->LoadSolution(Uatk); // Loading Uatk
        an->Assemble();
        an->Rhs() *= -1.0; //- [R(U0)];
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            an->Rhs().Print("Res = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        // Computing ||[R(Uatk)]||
        double ResidualNorm = Norm(an->Rhs());
        double norm = NormOfDeltaU; //ResidualNorm;
        out << "Iteration n : " << (iter+1) << " : norms ||DeltaU|| e ||[R(Uatk)]|| : " << NormOfDeltaU << " / " << ResidualNorm << std::endl;
        
        if(norm < tol /*|| NormResLambda < tol*/) {
            out << "\nNewton Converged! Tolerance Of Norm(DeltaU) at n : " << (iter+1) << std::endl;
            out << "Norm ||DeltaU|| - USED : " << NormOfDeltaU << std::endl;
            out << "Norm ||[R(Uatk)]||  : " << ResidualNorm << std::endl;
            converged = true;
        }
        else if( (ResidualNorm - NormResLambdaLast) > 1.e-4 ) {
            out << "\nDivergent Method\n" << "Implement Line Search Please!!!!!!!!!!!!!!!!!!!!" << std::endl;
        }
        
        NormResLambdaLast = ResidualNorm;
        error = norm;
        iter++;
        Uatk0 = Uatk;
        out.flush();
    }
    
    if (error > tol) {
        DebugStop(); // Something is very wrong
    }
    
}

void ReservoirPressure(const TPZVec<STATE> &x, TPZVec<STATE> &p,  TPZFMatrix<STATE> &gradp)
{
    p[0] = 1.0e7;
}

/// Transverse grid
TPZAutoPointer<TPZGeoMesh> TransverseGrid(REAL halfplatewidth, REAL width, REAL depth, TPZVec<int> &nx, TPZVec<REAL> &dx)
{
    TPZManVector<REAL> x0(3,0.), x1(3,0.);
    x1[0] = width;
    x1[1] = -depth;
    TPZGenGrid gen(nx, x0, x1);
    TPZManVector<REAL,2> geometric_prog(2);
    geometric_prog[0] = gen.GeometricProgression(dx[0], width, nx[0]);
    geometric_prog[1] = gen.GeometricProgression(dx[1], width, nx[1]);
    gen.SetGeometricProgression(geometric_prog);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    gen.Read(gmesh);
    gen.SetBC(gmesh, 6, EBottom);
    gen.SetBC(gmesh, 5, ELateral);
    gen.SetBC(gmesh, 7, ELateral);
    TPZManVector<REAL> first(3,0.), last(3,0.);
    last[0] = halfplatewidth;
    gen.SetBC(gmesh, first, last, EMatSteel);
    gen.SetPointBC(gmesh.operator->(), first, EBeam);
    
    std::ofstream out("../transverse.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), out);
    return gmesh;
}

/// Extend the 2d mesh to volumetric elements
TPZAutoPointer<TPZGeoMesh> ExtendGrid(TPZAutoPointer<TPZGeoMesh> mesh2d, REAL platelength, REAL totalextent, int nelements)
{
    REAL layerthickness = totalextent/nelements;
    TPZExtendGridDimension extend(mesh2d, layerthickness);
    /**
     * @brief It reads the mesh since the archive of entrance finemesh, or since the fFineGeoMesh
     * passed in the constructor, and returns extended mesh. \n
     * The extension is from n (=1,2) dimensional mesh to (n+1) (=2,3) dimensional.
     * @param numlayers Numbers of layers to be incremented.
     * @param matidbottom Material id to bottom boundary surface after to extrude process.
     * @param matidtop Material id to top boundary surface after to extrude process.
     */
    TPZGeoMesh* newgeo = extend.ExtendedMesh(nelements, EFront, EFront);
    int64_t nel = newgeo->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = newgeo->Element(el);
        if (gel->MaterialId() != EMatSteel) {
            continue;
        }
        TPZManVector<REAL,3> xi(gel->Dimension(),0.),x(3,0.);
        gel->CenterPoint(gel->NSides()-1, xi);
        gel->X(xi, x);
        if (x[2] > platelength) {
            newgeo->DeleteElement(gel);
        }
        else
            if (gel->Dimension() == 1) {
                newgeo->DeleteElement(gel);
            }
    }
    
    std::ofstream out("../volumemesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(newgeo, out);
    
    return newgeo;
}

void OldTrial()
{
    
    TPZManVector<int,2> nx(2,8);
    REAL nelz = 10;
    REAL totalwidth = 5.;
    REAL totaldepth = 5.;
    REAL halfplatewith = 0.2;
    REAL platelength = 0.5;
    REAL platethickness = 0.01;
    REAL totallength = 6;
    TPZManVector<REAL,2> dx(2,0.);
    dx[0] = halfplatewith;
    // this is to run the 2d problem
    dx[1] = platethickness;
    
    TPZAutoPointer<TPZGeoMesh> gmesh2d = TransverseGrid(halfplatewith, totalwidth, totaldepth, nx,dx);
    
    /// Extend the 2d mesh to volumetric elements
    TPZAutoPointer<TPZGeoMesh> gmesh3d = ExtendGrid(gmesh2d, platelength, totallength, nelz);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMesh.txt");
        gmesh2d->Print(argument);
        std::ofstream Dummyfile("GeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2d,Dummyfile, true);
    }
    
    int Href = 1;
    int PElasticity = 2;
    //	UniformRefinement(gmesh2d, Href);
    
#ifdef LOG4CXX
    {
        //	Print Geometrical refined Base Mesh
        std::ofstream argument("RefinedGeometricMesh.txt");
        gmesh2d->Print(argument);
        std::ofstream Dummyfile("RefinedGeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh2d,Dummyfile, true);
    }
#endif
    
    TPZCompMesh * ComputationalMeshElasticity2D = ComputationalElasticityMesh2D(gmesh2d, PElasticity);
    
    
    TPZCompMesh * ComputationalMeshElasticity3D = ComputationalElasticityMesh3D(gmesh3d.operator->(), PElasticity);
    //	Print First computational mesh
    std::ofstream ArgumentElasticity2D("../ComputationalMeshForElasticity2D.txt");
    ComputationalMeshElasticity2D->Print(ArgumentElasticity2D);
    
    
    // Visualization of computational meshes
    bool mustOptimizeBandwidth = false;
    TPZAnalysis * ElasticAnalysis = new TPZAnalysis(ComputationalMeshElasticity2D,mustOptimizeBandwidth);
    
    SolveSist(ElasticAnalysis, ComputationalMeshElasticity2D);
    
    std::string ElasticityOutput;
    ElasticityOutput = "ComputationalMeshElasticity2D";
    std::stringstream ElasticityOutputfiletemp;
    ElasticityOutputfiletemp << ElasticityOutput << ".vtk";
    std::string ElasticityPlotfile = ElasticityOutputfiletemp.str();
    PostProcessElasticity(*ElasticAnalysis, ElasticityPlotfile);
    
    
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    
}

void CutByPlane(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &center, TPZVec<REAL> &dir)
{
    // identify all elements that have 4 ribs cut by the plane
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el< nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZManVector<int> sidescut(nsides,0);
        int nsidescut = 0;
        
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZManVector<REAL,3> point1(3), point2(3);
            TPZGeoNode *node1 = gel->SideNodePtr(is, 0);
            TPZGeoNode *node2 = gel->SideNodePtr(is, 1);
            node1->GetCoordinates(point1);
            node2->GetCoordinates(point2);
            REAL fac1 = 0., fac2 = 0.;
            for (int i=0; i<3; i++) {
                fac1 += (point1[i]-center[i])*dir[i];
                fac2 += (point2[i]-center[i])*dir[i];
            }
            if (fac1 * fac2 < 0.) {
                sidescut[is] = 1;
                nsidescut++;
            }
        }
        if (nsidescut == 4) {
            TPZAutoPointer<TPZRefPattern> match = TPZRefPatternTools::PerfectMatchRefPattern(gel, sidescut);
            gel->SetRefPattern(match);
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    AdjustBoundary(gmesh);
}

void PutEnvelope(TPZGeoMesh *gmesh, int matid, REAL zco, REAL yco)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() || gel->Dimension() != 3) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != 2) {
                continue;
            }
            TPZGeoElSide gelside(gel,is);
            TPZManVector<REAL,3> centerx(3);
            gelside.CenterX(centerx);
            bool candidate = false;
            if (fabs(centerx[2]-zco) < 1. && fabs(centerx[1]) < yco) {
                candidate = true;
            }
            if (fabs(fabs(centerx[1])-yco) < 1. &&  centerx[2] > -30.) {
                candidate = true;
            }
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == matid) {
                    candidate = false;
                }
                neighbour = neighbour.Neighbour();
            }
            if (candidate) {
                gel->CreateBCGeoEl(is, matid);
            }
        }
    }
    AdjustBoundary(gmesh);
}

void CutByPlaneAbove(TPZGeoMesh *gmesh, int matid, TPZVec<REAL> &center, TPZVec<REAL> &dir, REAL zco)
{
    // identify all elements that have 4 ribs cut by the plane
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el< nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide vol(gel,nsides-1);
        TPZManVector<REAL,3> centergel(3);
        vol.CenterX(centergel);
        if (centergel[2] < zco) {
            continue;
        }
        TPZManVector<int> sidescut(nsides,0);
        int nsidescut = 0;
        
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZManVector<REAL,3> point1(3), point2(3);
            TPZGeoNode *node1 = gel->SideNodePtr(is, 0);
            TPZGeoNode *node2 = gel->SideNodePtr(is, 1);
            node1->GetCoordinates(point1);
            node2->GetCoordinates(point2);
            REAL fac1 = 0., fac2 = 0.;
            for (int i=0; i<3; i++) {
                fac1 += (point1[i]-center[i])*dir[i];
                fac2 += (point2[i]-center[i])*dir[i];
            }
            if (fac1 * fac2 < 0.) {
                sidescut[is] = 1;
                nsidescut++;
            }
        }
        if (nsidescut == 4) {
            TPZAutoPointer<TPZRefPattern> match = TPZRefPatternTools::PerfectMatchRefPattern(gel, sidescut);
            gel->SetRefPattern(match);
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    //    return;
    nel = gmesh->NElements();
    for (int64_t el=0; el< nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            continue;
        }
        TPZAutoPointer<TPZRefPattern> match = TPZRefPatternTools::PerfectMatchRefPattern(gel);
        if (match) {
            gel->SetRefPattern(match);
            TPZManVector<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
    AdjustBoundary(gmesh);
}

void SwitchBCTrilho2(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            DebugStop();
        }
        if (gel->MaterialId() == bcsidex) {
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour(gelside.Neighbour());
            while (neighbour != gelside) {
                if (neighbour.Element()->MaterialId() == mattrilho2) {
                    gel->SetMaterialId(bctrilho2);
                    break;
                }
                if (neighbour.Element()->MaterialId() == mattrilho1) {
                    gel->SetMaterialId(bctrilho1);
                    break;
                }
                neighbour= neighbour.Neighbour();
            }
        }
    }
}

void AdaptPOrders(TPZCompMesh *cmesh, REAL zsmall, int porder)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZManVector<REAL,3> xcenter(3);
        gelside.CenterX(xcenter);
        if (xcenter[2] < zsmall) {
            intel->PRefine(porder);
        }
    }
    cmesh->ExpandSolution();
}

void AddZZeroFaces(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement()) {
            DebugStop();
        }
        if (gel->Dimension() != 3) {
            continue;
        }
        if (gel->MaterialId() != matgraut && gel->MaterialId() != matenchimento) {
            continue;
        }
        for (int iface = 8; iface < 26; iface++) {
            TPZGeoElSide gelface(gel,iface);
            if (gelface.Dimension() != 2) {
                continue;
            }
            TPZManVector<REAL,3> centerx(3);
            gelface.CenterX(centerx);
            if (fabs(centerx[2]) < 1. && fabs(centerx[0]) < 200 && fabs(centerx[1]) < 136) {
                //                TPZGeoElSide neighbour = gelface.Neighbour();
                std::cout << "el " << el << " face " << iface << " center = " << centerx << std::endl;
                gel->CreateBCGeoEl(iface, bcloadtopTESTE);
            }
        }
        
    }
}
