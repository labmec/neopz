#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzbfilestream.h"

#include "TPZGmshReader.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
//#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "mixedpoisson.h"
#include "pzelasmat.h"
#include "pzelasthybrid.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZVecL2.h"
#include "TPZMatLaplacianHybrid.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZMixedPoissonParabolic.h"



#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "TPZCompMeshTools.h"
#include "pzcondensedcompel.h"
#include "pzfunction.h"
#include "pzgraphmesh.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"
#include "pzvisualmatrix.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzcheckgeom.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedHybridMeshControl.h"

#include "TPZFracSet.h"
#include "TPZFracSimulation.h"

//#include "meshgen.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

struct TRunConfig
{
    int numHDivisions = 0;
    int pOrderInternal = 1;
    int numDivSkeleton = 0;
    int pOrderSkeleton = 1;
    int Hybridize = 0;
    int Condensed = 1;
    int LagrangeMult = 0;
    int newline = 0;
    
    /// number of equations when not condensing anything
    int64_t fGlobalSystemSize = -1;
    /// number of equations considering local condensation
    int64_t fGlobalSystemWithLocalCondensationSize = -1;
    /// number of equations of the global system
    int64_t fNumeq = -1;
    
    REAL fDeltaT = 1.;
    
    /// number of timesteps
    int64_t nTimeSteps = 10;
    
    std::ostream &InlinePrint(std::ostream &out)
    {
        out <<  " numHDiv " << numHDivisions << " porderInternal " << pOrderInternal << " numDivSkeleton " << numDivSkeleton
        << " porderSkeleton " << pOrderSkeleton << " Hybridize " << Hybridize << " Condensed " << Condensed << " LagrangeMult " << LagrangeMult
        << " sysnocondense " << fGlobalSystemSize << " syslocalcondense " << fGlobalSystemWithLocalCondensationSize << " neq " << fNumeq;
        return out;
    }
    std::ostream &MathematicaInlinePrint(std::ostream &out)
    {
        out  << " numHDiv, " << numHDivisions << " ,porderInternal, " << pOrderInternal << " ,numDivSkeleton, " << numDivSkeleton
        << " ,porderSkeleton, " << pOrderSkeleton << " ,Hybridize, " << Hybridize << " ,Condensed, " << Condensed << " ,LagrangeMult, " << LagrangeMult
        << " ,sysnocondense, " << fGlobalSystemSize << " ,syslocalcondense, " << fGlobalSystemWithLocalCondensationSize << " ,neq, " << fNumeq;
        return out;
    }
    
    std::ostream &ConfigPrint(std::ostream &out)
    {
        out  << "_HSkel" << numDivSkeleton << "_pSkel" << pOrderSkeleton << "_HDiv" << numHDivisions << "_pInt" << pOrderInternal;
        return out;
    }
};


void SolveProblem(TPZMHMeshControl &MHM, std::string prefix, TRunConfig config);


using namespace std;



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif



int main(int argc, char *argv[])
{
    TExceptionManager except;
    
    TRunConfig Configuration;
    
    /// h refinement is not implemented now
    Configuration.numHDivisions = 0;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 1;
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    Configuration.Hybridize = 0;
    Configuration.Condensed = 0;

    if (argc == 8)
    {
        std::cout << "Executing using command line arguments\n";
    }

    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    int dimension = 2;
    TPZAutoPointer<TPZMHMixedHybridMeshControl> MHM = new TPZMHMixedHybridMeshControl(dimension);
    
    TPZFracSimulation frac(MHM);
    
    std::stringstream MHMPref;
    MHMPref << "HorizontalFrac";

    frac.ReadDataFile(MHMPref.str());
    
    
    MHM->SwitchLagrangeMultiplierSign(true);

    
    MHM->SetInternalPOrder(Configuration.pOrderInternal);
    MHM->SetSkeletonPOrder(Configuration.pOrderSkeleton);
    MHM->DivideSkeletonElements(Configuration.numDivSkeleton);
    if (Configuration.Hybridize)
    {
        MHM->Hybridize();
    }
    bool substructure = (bool) Configuration.Condensed;
    MHM->BuildComputationalMesh(substructure);

    std::cout << "MHM Computational meshes created\n";
    std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
    
    
    std::string configuration;
    {
        std::stringstream sout;
        sout << "H" << Configuration.numHDivisions << "-P" << Configuration.pOrderInternal;
        configuration = sout.str();
    }
    if(Configuration.LagrangeMult)
    {
        MHMPref << "_Lagr";
    }
    if (Configuration.Hybridize) {
        MHMPref << "_Hybr";
    }
    // compute the MHM solution
    Configuration.fGlobalSystemWithLocalCondensationSize = MHM->fGlobalSystemWithLocalCondensationSize;
    Configuration.fGlobalSystemSize = MHM->fGlobalSystemSize;
    Configuration.fNumeq = MHM->fNumeq;
    
    std::set<int> matids;
    TPZManVector<STATE,3> fluxval;
    matids.clear();
//    matids.insert(frac.MaterialId("BCIN"));
//    fluxval = MHM->CMesh()->Integrate("State", matids);

    
    SolveProblem(MHM, MHMPref.str(), Configuration);
    
    matids.clear();
    matids.insert(frac.MaterialId("f1"));
    fluxval = MHM->CMesh()->Integrate("Flux", matids);
    std::cout << "Flux integrated over " << frac.MaterialId("f1") << " " << fluxval << std::endl;
    matids.clear();
    matids.insert(frac.MaterialId("f2"));
    fluxval = MHM->CMesh()->Integrate("Flux", matids);
    std::cout << "Flux integrated over f2 " << fluxval << std::endl;
    matids.clear();
//    matids.insert(35);
//    fluxval = MHM->CMesh()->Integrate("Flux", matids);
//    std::cout << "Flux integrated over 35 " << fluxval << std::endl;
    matids.clear();
    matids.insert(frac.MaterialId("BCIN"));
    fluxval = MHM->CMesh()->Integrate("State", matids);
    std::cout << "Flux through inflow " << fluxval << std::endl;
    matids.clear();
    matids.insert(frac.MaterialId("BCOUT"));
    fluxval = MHM->CMesh()->Integrate("State", matids);
    std::cout << "Flux through outflow " << fluxval << std::endl;
    matids.clear();
    matids.insert(frac.MaterialId("BCNOFLOW"));
    fluxval = MHM->CMesh()->Integrate("State", matids);
    std::cout << "Flux through noflow " << fluxval << std::endl;
    matids.clear();

    return 0;
}

void SolveProblem(TPZMHMeshControl &MHM, std::string prefix, TRunConfig config)
{
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(MHM.CMesh(),shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(MHM.CMesh().operator->());
    strmat.SetNumThreads(0);
    
#else
    TPZSkylineStructMatrix strmat(MHM.CMesh().operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    TPZManVector<TPZAutoPointer<TPZCompMesh>,3> meshvec(2);
    meshvec = MHM.GetMeshes();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, MHM.CMesh());
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        MHM.CMesh()->Print(out);
    }
#endif
    
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::string plotfile1,plotfile2;
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim1.vtk";
        plotfile1 = sout.str();
    }
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim2.vtk";
        plotfile2 = sout.str();
    }
    std::cout << "plotfiles " << plotfile1.c_str() << " " << plotfile2.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Pressure");
    vecnames.Push("Flux");
    vecnames.Push("Derivative");
    int dim = MHM.CMesh()->Dimension();
    //    an.DefineGraphMesh(cmesh->Dimension()-1, scalnames, vecnames, plotfile1);
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile2);
    int resolution = 0;
    //    an.PostProcess(resolution,cmesh->Dimension()-1);
    an.PostProcess(resolution,MHM.CMesh()->Dimension());
}

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveParabolic(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config)
{
    
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }
#endif
    
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::string plotfile1,plotfile2;
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim1.vtk";
        plotfile1 = sout.str();
    }
    {
        std::stringstream sout;
        sout << prefix << "Approx-";
        config.ConfigPrint(sout) << "_dim2.vtk";
        plotfile2 = sout.str();
    }
    std::cout << "plotfiles " << plotfile1.c_str() << " " << plotfile2.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmesh->Dimension()-1, scalnames, vecnames, plotfile1);
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile2);
    int resolution = 0;
    an.PostProcess(resolution,cmesh->Dimension()-1);
    an.SetStep(0);
    an.PostProcess(resolution,cmesh->Dimension());
    
    for (int is=0; is<config.nTimeSteps; is++) {
        an.AssembleResidual();
        an.Solve();
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
        an.PostProcess(resolution,cmesh->Dimension()-1);
        an.SetStep(is+1);
        an.PostProcess(resolution, cmesh->Dimension());
        
    }
}


