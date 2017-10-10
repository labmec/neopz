#ifdef HAVE_CONFIG_H
#include <config.h>
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
#include "pzintel.h"
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
    long fGlobalSystemSize = -1;
    /// number of equations considering local condensation
    long fGlobalSystemWithLocalCondensationSize = -1;
    /// number of equations of the global system
    long fNumeq = -1;
    
    REAL fDeltaT = 1.;
    
    /// number of timesteps
    long nTimeSteps = 10;
    
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

struct TPZTimeStepResults
{
    TPZStack<REAL> fDelt, fTime;
    
    std::map<std::string,TPZStack<TPZManVector<STATE,3> > > fTimeValues, fAccumulatedValues;
    
    void Print(std::ostream &out);
    
};

void SolveProblem(TPZMHMeshControl &MHM, std::string prefix, TRunConfig config);

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveParabolic(TPZFracSimulation &frac, std::string prefix, TRunConfig config, TPZTimeStepResults &results, bool vtk);

/// initialize the pressure field of the simulation
void InitializePressureField(TPZMHMeshControl &MHM, REAL pressure);

/// adjusts the timestep in all darcy material objects
void SetTimeStep(TPZMHMeshControl &MHM, REAL delt);

/// accumulates the post processing results
void PostProcess(TPZFracSimulation &frac, REAL delt, TPZTimeStepResults &results);

using namespace std;



#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif
#ifdef LOG4CXX
static LoggerPtr loggerRes(Logger::getLogger("pz.submesh.residual"));
#endif



int main(int argc, char *argv[])
{
//    TExceptionManager except;
    
    TRunConfig Configuration;
    
    /// h refinement is not implemented now
    Configuration.numHDivisions = 1;
    /// PolynomialOrder - p-order
    Configuration.pOrderInternal = 1;
    Configuration.pOrderSkeleton = 1;
    Configuration.numDivSkeleton = 0;
    Configuration.Hybridize = 0;
    Configuration.Condensed = 1;

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
    MHMPref << "../FracMeshes/Fracture22";

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
    TPZTimeStepResults results;

    if (frac.SimulationType() == 0)
    {
        SolveProblem(MHM, MHMPref.str(), Configuration);
        PostProcess(frac, 0., results);
    }
    else
    {
        InitializePressureField(MHM, frac.InitialPressure());
        int numseries = frac.fTimeSteps.size();
        for (int step = 0; step < numseries; step++) {
            REAL delt = frac.fTimeSteps[step].first;
            SetTimeStep(MHM, delt);
            Configuration.fDeltaT = delt;
            Configuration.nTimeSteps = frac.fTimeSteps[step].second;
            SolveParabolic(frac, MHMPref.str(), Configuration, results, false);
        }
    }
    results.Print(std::cout);

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
    an.DefineGraphMesh(dim-1, scalnames, vecnames, plotfile1);
    an.DefineGraphMesh(dim, scalnames, vecnames, plotfile2);
    int resolution = 0;
    an.PostProcess(resolution,dim-1);
    an.PostProcess(resolution,dim);
}

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
void SolveParabolic(TPZFracSimulation &frac, std::string prefix, TRunConfig config, TPZTimeStepResults &results, bool vtk)

/// Solve the problem composed of a multiphysics mesh composed of compmeshes - applies to MHM and MHM-H(div)
//void SolveParabolic(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, TRunConfig config)
{
    TPZAutoPointer<TPZMHMixedHybridMeshControl> MHM = frac.MHM();
    TPZAutoPointer<TPZCompMesh> cmesh = MHM->CMesh();
    //calculo solution
    bool shouldrenumber = true;
    TPZAnalysis an(cmesh,shouldrenumber);
    
    int dim = cmesh->Dimension();
    int resolution = 0;
    
    if(vtk && results.fDelt.size() == 0)
    {
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
        std::cout << "plotfiles \n" << plotfile1.c_str() << "\n" << plotfile2.c_str() << std::endl;
        TPZStack<std::string> scalnames,vecnames;
        int matid = *(MHM->fMaterialIds.begin());
        TPZMaterial *mat = cmesh->FindMaterial(matid);
        if (!mat) {
            DebugStop();
        }
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
        an.DefineGraphMesh(dim-1, scalnames, vecnames, plotfile1);
        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile2);
        an.SetStep(results.fDelt.size());
        an.PostProcess(resolution,dim-1);
        an.SetStep(results.fDelt.size());
        an.PostProcess(resolution,dim);
    }

    
    
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(8);
#endif
    strmat.SetNumThreads(8);

    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
#ifdef LOG4CXX
    if(loggerRes->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "CalcStiff Analysis Rhs \n";
        an.Rhs().Print(sout);
        LOGPZ_DEBUG(loggerRes, sout.str())
    }
#endif
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
    TPZManVector<TPZAutoPointer<TPZCompMesh> ,5> compmeshes = MHM->GetMeshes();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    if(0)
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

    if(vtk)
    {
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
        std::cout << "plotfiles \n" << plotfile1.c_str() << "\n" << plotfile2.c_str() << std::endl;
        TPZStack<std::string> scalnames,vecnames;
        int matid = *(MHM->fMaterialIds.begin());
        TPZMaterial *mat = cmesh->FindMaterial(matid);
        if (!mat) {
            DebugStop();
        }
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
        an.DefineGraphMesh(dim-1, scalnames, vecnames, plotfile1);
        an.DefineGraphMesh(dim, scalnames, vecnames, plotfile2);
        an.SetStep(results.fDelt.size()+1);
        an.PostProcess(resolution,dim-1);
        an.SetStep(results.fDelt.size()+1);
        an.PostProcess(resolution,dim);
    }
    PostProcess(frac, config.fDeltaT, results);
    
    for (int is=1; is<config.nTimeSteps; is++) {
        an.AssembleResidual();
#ifdef LOG4CXX
        if(loggerRes->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "AssembleResidual Analysis Rhs \n";
            an.Rhs().Print(sout);
            LOGPZ_DEBUG(loggerRes, sout.str())
        }
#endif
        an.Solve();
        if(vtk)
        {
            an.SetStep(results.fDelt.size()+1);
            an.PostProcess(resolution,dim);
            an.SetStep(results.fDelt.size()+1);
            an.PostProcess(resolution, dim-1);
        }
        PostProcess(frac, config.fDeltaT, results);
    }
}

/// initialize the pressure field of the simulation
void InitializePressureField(TPZMHMeshControl &MHM, REAL pressure)
{
    TPZAutoPointer<TPZCompMesh> pressuremesh = MHM.PressureMesh();
    int dim = pressuremesh->Dimension();
    long nel = pressuremesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            continue;
        }
        int nc = intel->NCornerConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = intel->Connect(ic);
            long seqnum = c.SequenceNumber();
            pressuremesh->Block()(seqnum,0,0,0) = pressure;
        }
    }
    TPZBuildMultiphysicsMesh::TransferFromMeshes(MHM.GetMeshes(), MHM.CMesh());

}

/// adjusts the timestep in all darcy material objects
void SetTimeStep(TPZMHMeshControl &MHM, REAL delt)
{
    for (auto it = MHM.fMaterialIds.begin(); it != MHM.fMaterialIds.end(); it++) {
        TPZMaterial *mat = MHM.CMesh()->FindMaterial(*it);
        TPZMixedPoissonParabolic *mix = dynamic_cast<TPZMixedPoissonParabolic *>(mat);
        if(mix)
        {
            mix->SetDeltaT(delt);
        }
    }
}

/// accumulates the post processing results
void PostProcess(TPZFracSimulation &frac, REAL delt, TPZTimeStepResults &results)
{
    int nacc = results.fDelt.size();
    results.fDelt.Push(delt);
    results.fTime.Push(delt);
    TPZAutoPointer<TPZMHMixedHybridMeshControl> MHM = frac.MHM();
    
    for(int ip=0; ip < frac.fPostProcnames.size(); ip++)
    {
        std::string matname = frac.fPostProcnames[ip].first;
        int matid = frac.MaterialId(matname);
        std::set<int> matids;
        matids.insert(matid);
        TPZManVector<STATE> result = MHM->CMesh()->Integrate(frac.fPostProcnames[ip].second, matids);
        std::cout << "Flux integrated over " << frac.fPostProcnames[ip].first << " " << result << std::endl;
        results.fTimeValues[matname].Push(result);
        for (int i=0; i<result.size(); i++) {
            result[i] *= delt;
        }
        results.fAccumulatedValues[matname].Push(result);
    }
    if (nacc > 0) {
        results.fTime[nacc] += results.fTime[nacc-1];
        for (auto it = results.fAccumulatedValues.begin(); it != results.fAccumulatedValues.end(); it++) {
            TPZVec<STATE> &accumulate = (it->second)[nacc];
            TPZVec<STATE> &prev = (it->second)[nacc-1];
            int nv = accumulate.size();
            for (int iv = 0; iv<nv; iv++) {
                accumulate[iv] += prev[iv];
            }
        }
    }
}

void TPZTimeStepResults::Print(std::ostream &out)
{
//    TPZStack<REAL> fDelt, fTime;
//
//    std::map<std::string,TPZStack<TPZManVector<STATE,3> > > fTimeValues, fAccumulatedValues;
    int nsteps = fDelt.size();
    if (nsteps == 1) {
        for (auto it = fTimeValues.begin(); it != fTimeValues.end(); it++) {
            out << "Object name " << it->first << " values " << (it->second)[0] << std::endl;
        }
    }
    else
    {
        for (auto it = fTimeValues.begin(); it != fTimeValues.end(); it++) {
            std::string name = it->first;
            TPZStack<TPZManVector<STATE,3> > &timeval = it->second;
            TPZStack<TPZManVector<STATE,3> > &accumval = fAccumulatedValues[name];
            for (int step =0; step<nsteps; step++) {
                out << "Name " << it->first << " Time " << fTime[step] << " Delt " << fDelt[step] << " Values " << timeval[step] <<
                " Accumulated values " << accumval[step] << std::endl;
            }
        }
    }
    
    
};
