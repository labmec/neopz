#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common.h"

#include "TPZSBFemElementGroup.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void SubstituteMaterialObjects(TPZCompMesh *cmesh);

void InitializeSolution(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;

    int maxnelxcount = 2;
    int numrefskeleton = 1;
    int maxporder = 2;
    int counter = 1;
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::ELoadedBeam;
#endif
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 2 << (nelxcount-1);
                bool useexact = true;
                
                TPZCompMesh *SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
#ifdef LOG4CXX
                if(logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    SBFem->Print(sout);
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * Analysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                Analysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                REAL delt = 1;
                REAL nsteps = 5;
                int numthreads = 0;
                SolveParabolicProblem(Analysis, delt, nsteps, numthreads);
                
                
                
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
                
                TPZManVector<STATE> errors(3,0.);
                
                long neq = SBFem->Solution().Rows();
                
                if(scalarproblem)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularSolution.vtk");
                    Analysis->PostProcess(3);
                }
                else
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    vecnames.Push("Displacement");
                    scalnames.Push("SigmaX");
                    scalnames.Push("TauXY");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularElasticity2DSolution.vtk");
                    Analysis->PostProcess(3);
                }

                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
                
                
                
                std::stringstream sout;
                sout << "../RegularSolution";
                if (scalarproblem) {
                    sout << "Scalar.txt";
                }
                else
                    sout << "Elastic2D.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(* nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << nelxcount << "]][[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                if(0)
                {
                    std::cout << "Plotting shape functions\n";
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    Analysis->ShowShape("Heterogeneous.vtk", eqindex);
                }
                
                delete Analysis;
                delete SBFem;
                //                exit(-1);
            }
            //            exit(-1);
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}






void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void InitializeSolution(TPZCompMesh *cmesh)
{
    // create a H1 projection material object
    // set the exact solution
    // project the solution
    
}

void SwitchComputationMode(TPZCompMesh *cmesh, TPZSBFemElementGroup::EComputationMode mode, REAL delt)
{
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(!elgr) continue;
        switch (mode) {
            case TPZSBFemElementGroup::EStiff:
                elgr->SetComputeStiff();
                break;
            case TPZSBFemElementGroup::EMass:
                elgr->SetComputeTimeDependent(delt);
                break;
            case TPZSBFemElementGroup::EOnlyMass:
                elgr->SetComputeOnlyMassMatrix();
                break;
            default:
                DebugStop();
                break;
        }
    }
}
