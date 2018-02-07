#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common3D.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int maxnelx = 25;
    int numrefskeleton = 2;
    int maxporder = 5;
    int counter = 1;
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            for(int nelx = 1; nelx < maxnelx; nelx *=2)
            {
            
                TPZCompMesh *SBFem = SetupSquareMesh3D(nelx,irefskeleton,POrder, false);
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
                SolveSist(Analysis, SBFem);
                
                
                
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
                
                Analysis->SetExact(Harmonic_exact);
                //                ElasticAnalysis->SetExact(Singular_exact);
                
                TPZManVector<STATE> errors(3,0.);
                
                long neq = SBFem->Solution().Rows();
                
                if(1)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    scalnames.Push("State");
                    Analysis->DefineGraphMesh(3, scalnames, vecnames, "../Scalar3DSolution.vtk");
                    Analysis->PostProcess(3);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                Analysis->PostProcessError(errors);
                
                
                
                std::stringstream sout;
                sout << "../Scalar3DSolution.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat_" << nelx << "_" << irefskeleton << "_" << POrder << " = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                std::cout << "Plotting shape functions\n";
                if(1 && nelx==1 && POrder == 1 && irefskeleton == 0)
                {
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    Analysis->ShowShape("Scalar3D.vtk", eqindex);
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


