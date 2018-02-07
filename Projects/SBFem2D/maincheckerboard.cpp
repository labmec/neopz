#include "Common.h"

void rect_mesh(int numEleVer, double vert_domainsize);

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


int main(int argc, char *argv[])
{
    
    
    //tutorial();
    int numvert_elements = 5;
    double domain_vertsize = 25.;
    
    
    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minnelx = 2;
    int maxnelx = 3;
    int minrefskeleton = 2;
    int maxrefskeleton = 3;
    int minporder = 3;
    int maxporder = 4;
    int counter = 1;
    for(int nelx = minnelx; nelx < maxnelx; nelx *=2)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for ( int POrder = minporder; POrder < maxporder; POrder += 1)
            {
                rect_mesh(nelx, domain_vertsize);
                TPZCompMesh *SBFem = ReadJSonFile("rect.json",irefskeleton,POrder);
                
                std::cout << "nelx = " << nelx << std::endl;
                std::cout << "irefskeleton = " << irefskeleton << std::endl;
                std::cout << "POrder = " << POrder << std::endl;
                
                // Visualization of computational meshes
                bool mustOptimizeBandwidth = true;
                TPZAnalysis * ElasticAnalysis = new TPZAnalysis(SBFem,mustOptimizeBandwidth);
                ElasticAnalysis->SetStep(counter++);
                std::cout << "neq = " << SBFem->NEquations() << std::endl;
                SolveSist(ElasticAnalysis, SBFem);
                
                std::cout << "Post processing\n";
                
                TPZManVector<STATE> errors(3,0.);
                
                long neq = SBFem->Solution().Rows();
                
                std::cout << "Neq = " << neq << std::endl;
                
                if(1)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    vecnames.Push("State");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../elasticityRegular.vtk");
                    ElasticAnalysis->PostProcess(3);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                //                ElasticAnalysis->PostProcessError(errors);
                
                
                
//                std::stringstream sout;
//                sout << "../CheckerboardDiagnostic.txt";
                
//                std::ofstream results(sout.str(),std::ios::app);
                //                results.precision(15);
                //                results << "nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << std::endl;
                //                TPZFMatrix<double> errmat(1,3);
                //                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                //                std::stringstream varname;
                //                varname << "Errmat_" << POrder << "_" << irefskeleton << " = (1/1000000)*";
                //                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                
                std::cout << "Plotting shape functions\n";
                if(0 && irefskeleton == 0)
                {
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    ElasticAnalysis->ShowShape("ElasticShape.vtk", eqindex);
                }
                
                delete ElasticAnalysis;
                delete SBFem;
                //                exit(-1);
            }
            //            exit(-1);
        }
    }
    std::cout << "Check:: Calculation finished successfully" << std::endl;
    return EXIT_SUCCESS;
}

