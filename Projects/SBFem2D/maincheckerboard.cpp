#include "Common.h"

void rect_mesh(int numEleVer, double vert_domainsize);

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void AnalyseSolutionOfHeterogeneousGroup(TPZCompMesh *cmesh, int POrder, int irefskeleton, REAL contrast, TPZVec<REAL> &scalingcenter);


int main(int argc, char *argv[])
{
    
    
    //tutorial();
    int numvert_elements = 50;
    double domain_vertsize = 25.;
    REAL contrast = 1./20.;
    TPZManVector<REAL,3> scalingcenter(3,0.);
    scalingcenter[0] = domain_vertsize/(2.*numvert_elements);
    scalingcenter[1] = domain_vertsize/(2.*numvert_elements);

    std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minnelx = 4;
    int maxnelx = 5;
    int minrefskeleton = 0;
    int maxrefskeleton = 3;
    int minporder = 1;
    int maxporder = 5;
    int counter = 1;
    for(int nelx = minnelx; nelx < maxnelx; nelx *=2)
    {
        scalingcenter[0] = domain_vertsize/(nelx);
        scalingcenter[1] = domain_vertsize/(nelx);

        for (int irefskeleton = minrefskeleton; irefskeleton < maxrefskeleton; irefskeleton++)
        {
            for ( int POrder = minporder; POrder < maxporder; POrder += 1)
            {
                rect_mesh(nelx, domain_vertsize);
                TPZCompMesh *SBFem = ReadJSonFile("rect.json",irefskeleton,POrder,contrast);
                
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
                    vecnames.Push("Displacement");
                    scalnames.Push("SigmaX");
                    scalnames.Push("SigmaY");
                    ElasticAnalysis->DefineGraphMesh(2, scalnames, vecnames, "../elasticityCheckerboard.vtk");
                    ElasticAnalysis->PostProcess(1);
                }
                
                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                //                ElasticAnalysis->PostProcessError(errors);
                
                AnalyseSolutionOfHeterogeneousGroup(SBFem, POrder, irefskeleton, contrast,scalingcenter);
                
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

#include "TPZSBFemElementGroup.h"

TPZSBFemElementGroup *FindHeterogeneous(TPZCompMesh *cmesh, TPZVec<REAL> &scalingcenter)
{
    long nel = cmesh->NElements();
    for(long el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!elgr) {
            continue;
        }
        REAL dist = 0.;
        TPZStack<TPZCompEl *,5> elements = elgr->GetElGroup();
        std::set<int> matids;
        for (int subel=0; subel<elements.size(); subel++) {
            TPZCompEl *subcel = elements[subel];
            TPZGeoEl *subgel = subcel->Reference();
            if (!subgel) {
                DebugStop();
            }
            int nnodes = subgel->NNodes();
            TPZManVector<REAL,3> x(3);
            subgel->Node(nnodes-1).GetCoordinates(x);
            for (int i=0; i<3; i++) {
                dist += (x[i]-scalingcenter[i])*(x[i]-scalingcenter[i]);
            }
            matids.insert(subgel->MaterialId());
        }
        if (matids.size() > 1 && dist < 1.e-3) {
            return elgr;
        }
    }
    return 0;
}

void AnalyseSolutionOfHeterogeneousGroup(TPZCompMesh *cmesh, int POrder, int irefskeleton, REAL contrast, TPZVec<REAL> &scalingcenter)
{
    TPZSBFemElementGroup *celgrp = FindHeterogeneous(cmesh, scalingcenter);
    long neq = cmesh->NEquations();
    if (!celgrp) {
        DebugStop();
    }
    std::stringstream sout;
    sout << "../CheckerBoardEigenvalues.txt";
    
    std::ofstream results(sout.str(),std::ios::app);
    results.precision(15);
    // for circular domain with contrast
    results << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << " contrast " << contrast << std::endl;
    
    if(1)
    {
        std::multimap<REAL,REAL> eigmap;
        TPZManVector<double> eigval = celgrp->EigenvaluesReal();
        TPZFMatrix<double> coef = celgrp->CoeficientsReal();
        for (int i=0; i<eigval.size(); i++) {
            eigmap.insert(std::pair<double,double>(eigval[i],coef(i,0)));
        }
        for (std::multimap<double, double>::reverse_iterator it = eigmap.rbegin(); it!=eigmap.rend(); it++) {
            results << it->first << "|" << it->second << " ";
        }
    }
    results << std::endl;
    results << celgrp->EigenValues() << std::endl;

}
