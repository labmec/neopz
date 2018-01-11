#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common3D.h"
#include "TPZSBFemElementGroup.h"
#include "TPZSBFemVolume.h"
#include "pzaxestools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
void AnalyseSolution(TPZCompMesh *cmesh);
#endif

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    int minnelx = 1;
    int maxnelx = 2;
    int minrefskeleton = 0;
    int numrefskeleton = 1;
    int minporder = 4;
    int maxporder = 5;
    int counter = 1;
	int numthreads = 1;

#ifdef _AUTODIFF
    ExactElast.fProblemType = TElasticity3DAnalytic::ETestShearMoment;
    ExactElast.fE = 1.;
    ExactElast.fPoisson = 0.;
#endif
    
    for ( int POrder = minporder; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = minrefskeleton; irefskeleton < numrefskeleton; irefskeleton++)
        {
            for(int nelx = minnelx; nelx < maxnelx; nelx *=2)
            {
            
                bool elast = true;
                TPZCompMesh *SBFem = SetupSquareMesh3D(nelx,irefskeleton,POrder, elast);
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
                SolveSist(Analysis, SBFem, numthreads);
                
                
//                AnalyseSolution(SBFem);
                
                std::cout << "Post processing\n";
                //        ElasticAnalysis->Solution().Print("Solution");
                //        mphysics->Solution().Print("expandec");
                
#ifdef _AUTODIFF
                Analysis->SetExact(Elasticity_exact);

                //                ElasticAnalysis->SetExact(Singular_exact);
#endif

                
                long neq = SBFem->Solution().Rows();
                
                std::string vtkfilename;
                if (elast) {
                    vtkfilename = "../Elast3DSolution.vtk";
                }
                else
                {
                    vtkfilename = "../Scalar3DSolution.vtk";
                }
                
                if(1)
                {
                    TPZStack<std::string> vecnames,scalnames;
                    // scalar
                    if(elast)
                    {
                        vecnames.Push("State");
                        scalnames.Push("StressX");
                        scalnames.Push("StressY");
                        scalnames.Push("StressZ");
                    } else
                    {
                        scalnames.Push("State");
                    }
                    Analysis->DefineGraphMesh(3, scalnames, vecnames, vtkfilename);
                    Analysis->PostProcess(3);
                }
                
                if(1)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                TPZManVector<REAL> errors(3,0.);

                Analysis->PostProcessError(errors);
                
                
                
                std::stringstream sout;
                if (elast) {
                    
                    sout << "../Elast3DSolution.txt";
                }
                else
                {
                    sout << "../Scalar3DSolution.txt";
                }
                
                {
                    std::ofstream results(sout.str(),std::ios::app);
                    results.precision(15);
                    results << "nx " << nelx << " numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << std::endl;
                    TPZFMatrix<double> errmat(1,3);
                    for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                    std::stringstream varname;
                    varname << "Errmat_" << nelx << "_" << irefskeleton << "_" << POrder << " = (1/1000000)*";
                    errmat.Print(varname.str().c_str(),results,EMathematicaInput);
                }
                std::cout << "Plotting shape functions\n";
                if(0 && nelx==minnelx && POrder == minporder && irefskeleton == minrefskeleton)
                {
                    std::string vtkfilename;
                    if (elast) {
                        vtkfilename = "../Elast3DShape.vtk";
                    }
                    else
                    {
                        vtkfilename = "../Scalar3DShape.vtk";
                    }
                    int numshape = 25;
                    if (numshape > SBFem->NEquations()) {
                        numshape = SBFem->NEquations();
                    }
                    TPZVec<long> eqindex(numshape);
                    for (int i=0; i<numshape; i++) {
                        eqindex[i] = i;
                    }
                    Analysis->ShowShape(vtkfilename, eqindex);
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

long SBFemGroup(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *grp = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(grp) return el;
    }
    return -1;
}

#ifdef _AUTODIFF
void AnalyseSolution(TPZCompMesh *cmesh)
{
    long el = SBFemGroup(cmesh);
    TPZSBFemElementGroup *elgrp = dynamic_cast<TPZSBFemElementGroup *>(cmesh->Element(el));
    auto subels = elgrp->GetElGroup();
    int nsub = subels.size();
    for (int isub = 0; isub < nsub; isub++) {
        TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(subels[isub]);
        TPZGeoEl *gel = vol->Reference();
        std::cout << "\n\n\n*********** ELEMENT " << isub << " ****************\n\n\n\n";
        for (int i=-1; i<2; i+=2) {
            for (int j=-1; j<2; j+=2) {
                TPZManVector<REAL,3> x(3,-1.), xco(3);
                x[0] = i;
                x[1] = j;
                gel->X(x, xco);
                TPZManVector<STATE,3> solex(3);
                TPZFNMatrix<9,STATE> dsolex(3,3);
                ExactElast.Solution(xco, solex, dsolex);
                
                TPZSolVec sol;
                TPZGradSolVec dsolax;
                TPZFNMatrix<9,REAL> axes(3,3);
                TPZFNMatrix<9,STATE> dsol(3,3), diff(3,3);
                vol->ComputeSolution(x, sol, dsolax, axes);
//                static void Axes2XYZ(const TPZFMatrix<TVar> &dudaxes, TPZFMatrix<TVar> &dudx, const TPZFMatrix<REAL> &axesv, bool colMajor = true){
                TPZAxesTools<STATE>::Axes2XYZ(dsolax[0], dsol, axes);
                diff = dsol-dsolex;
                REAL err = Norm(diff);
                if (err > 1.e-8)
                {
                    std::cout << "xco = " << x << " sol fem " << sol[0] << " dsol fem " << dsol << " axes " << axes << std::endl;
                    std::cout << "xco = " << x << " sol exa " << solex <<  " dsol exa " << dsolex << std::endl;
                    std::cout << "diff " << diff << std::endl;
                }
            }
        }
    }
}
#endif
