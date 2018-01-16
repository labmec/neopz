#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Common.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = false;

    int maxnelxcount = 7;
    int numrefskeleton = 2;
    int maxporder = 2;
    int counter = 1;
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::EStretchx;
#endif
    for ( int POrder = 1; POrder < 4; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
            if (POrder == 3 && !scalarproblem) {
                maxnelxcount = 3;
            }
            for(int nelxcount = 1; nelxcount < maxnelxcount; nelxcount += 1)
            {
                int nelx = 1 << (nelxcount-1);
                bool useexact = true;
                if(!scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.fE = 10;
                    ElastExact.fPoisson = 0.3;
                    ElastExact.fPlaneStress = 0;
#endif
                }
                TPZCompMesh *SBFem = SetupSquareMesh(nelx,irefskeleton,POrder, scalarproblem,useexact);
                if(0 && !scalarproblem)
                {
#ifdef _AUTODIFF
                    ElastExact.fProblemType = TElasticity2DAnalytic::EBend;
                    TPZManVector<REAL,3> x(3,0.);
                    TPZFNMatrix<4,STATE> tensor(2,2);
                    for(int i=-1; i<3; i+=2)
                    {
                        for (int j=-1; j<3; j+=2) {
                                x[0] = i;
                                x[1] = j;
                                ElastExact.Sigma(x, tensor);
                                std::cout << "x = " << x << " tensor " << tensor << std::endl;
                        }
                    }
#endif
				}
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
#ifdef _AUTODIFF
                if(scalarproblem)
                {
                    Analysis->SetExact(Harmonic_exact);
                }
                else
                {
                    Analysis->SetExact(Elasticity_exact);
                }
#endif
                //                ElasticAnalysis->SetExact(Singular_exact);
                
                
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
                    scalnames.Push("SigmaY");
                    scalnames.Push("TauXY");
                    scalnames.Push("EpsX");
                    scalnames.Push("EpsY");
                    scalnames.Push("EpsXY");
                    Analysis->DefineGraphMesh(2, scalnames, vecnames, "../RegularElasticity2DSolution.vtk");
                    Analysis->PostProcess(3);
                }

                if(0)
                {
                    std::ofstream out("../CompMeshWithSol.txt");
                    SBFem->Print(out);
                }
                
                std::cout << "Compute errors\n";
                
                TPZManVector<REAL> errors(3,0.);
                Analysis->PostProcessError(errors);
                
//                VerifyShapeFunctionIntegrity(Analysis->Mesh());
                
//                IntegrateDirect(Analysis->Mesh());
                
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

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"

void IntegrateDirect(TPZCompMesh *cmesh)
{
    long nel = cmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (elgr) {
            TPZStack<TPZCompEl *,5> elstack = elgr->GetElGroup();
            int nvol = elstack.size();
            TPZElementMatrix ekvol, efvol, ekgrp, efgrp;
            elgr->CalcStiff(ekgrp, efgrp);
            for (int iv=0; iv<nvol; iv++) {
                TPZCompEl *vcel = elstack[iv];
                TPZSBFemVolume *elvol = dynamic_cast<TPZSBFemVolume *>(vcel);
                TPZElementMatrix ek,ef;
                elvol->CalcStiff(ek, ef);
                if (iv==0) {
                    ekvol = ek;
                    efvol = ef;
                }
                else
                {
                    ekvol.fMat += ek.fMat;
                    efvol.fMat += ef.fMat;
                }
            }
            ekgrp.fMat.Print("EKGRP = ",std::cout,EMathematicaInput);
            ekvol.fMat.Print("EKVOL = ",std::cout,EMathematicaInput);
            break;
        }
    }

    
}

