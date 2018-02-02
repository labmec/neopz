#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "Common.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif

#ifdef _AUTODIFF
TElasticity2DAnalytic ElastExactLower;
TElasticity2DAnalytic ElastExactUpper;
#endif

void IntegrateDirect(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    bool scalarproblem = true;
    bool hasexact = false;

    int numrefskeleton = 4;
    int maxporder = 4;
    int counter = 1;
#ifdef _AUTODIFF
    ElastExact.fProblemType = TElasticity2DAnalytic::ESquareRoot;
#endif
    for ( int POrder = 1; POrder < maxporder; POrder += 1)
    {
        for (int irefskeleton = 0; irefskeleton < numrefskeleton; irefskeleton++)
        {
#ifdef _AUTODIFF
                ElastExact.fE = 10;
                ElastExact.fPoisson = 0.3;
                ElastExact.fPlaneStress = 0;
                ElastExact = ElastExact;
                ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRoot;
                
#endif
            bool elastic = !scalarproblem;
            TPZCompMesh *SBFem = SetupCrackedOneElement(irefskeleton, POrder, hasexact, elastic);
#ifdef _AUTODIFF
            ElastExact.fE = 10;
            ElastExact.fPoisson = 0.3;
            ElastExact.fPlaneStress = 0;
            ElastExactLower = ElastExact;
            ElastExactUpper = ElastExact;
            ElastExactLower.fProblemType = TElasticity2DAnalytic::ESquareRootLower;
            ElastExactUpper.fProblemType = TElasticity2DAnalytic::ESquareRootUpper;
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat1);
                mat->SetForcingFunction(ElastExactLower.ForcingFunction());
            }
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat2);
                mat->SetForcingFunction(ElastExact.ForcingFunction());
            }
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Emat3);
                mat->SetForcingFunction(ElastExactUpper.ForcingFunction());
            }
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc1);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(ElastExactLower.TensorFunction());
            }
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc2);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(ElastExact.TensorFunction());
            }
            if (hasexact)
            {
                TPZMaterial *mat = SBFem->FindMaterial(Ebc3);
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(mat);
                bc->SetType(0);
                mat->SetForcingFunction(ElastExactUpper.TensorFunction());
            }
#endif
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::ofstream gout("gmesh.vtk");
                TPZVTKGeoMesh::PrintGMeshVTK(SBFem->Reference(), gout,true);
                std::stringstream sout;
                SBFem->Print(sout);
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
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
            Analysis->SetExact(Elasticity_exact);
#endif
            //                ElasticAnalysis->SetExact(Singular_exact);
            
            TPZManVector<STATE> errors(3,0.);
            
            long neq = SBFem->Solution().Rows();
            
            if(!scalarproblem)
            {
                std::stringstream filename;
                filename << "SquareRootOneElement_NR_" << irefskeleton << "_P_" << POrder << ".vtk";
                TPZStack<std::string> vecnames,scalnames;
                // scalar
                vecnames.Push("Displacement");
                scalnames.Push("SigmaX");
                scalnames.Push("SigmaY");
                scalnames.Push("TauXY");
                scalnames.Push("EpsX");
                scalnames.Push("EpsY");
                scalnames.Push("EpsXY");
                Analysis->DefineGraphMesh(2, scalnames, vecnames, filename.str());
                Analysis->PostProcess(3);
            }

            if(0)
            {
                std::ofstream out("../CompMeshWithSol.txt");
                SBFem->Print(out);
            }
            
            if(hasexact)
            {
            
                std::cout << "Compute errors\n";
                
                Analysis->PostProcessError(errors);
                
    //                VerifyShapeFunctionIntegrity(Analysis->Mesh());
                
    //                IntegrateDirect(Analysis->Mesh());
                
                std::stringstream sout;
                sout << "../CrackRestrainedShape";
                if (scalarproblem) {
                    sout << "Scalar.txt";
                }
                else
                    sout << "Elastic2D.txt";
                
                std::ofstream results(sout.str(),std::ios::app);
                results.precision(15);
                results << "(*  numrefskel " << irefskeleton << " " << " POrder " << POrder << " neq " << neq << "*)" << std::endl;
                TPZFMatrix<double> errmat(1,3);
                for(int i=0;i<3;i++) errmat(0,i) = errors[i]*1.e6;
                std::stringstream varname;
                varname << "Errmat[[" << irefskeleton+1 << "]][[" << POrder << "]] = (1/1000000)*";
                errmat.Print(varname.str().c_str(),results,EMathematicaInput);
            }
            if(1)
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
                Analysis->ShowShape("OneElementCracked.vtk", eqindex);
            }
            
            delete Analysis;
            delete SBFem;
            //                exit(-1);
        }
        //            exit(-1);
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

